library(data.table)
library(LDlinkR)
library(parallel)
library(XML)
'%ni%' <- Negate('%in%')

# Default arguments
args <- list()
args$scratch_folder <- 'scratch'
args$ldlink_winsize <- 100000
args$ldlink_genome <- 'grch38_high_coverage'
args$ldlink_pop <- 'no --ldlink_pop argument provided!'
args$ldlink_r2 <- 0.8
args$threads <- 1
args$output_fnm <- 'my_results.txt'
# TODO Gotta be a cleaner way to do this no?
args$bcftools_exe <- if(file.exists('tool/bcftools/bcftools'))      'tool/bcftools/bcftools'      else if(file.exists('../tool/bcftools/bcftools'))      '../tool/bcftools/bcftools'      else Sys.which('bcftools')
args$bgenix_exe   <- if(file.exists('tool/bgen/build/apps/bgenix')) 'tool/bgen/build/apps/bgenix' else if(file.exists('../tool/bgen/build/apps/bgenix')) '../tool/bgen/build/apps/bgenix' else Sys.which('bgenix')
args$plink2_exe   <- if(file.exists('tool/plink2'))                 'tool/plink2'                 else if(file.exists('../tool/plink2'))                 '../tool/plink2'                 else Sys.which('plink2')

# Parse command-line args
pieces <- tstrsplit(split='--', paste(gsub('\n','',commandArgs(T)), collapse=' '))[-1]
for(piece in pieces) { parts <- unlist(tstrsplit(piece,' |=')); args[[parts[[1]]]] <- parts[-1] }
args <- lapply(args, \(x) x[x!=''])

# --- Input validation ---
recognized_args <- c('geno_files','sample_file','score_file','fill_vcf_ids_with','ldlink_token','ldlink_pop','ldlink_r2','ldlink_genome','ldlink_winsize','scratch_folder','threads','output_fnm','bcftools_exe','bgenix_exe','plink2_exe', paste0('score_file_',c('chr','pos','id','ref','alt','ea'),'_col'), 'score_file_weight_cols')
if(any(names(args) %ni% recognized_args)) stop('Unrecognized argument(s):', paste0(' --',setdiff(names(args),recognized_args)))

## Required args provided? Not too many? They exist?
if(is.null(args$geno_files            )) stop('--geno_files is required!\nThis may be .vcf[.gz]/.bcf[.gz], .bgen, .bed+.bim+.fam, or .pgen+.pvar+.psam file(s).')
if(is.null(args$score_file            )) stop('--score_file is required!')
if(is.null(args$score_file_chr_col    )) stop('--score_file_chr_col is required! This refers to the name or index of the chromosome column in your score file.\nThe chromosome format (e.g. "1" vs. "chr1") should match that in your --geno_files input.')
if(is.null(args$score_file_pos_col    )) stop('--score_file_pos_col is required! This refers to the name or index of the position column in your score file.')
if(is.null(args$score_file_id_col     )) stop('--score_file_id_col  is required! This refers to the name or index of the ID column in your score file.\nThe ID format (e.g. "rs1234" vs. "1:1234" vs. "chr1:1234:A:G" etc.) should match that of your --geno_files input.')
if(is.null(args$score_file_ref_col    )) stop('--score_file_ref_col is required! This refers to the name or index of the reference allele column in your score file.')
if(is.null(args$score_file_alt_col    )) stop('--score_file_alt_col is required! This refers to the name or index of the alternate allele column in your score file.')
if(is.null(args$score_file_ea_col     )) stop('--score_file_ea_col  is required! This refers to the name or index of the effect allele column in your score file.')
if(is.null(args$score_file_weight_cols)) stop('--score_file_weight_cols  is required! This refers to the name(s) or index(es) of the weights column(s) in your score file.')

## Resolve patterns to actual files (e.g. *.bgen -> chr1.bgen, chr2.bgen)
### Resolve pattern over URL
if( all(grepl('^https://', args$geno_files)) & any(grepl('\\*', args$geno_files)) ) stop('Using a pattern with https:// URLs is not supported. Manually specifying all the URLS without using a pattern WILL still work. But first try editing your URL pattern to use http:// instead of https://, that usually still works.')
if(!all(grepl(      '://', args$geno_files)) & any(grepl('://', args$geno_files)) ) stop('Cannot use both local genotype files and genotype files over URL') # Technically could but complicates code a bit
if( all(grepl(      '://', args$geno_files)) & any(grepl('\\*', args$geno_files)) ) {
  message('Detected a pattern of URLs. Scraping the site for links to the files...')
  args$geno_files <-
    htmlParse(dirname(args$geno_files)) |> # Scrape HTML
    xpathSApply('//a/@href') |> # Get only hrefs
    grep(pattern=basename(args$geno_files),value=T) |> # Only keep matching user's pattern (probably *.vcf.gz)
    (\(link) file.path(dirname(args$geno_files),link))() # Construct full URL+filename path
  args$geno_files <- args$geno_files[!grepl('tbi$|csi$',args$geno_files)] # Don't include index files
  if(!all(grepl('bcf(.gz)?$|vcf(.gz)?$',args$geno_files))) stop('Only vcf[.gz]/bcf[.gz] files can be accessed over URL.')
}
### Resolve patterns for local files
                                         args$score_file  <- list.files(dirname(args$score_file),  pattern=paste(collapse='|',basename(args$score_file )), full.names=T)
if(!any(grepl('://',args$geno_files))) { args$geno_files  <- list.files(dirname(args$geno_files),  pattern=paste(collapse='|',basename(args$geno_files )), full.names=T) }
if(!is.null(args$sample_file))         { args$sample_file <- list.files(dirname(args$sample_file), pattern=paste(collapse='|',basename(args$sample_file)), full.names=T) }

arg_lengths <- lengths(args)[names(args) %ni% c('geno_files','score_file_weight_cols')]
args_too_long <- names(arg_lengths)[arg_lengths>1]
if(length(args_too_long)>0) stop('Only one parameter is allowed for each argument, except for --geno_files and --score_file_weight_cols\n', paste(args_too_long,collapse=' and '), ' have too many arguments!')
if(length(args$geno_files)==0) stop('Provided --geno_files do not exist!  \n(If you\'re submitting a job on a compute cluster, you may need to specify absolute file paths.)')
if(length(args$score_file)==0) stop('Provided --score_file does not exist!\n(If you\'re submitting a job on a compute cluster, you may need to specify absolute file paths.)')

## Detect --geno_file type (and --sample_file)
if(all(grepl('bgen$|bgi$',args$geno_files))) {
  if(is.null(args$sample_file)   ) stop('--sample_file is required since your --geno_files are .bgen format!')
  if( length(args$sample_file)==0) stop('Provided --sample_file does not exist!')
  # TODO would be good to detect index files and warn if not found
  args$geno_files <- args$geno_files[!grepl('bgi$',args$geno_files)] # Don't include index files
  geno_files_type <-'bgen' 
} else if(all(grepl('bcf(.gz)?$|vcf(.gz)?$|tbi$|csi$',args$geno_files))) {
  # TODO would be good to detect index files and warn if not found
  if(all(grepl('bcf(.gz)?$|tbi$|csi$',args$geno_files))) geno_files_type <- 'bcf'
  if(all(grepl('vcf(.gz)?$|tbi$|csi$',args$geno_files))) geno_files_type <- 'vcf'
  args$geno_files <- args$geno_files[!grepl('tbi$|csi$',args$geno_files)] # Don't include index files
} else if(all(grepl('bed$|bim$|fam$',args$geno_files))) {
  args$geno_files <- unique(sub('.bed$|.bim$|.fam$','',args$geno_files)) # PLINK only takes the file prefix
  geno_files_type <- 'plink1'
} else if(all(grepl('pgen$|pvar$|psam$|bim$|fam$',args$geno_files))) {
  args$geno_files <- unique(sub('pgen$|pvar$|psam$|bim$|fam$','',args$geno_files)) # PLINK only takes the file prefix
  geno_files_type <- 'plink2'
} else if(all(grepl('gds$',args$geno_files))) {
  geno_files_type <- 'gds'
} else { stop(paste(c('Input geno_files format is unsupported, or is a mix of different formats. Below is the list of detected geno_files:',args$geno_files),collapse='\n')) }
message('Detected ',geno_files_type,' files:\n    ', paste(args$geno_files,collapse='\n    '), '\n')

if(geno_files_type %in% c('bcf','vcf') && !file.exists(args$bcftools_exe)) stop('Wasn\'t able to find bcftools executable, required because you input bcf/vcf --geno_files.\n Please provide a path to --bcftools_exe, or, convert the files to another format using PLINK and use those instead.')
if(geno_files_type  ==       'bgen'    && !file.exists(args$bgenix_exe)  ) stop('Wasn\'t able to find bgenix   executable, required because you input bgen --geno_files.\n Please provide a path to --bgenix_exe, or, convert the files to another format using PLINK and use those instead.')
if(                                       !file.exists(args$plink2_exe)  ) stop('Wasn\'t able to find plink2   executable. Please provide a path to --plink2_exe.')
if(geno_files_type  ==       'gds'     && !require(SeqArray,quietly=T)   ) {
  warning('Wasn\'t able to find SeqArray R package, required because you input gds --geno_files. Installing now.\n
           If you plan to run this pipeline repeatedly on the cloud or on a compute cluster, consider using an R environment with SeqArray pre-installed so that you don\'t waste time installing the package repeatedly.')
  if(!require("BiocManager",quietly=T)) install.packages("BiocManager"); BiocManager::install("SeqArray",ask=F)
}

## Make sure --score_file's format is correct.
score_dt <- fread(args$score_file)
score_col_specs <- list(chr=args$score_file_chr_col,
                        pos=args$score_file_pos_col,
                        id =args$score_file_id_col,
                        ref=args$score_file_ref_col,
                        alt=args$score_file_alt_col,
                        ea =args$score_file_ea_col,
                        weights=args$score_file_weight_cols) |> unlist()

score_col_idxs <- sapply( score_col_specs, \(col) if(grepl('^[0-9]+$',col)) as.numeric(col) else which(names(score_dt) == col) )
# TODO error if the same column is used twice. Separate ref/alt/ea cols should be enforced.

score_col_farthest <- max(unlist(score_col_idxs))
if(score_col_farthest  > ncol(score_dt)) stop('There are fewer than ',score_col_farthest,' columns in the score file!')

score_col_specs_not_found <- score_col_specs[lengths(score_col_idxs)==0]
if(length(score_col_specs_not_found)>0) stop('Some columns were not found in the score file! These columns were: ', paste(collapse=' ',score_col_specs_not_found))

score_og_weight_nms <- names(score_dt)[score_col_idxs[grepl('weights',names(score_col_idxs))]]
message(paste0('Using score file\'s "',names(score_dt)[score_col_idxs],'" column as ',names(score_col_specs),'\n'), '\n\x1b[36mMake sure the correct score file columns were selected!!\x1b[m\n')

score_dt <- fread(args$score_file) |> setnames(old=score_col_idxs, new=names(score_col_specs))

if(score_dt[, !is.integer(chr) & !is.character(chr)]) message('Warning: score file\'s chromosome column is not character or integer type which is suspicious. Here is a sample: ',        paste(head(score_dt$chr),collapse=' '))
if(score_dt[, !is.integer(pos)                     ]) message('Warning: score file\'s position column is not integer type which is suspicious. Here is a sample: ',                       paste(head(score_dt$pos),collapse=' '))
if(score_dt[, any(grepl('[^ATCGNatcgn]', ref))     ]) message('Warning: score file\'s reference allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ', paste(head(score_dt$ref),collapse=' '))
if(score_dt[, any(grepl('[^ATCGNatcgn]', alt))     ]) message('Warning: score file\'s alternate allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ', paste(head(score_dt$alt),collapse=' '))
if(score_dt[, any(grepl('[^ATCGNatcgn]', ea ))     ]) message('Warning: score file\'s effect allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ',    paste(head(score_dt$ea ),collapse=' '))
if(score_dt[,!any(grepl('[0-9]',         id ))     ]) message('Warning: score file\'s id column has no numbers in it, which is suspicious, here is a sample: ',                           paste(head(score_dt$id ),collapse=' ')) 
if(score_dt[,              anyDuplicated(id)>0     ]) stop('Found duplicated IDs in the score file! Please merge the duplicates.') # TODO suggest data.table by= or tibble group_by code maybe
if(score_dt[,!all(sapply(.SD,is.numeric)), .SDcols=patterns('weights')]) stop('Score file weights columns do not contain all-numeric values! Here is a sample:\n', paste(collapse='\n',capture.output(score_dt[,.SD,.SDcols=patterns('weights')])))

## LDlink inputs
if(is.null(args$ldlink_token)) {
  message('\x1b[31mNo LDlink API token provided\x1b[m, so no proxies will be gotten for missing variants.')
} else {
  invalid_pop_codes <- args$ldlink_pop[args$ldlink_pop %ni% list_pop()$pop_code]
  if(length(invalid_pop_codes)>0) {
    message('Error: Invalid LDlink population code(s): ', paste(invalid_pop_codes,collapse=' '),'\n',
         'Available populations listed below. You can use more than one -- but LD calculation will be slower.\n',
         paste(capture.output(list_pop()), collapse = '\n'))
    stop('^^')
  }
}
if(args$ldlink_genome %ni% c('grch37','grch38','grch38_high_coverage')) stop('--ldlinkg_genome must be one of: grch37 grch38 grch38_high_coverage')

# --- Finished input validation ---

dir.create(args$scratch_folder, recursive=T, showWarnings=F)

# Get list of score file variant IDs present in the geno_files (not extracting genotype data yet). Used for getting proxies later.
# Why not use PLINK to read all formats? Because it cannot make use of .bgi (bgen) or .csi (vcf/bcf) index files. Instead it would try to convert everything to plink2 format first, (very slow for huge datasets).
filter_ranges_fnm <- file.path(args$scratch_folder,'filter_ranges.txt')
filter_ids_fnm    <- file.path(args$scratch_folder,'filter_variant_ids.txt'   )
writeLines(score_dt[,paste0(chr,'\t',pos)], filter_ranges_fnm)
writeLines(score_dt$id, filter_ids_fnm)

bcftools_annotation_part <- {
  if(is.null(args$fill_vcf_ids_with))
    ''
  else if(args$fill_vcf_ids_with=='chr:pos:ref:alt')
    paste(args$bcftools_exe, 'view -Ov | awk \'BEGIN {OFS="\t"} /^#/ {print; next} {$3=$1":"$2":"$4":"$5; print}\' | ')
  else if(grepl("bcf(.gz)?$|vcf(.gz)?$", args$fill_vcf_ids_with)) # Another file to get IDs from
    paste(args$bcftools_exe, 'annotate -c ID -a', args$fill_vcf_ids_with, '-Ou | ')
  else
    stop('Argument to --fill_vcf_ids must be \'chr:pos:ref:alt\' or a BCF/VCF file to get IDs from.')
}

id_extraction <- \(geno_file, filter_ranges_file, filter_ids_file, output_fnm) {
  cmd <-
    if(geno_files_type=='bgen') paste0(
      args$bgenix_exe, ' -list -g ', geno_file, # -list Only list variant info, no genotype data
      #' -incl-range ', filter_ranges_file, # Purposefully commented out. bgenix extracts the _union_ of -incl-range and -incl-rsids. This can lead to variants we don't want which are at the same position as the desired variants. Doesn't help speed because bgen files are already indexed by ID, not only chr+pos.
      ' -incl-rsids ', filter_ids_file,
      ' | cut -f1', # Keep only ID column
      ' | sed \'/# bgenix/d\' | sed \'/alternate_ids/d\'', # Remove comment and colname lines
      ' > ', output_fnm
    )
    else if(geno_files_type %in% c('vcf','bcf')) paste0(
      args$bcftools_exe, ' view -G ', geno_file, # -G means only list variant info, no genotype data nor header
      ' -R ', filter_ranges_file, # Need to filter by ranges otherwise will be very slow, because VCF/BCF files are only indexed by chr+pos, not ID. Unlike bgenix, -R combined with a -i filter will take the intersection of variants.
      ' -Ou | ',
      bcftools_annotation_part, # If no --fill_vcf_ids_with given, this is an empty string
      args$bcftools_exe, ' query',
      ' -i"ID=@',filter_ids_file,'"',
      ' -f"%ID"', # Output only IDs
      ' > ', output_fnm
    )
    else if(geno_files_type=='plink1') paste0(
      args$plink2_exe, ' --bfile ', geno_file,
      ' --rm-dup force-first', # If dups, PLINK errors unless this is here.
      ' --extract ', filter_ids_file,
      ' --write-snplist --out ', sub('.snplist$','',output_fnm) # PLINK automatically adds '.snplist' to the name
    )
    else if(geno_files_type=='plink2') paste0(
      args$plink2_exe, ' --pfile ', geno_file,
      ' --rm-dup force-first',
      ' --extract ', filter_ids_file,
      ' --write-snplist --out ', sub('.snplist$','',output_fnm)
    )
    else if(geno_files_type=='gds') {
      gds <- seqOpen(geno_file)
      seqSetFilterPos(gds, score_dt$chr, score_dt$pos, score_dt$ref, score_dt$alt)
      ids <- seqGetData(gds, '$chrom_pos_allele')
      ids <- gsub('_',':',ids) # chr:pos_ref_alt -> chr:pos:ref:alt
      writeLines(ids,output_fnm)
      seqClose(gds)
      ' ' # Dummy cmd
    }
  system(cmd, ignore.stderr=T)
}
found_id_fnms <- paste0(args$scratch_folder,'/',basename(args$score_file),'-',basename(args$geno_files),'.snplist')

message('Checking for variants in the score_file missing from the geno_files...')
## Don't rerun this if the found id files alredy exist, bottleneck.
if(!all(file.exists(found_id_fnms))) invisible(mcmapply(args$geno_files, filter_ranges_fnm, filter_ids_fnm, found_id_fnms, FUN=id_extraction, mc.cores=args$threads))

ids_found <- sapply(found_id_fnms, scan, what=character(), quiet=T) |> unlist()
ids_not_found <- setdiff(score_dt$id, ids_found)
ids_not_found_chrpos <- score_dt[id %in% ids_not_found, paste0('chr',chr,':',pos)] |> sub(pattern='chrchr',replacement='chr') # LDproxy only takes rsIDs or chr:pos IDs

if(length(ids_found)==0) { unlink(found_id_fnms); stop('No IDs from score_file were found in geno_files! This is probably a chromosome or ID format mismatch issue, e.g. "1" vs. "chr1", or rsID vs. chr:pos ID, etc.. Or, maybe your geno_files are VCFs without IDs (like 1000G files), and you forgot to specify --fill_vcf_ids_with?') } # TODO print out a couple Ids from each file.
message(length(ids_not_found),'/',nrow(score_dt),' variant IDs in score file not found in genotype file(s).')

# Find proxies for score file variants mising from the genotype data
if(!is.null(args$ldlink_token) & length(ids_not_found)>0) {

  old_wd <- setwd(args$scratch_folder) # LDproxy can only output in working directory, but we want it in the scratch dir
  ## Don't want to ask the busy LDproxy server the same thing twice, so give the output a unique filename based on the inputs and check if it already exists before calling LDproxy.
  proxy_output_fnm <- paste0('proxy-',basename(args$score_file),'-',args$ldlink_pop,'-',digest::digest(c(ids_not_found,args$ldlink_pop)),'.txt')
  if(!file.exists(proxy_output_fnm)) {
    LDproxy_batch(ids_not_found_chrpos, pop=args$ldlink_pop, token=args$ldlink_token, genome_build=args$ldlink_genome, append=T, win_size=args$ldlink_winsize)
    file.rename(paste0('combined_query_snp_list_',args$ldlink_genome,'.txt'), proxy_output_fnm)
  }
  proxy_dt <- suppressWarnings(fread(proxy_output_fnm))
  setwd(old_wd)
  
  proxy_dt <- proxy_dt[
      R2 > args$ldlink_r2
    ][!grepl('-',Alleles) # TODO ignoring these 1-bp indels b/c idk how to handle them
    ][Distance != 0
  
    # Rm unneeded cols
    ][, `:=`(V1=NULL, MAF=NULL, Dprime=NULL, FORGEdb=NULL, RegulomeDB=NULL, Function=NULL)
  
    ][, chr := sub(':.*','',Coord)
    ][, pos := sub('.*:','',Coord)
    ][, ref := sub('\\(','',tstrsplit(Alleles,'/')[[1]])
    ][, alt := sub('\\)','',tstrsplit(Alleles,'/')[[2]])
    ][, chr_n := sub('^chr','',chr)
  
    # Figure out which allele of the proxy corresponds to the effect allele of the original variant
    ][, id_from_score_file := sapply(query_snp, \(query_chrpos) ids_not_found[which(ids_not_found_chrpos==query_chrpos)])
    ][, ea_from_score_file := sapply(id_from_score_file, \(id2) score_dt[id==id2, ea])
    ][, ea := mapply(Correlated_Alleles, ea_from_score_file, FUN = \(as, s_ea) {
                as <- strsplit(as, ',|=')[[1]] # 'A=C,T=G' -> list(A,C,T,G). A&T would be query variant's alleles, C&G the proxy's.
                if(as[1]==s_ea) as[2] else as[4]
              })
    # TODO ^ ea code mega ugly surely more elgant way
  
    # Instead of detecting the format of the user input IDs in the score file, just throw every possible ID format at the wall, simpler.
    ][,       rs_id  := RS_Number
    ][,     `c:p_id` := paste0(chr,  ':',pos)
    ][, `c:p:r:a_id` := paste0(chr,  ':',pos,':',ref,':',alt)
    ][, `c:p_r_a_id` := paste0(chr,  ':',pos,'_',ref,'_',alt)
    ][,     `n:p_id` := paste0(chr_n,':',pos)
    ][, `n:p:r:a_id` := paste0(chr_n,':',pos,':',ref,':',alt)
    ][, `n:p_r_a_id` := paste0(chr_n,':',pos,'_',ref,'_',alt)
  ]

  writeLines(unlist(proxy_dt[,.SD,.SDcols=patterns('_id')]), filter_ids_fnm)
  writeLines(       proxy_dt[,paste0(chr,'\t',pos)],         filter_ranges_fnm)
  
  proxy_found_id_fnms <- paste0(args$scratch_folder,'/',basename(proxy_output_fnm),'-',basename(args$geno_files),'.snplist')
  
  if(!all(file.exists(proxy_found_id_fnms))) invisible(mcmapply(args$geno_files, filter_ranges_fnm, filter_ids_fnm, proxy_found_id_fnms, FUN=id_extraction, mc.cores=args$threads))
  
  proxy_ids_found <- sapply(proxy_found_id_fnms, scan, what=character(), quiet=T) |> unlist()
  
  if(length(proxy_ids_found)==0) { warning('Unfortunately \x1b[36mnone of the possible proxies were found in geno_files!\x1b[m Continuing without them.')
  } else { # Update score_dt with proxies
    where_matched_id <- Reduce('|', lapply(proxy_ids_found,'==',proxy_dt))
    
    proxy_dt$replacement_id <- ''
    invisible(mapply( # Will use the IDs that matched in the geno file to replace IDs in the score file of the variants that needed proxies.
      row(proxy_dt)[where_matched_id],
      col(proxy_dt)[where_matched_id],
      FUN= \(r,c) set(proxy_dt, r, 'replacement_id', proxy_dt[r,c,with=F])
    ))
    
    proxy_dt <- proxy_dt[rowSums(where_matched_id)>0] # Only proxies found in the geno file
    proxy_dt <- proxy_dt[proxy_dt[, .I[which.max(R2)], by=id_from_score_file]$V1] # Only the best-correlated proxy for each original variant
    
    ids_no_proxy <- setdiff(ids_not_found, proxy_dt$id_from_score_file)
    if(length(ids_no_proxy)>0) message('\x1b[36mCould not find proxies for ',length(ids_no_proxy),'/',length(ids_not_found),' of the variants which needed proxies\x1b[m:\n',paste(ids_no_proxy,collapse='\n'),'\nTo make sure this variant isn\'t omitted, you could try decreasing --ldlink_r2 or increasing --ldlink_winsize.\nFor now, continuing without these proxyless variants.\n')
    
    score_dt <- proxy_dt[score_dt, on=c(id_from_score_file='id')]
    score_dt <- score_dt[
      ][                      , id := id_from_score_file
      ][!is.na(replacement_id), id := replacement_id
      ][is.na(ea), ea := i.ea
      ][is.na(pos), pos := i.pos
      ][, chr := i.chr # TODO bug if proxy on different chr, but that never happens right?
    ]
    #score_dt[id!=id_from_score_file] # TMP: to ensure proxies look right
  }

} # END if(!is.null(args$ldlink_token))

score_dt_simple <- score_dt[, cbind(chr,pos,id,ea,.SD), .SDcols=patterns('weights[0-9]+')]
setnames(score_dt_simple, grep('weights[0-9]+',names(score_dt_simple), value=T), score_og_weight_nms) # TODO very messy way of restoring weight column nms
score_dt_simple_fnm <- file.path(args$scratch_folder,'score_file-formatted_for_plink.csv') # TODO better run-specific nm
fwrite(score_dt_simple, score_dt_simple_fnm, sep=' ')

# Extract full genotype data for the variants in the --score_file.  
geno_subset_file_paths <- paste0(args$scratch_folder,'/',basename(args$score_file),'-',basename(args$geno_files))

writeLines(score_dt_simple[,paste0(chr,'\t',pos)], filter_ranges_fnm)
writeLines(score_dt_simple$id,                     filter_ids_fnm)

geno_extraction <- \(geno_file, geno_subset_file) {
  cmd <-
    if(geno_files_type=='bgen') paste0(
      args$bgenix_exe, ' -g ', geno_file,
      #' -incl-range ', filter_ranges_file, # Purposefully commented out. bgenix extracts the _union_ of -incl-range and -incl-rsids. This can lead to variants we don't want which are at the same position as the desired variants. Doesn't help speed because bgen files are already indexed by ID, not only chr+pos. 
      ' -incl-rsids ', filter_ids_fnm,
      ' > ', geno_subset_file
    )
    else if(geno_files_type %in% c('vcf','bcf')) paste0(
      args$bcftools_exe, ' view ', geno_file,
      ' -R ', filter_ranges_fnm, # Need to filter by ranges otherwise will be very slow, because VCF/BCF files are only indexed by chr+pos, not ID. Unlike bgenix, -R combined with a -i filter will take the intersection of variants.
      ' -Ob | ', # Uncompressed BCF while piping, for speed
      bcftools_annotation_part, # If no --fill_vcf_ids_with given, this is just ''
      args$bcftools_exe, ' view ',
      ' -i"ID=@',filter_ids_fnm,'"',
      ' -Ob', # Output compressed BCF format
      ' > ', geno_subset_file
    )
    else if(geno_files_type=='plink1') paste0(
      args$plink2, ' --bfile', geno_file,
      ' --extract ', filter_ids_fnm,
      ' --make-pgen --out ', geno_subset_file # Might as well conver to PLINK2 .pgen now
    )
    else if(geno_files_type=='plink2') paste0(
      args4plink2, ' --pfile', geno_file,
      ' --extract ', filter_ids_fnm,
      ' --make-pgen --out ', geno_subset_file
    )
    else if(geno_files_type=='gds') {
      gds <- seqOpen(geno_file)
      seqSetFilterPos(gds, score_dt$chr, score_dt$pos, score_dt$ref, score_dt$alt, verbose=F)
      seqGDS2BED(gds, geno_subset_file, write.rsid='chr_pos_ref_alt', verbose=F) # PLINK2 cannot handle GDS files, so have to convert
      system(paste0('awk \'{gsub(/_/,":"); print}\' ',geno_subset_file,'.bim > tmp && mv tmp ',geno_subset_file,'.bim')) # Why not sed -i 's/_/:/g'? sed -i syntax is different between Linux & Mac/BSD. TODO: figure out something better, this sucks
      seqClose(gds)
      ' '
    }
  system(cmd, ignore.stderr=T)
}
 
message('Extracting genotype data...')
# Don't rerun this if the subsetted files alredy exist, big bottleneck.
  # TODO Unfortunately it's trickier to check this for GDS files since we chang ethe file type. Maybe move the file.exists check INSIDE the function then?
if(!all(file.exists(geno_subset_file_paths))) invisible(mcmapply(args$geno_files, geno_subset_file_paths, FUN = geno_extraction, mc.cores=args$threads))

if(geno_files_type=='gds') geno_files_type <- 'plink1' # As mentioned before, PLINK2 cannot handle GDS files
unlink(paste0(basename(args$geno_files),'.tbi')) # Clean up index files that bcftools pulls in
                                                                                                        # v TODO v: is there any way to be sure which is the reference alle in a user-provided .bgen?
if(geno_files_type  ==       'bgen'   ) plink_input_flags <- paste0('--bgen  ',geno_subset_file_paths,' ref-unknown --sample ',args$sample_file)
if(geno_files_type %in% c('bcf','vcf')) plink_input_flags <- paste0('--bcf   ',geno_subset_file_paths) # No matter what, the subset of genotype data we extract will be in BCF format.
if(geno_files_type  ==      'plink1'  ) plink_input_flags <- paste0('--bfile ',geno_subset_file_paths)
if(geno_files_type  ==      'plink2'  ) plink_input_flags <- paste0('--pfile ',geno_subset_file_paths)
stopifnot(!is.null(plink_input_flags))
plink_output_flags <- paste0('--allow-misleading-out-arg --out ', geno_subset_file_paths)
# PRS result files will be '<scratch_folder>/<score_filename>-<geno_filename>.sscore'

# Finally calculate PRSes!!
plink_prs_cmds <- paste(
  args$plink2_exe, plink_input_flags,
  '--score', score_dt_simple_fnm, '3 4',
    'ignore-dup-ids',
    'header-read', # Use row as names for clusters
    'cols=scoresums', # Output plain sum(dosages*weights), without averaging, so that we can sum scores across chromosomes. Then we can take the average.
    'list-variants',
  '--score-col-nums', paste0('5-',ncol(score_dt_simple)), # Skip ID & effect allele columns
  #'--rm-dup force-first',
  plink_output_flags
)
#message('\nRunning the following commands to calculate PRSes...:\n    ', paste(plink_prs_cmds,collapse='\n    '))

message('Calculating PRSes...')
err_codes <- mclapply(plink_prs_cmds,system, ignore.stdout=T,ignore.stderr=T, mc.cores=args$threads)

if(any(err_codes!=0)) message('\nPLINK PRS calculation errors happened for some files. Check these logs:\n    ', paste0(geno_subset_file_paths[err_codes!=0],'.log',collapse='\n    '), '\n(A common "error" is that no variants were in the input file, but that may just be because you provided files for all chromosomes but your score file didn\'t have a variant in every chromosome. If so, this is not a concern.)')
# TODO: No variants found should not be considered an error, because it is normal that a user provides separate chr 1-22 files, but there is not a weighted variant in every chr. To avoid this, right after the LDlink step, geno_subset_files with no score file variants should be removed. Something like that.


# Sum scores from each run
prs_files      <- list.files(path=args$scratch_folder, pattern=paste0(basename(geno_subset_file_paths),'.sscore$',collapse='|')) # Why list.files and not just add .sscore? Because some .sscore files might not exist if a chr had 0 of the score file's variants in it.
prs_file_paths <- file.path(args$scratch_folder, prs_files)
prs_sums <- rbindlist(lapply(prs_file_paths,fread))[, lapply(.SD, sum), by='#IID']

fwrite(prs_sums, args$output_fnm)
message('Done! PRS results file: ', args$output_fnm)

# TODO: --q-score-range option useful for global PRS, see Kenny script for referene
