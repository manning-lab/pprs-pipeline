deps <- c('data.table','LDlinkR','parallel','XML') 
install.packages(deps[lengths(sapply(deps,find.package,quiet=T))==0], repos='https://cloud.r-project.org')

library(data.table)
library(LDlinkR)
library(parallel)
library(XML)

'%ni%' <- Negate('%in%')

get_bcftools <- \() {
  message('Installing bcftools to ./bcftools/bcftools...')
  system('git clone --depth 1 --recurse-submodules https://github.com/samtools/htslib.git')
  system('git clone --depth 1 https://github.com/samtools/bcftools.git')
  err_code <- system('cd bcftools && make -j')
  if(err_code!=0) stop('Failed to build bcftools')
}
get_bgenix <- \() {
  message('Installing bgenix to ./bgen/build/apps/bgenix... (this takes a long time)')
  if(file.exists('bgen.tgz')) unlink('bgen.tgz') # In case user got impatient and leaves a half-untarred directory
  system('curl -LO http://code.enkre.net/bgen/tarball/release/bgen.tgz')
  system('tar zxf bgen.tgz')
  system('mv bgen.tgz bgen') # Extracted directory retains .tgz in its name for some reason
  system('cd bgen && ./waf configure && ./waf')
}
get_plink2 <- \() {
  message('Installing plink2 to ./plink2...')
  # FIXME: Windows and Mac don't actually have lscpu so this will break there. In the meantime Mac users can provide --thing_exe.
  win   <- Sys.info()['sysname']=='Windows'
  mac   <- Sys.info()['sysname']=='Darwin'
  linux <- Sys.info()['sysname']=='Linux'
  arm   <- Sys.info()['architecture']=='arm64'
  intel <- system('lscpu | grep GenuineIntel', ignore.stdout=T)==0
  amd   <- system('lscpu | grep AuthenticAMD', ignore.stdout=T)==0
  avx2  <- system('lscpu | grep         avx2', ignore.stdout=T)==0

  url <- if(win   &         avx2) 'https://s3.amazonaws.com/plink2-assets/alpha6/plink2_win_avx2_20241222.zip'
    else if(win                 ) 'https://s3.amazonaws.com/plink2-assets/alpha6/plink2_win64_20241222.zip'
    else if(mac   &         arm ) 'https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_arm64_20241222.zip'
    else if(mac   &         avx2) 'https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_avx2_20241222.zip'
    else if(mac                 ) 'https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_20241222.zip'
    else if(linux & intel & avx2) 'https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20241222.zip'
    else if(linux &  amd  & avx2) 'https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_amd_avx2_20241222.zip'
    else if(linux               ) 'https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241222.zip'
    else stop('PLINK2 installation failed; machine not supported.')

  system(paste('curl -L', url, '-o plink2.zip && unzip -o plink2.zip && rm plink2.zip'))
}
get_seqarray <- \() { if(!require("BiocManager",quietly=T)) install.packages("BiocManager"); BiocManager::install("SeqArray",ask=F) }


# Default arguments
args <- list()
args$scratch_folder <- 'scratch'
args$ldlink_winsize <- 100000
args$ldlink_genome <- 'grch38_high_coverage'
args$ldlink_pop <- 'no --ldlink_pop argument provided!'
args$ldlink_r2 <- 0.8
args$threads <- 1
args$output_fnm <- 'my_results.txt'
args$bcftools_exe <- if(Sys.which('bcftools')=='') './bcftools/bcftools'      else Sys.which('bcftools')
args$bgenix_exe   <- if(Sys.which( 'bgenix' )=='') './bgen/build/apps/bgenix' else Sys.which( 'bgenix' )
args$plink2_exe   <- if(Sys.which( 'plink2' )=='') './plink2'                 else Sys.which( 'plink2' )

# Parse command-line args
pieces <- tstrsplit(split='--', paste(gsub('\n','',commandArgs(T)), collapse=' '))[-1]
for(piece in pieces) { parts <- unlist(tstrsplit(piece,' |=')); args[[parts[[1]]]] <- parts[-1] }
args <- lapply(args, \(x) x[x!=''])

# --- Input validation ---
recognized_args <- c('geno_files','sample_file','score_file','ldlink_token','ldlink_pop','ldlink_r2','ldlink_genome','ldlink_winsize','ldlink_yes_really','scratch_folder','threads','output_fnm','bcftools_exe','bgenix_exe','plink2_exe', paste0('score_file_',c('chr','pos','ref','alt','ea'),'_col'), 'score_file_weight_cols')
if(any(names(args) %ni% recognized_args)) stop('Unrecognized argument(s):', paste0(' --',setdiff(names(args),recognized_args)))

## Required args provided? Not too many? They exist?
if(is.null(args$geno_files            )) stop('--geno_files is required!\nThis may be .vcf[.gz]/.bcf[.gz], .bgen, .bed+.bim+.fam, or .pgen+.pvar+.psam file(s).')
if(is.null(args$score_file            )) stop('--score_file is required!')
if(is.null(args$score_file_chr_col    )) stop('--score_file_chr_col is required! This refers to the name or index of the chromosome column in your score file.\nThe chromosome format (e.g. "1" vs. "chr1") should match that in your --geno_files input.')
if(is.null(args$score_file_pos_col    )) stop('--score_file_pos_col is required! This refers to the name or index of the position column in your score file.')
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
if(!any(grepl('://',args$geno_files))) { args$geno_files  <- list.files(dirname(args$geno_files),  pattern=paste(collapse='|',basename(args$geno_files )), full.names=T) |> unique() }
if(!is.null(args$sample_file))         { args$sample_file <- list.files(dirname(args$sample_file), pattern=paste(collapse='|',basename(args$sample_file)), full.names=T)             }

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

if(geno_files_type  ==       'gds'     && !require(SeqArray,quietly=T)   ) { message('Wasn\'t able to find SeqArray R package.  Installing now. (Required because you input   gds   --geno_files.)'); get_seqarray() }
if(geno_files_type %in% c('bcf','vcf') && !file.exists(args$bcftools_exe)) { message('Wasn\'t able to find bcftools executable. Installing now. (Required because you input bcf/vcf --geno_files.)'); get_bcftools() }
if(geno_files_type  ==       'bgen'    && !file.exists(args$bgenix_exe)  ) { message('Wasn\'t able to find  bgenix  executable. Installing now. (Required because you input  bgen   --geno_files.)'); get_bgenix()   }
if(                                       !file.exists(args$plink2_exe)  ) { message('Wasn\'t able to find  plink2  executable. Installing now.');                                                    get_plink2()   }

## Make sure --score_file's format is correct.
score_dt <- fread(args$score_file)
score_colnms <- c(chr=args$score_file_chr_col,
                  pos=args$score_file_pos_col,
                  ref=args$score_file_ref_col,
                  alt=args$score_file_alt_col,
                  ea =args$score_file_ea_col,
                  weights=args$score_file_weight_cols)

score_colnms_not_found <- score_colnms[score_colnms %ni% names(score_dt)]
if(length(score_colnms_not_found)>0) stop('Some columns were not found in the score file! These columns were: ', paste(collapse=' ',score_colnms_not_found))
message('Using score file columns:\n',sprintf('%10s: "%s"\n',names(score_colnms),score_colnms))

score_dt <- fread(args$score_file)[,..score_colnms] |> setnames(names(score_colnms))

if(score_dt[, !is.integer(chr) & !is.character(chr)]) message('Warning: score file\'s chromosome column is not character or integer type which is suspicious. Here is a sample: ',        paste(head(score_dt$chr),collapse=' '))
if(score_dt[, !is.integer(pos)                     ]) message('Warning: score file\'s position column is not integer type which is suspicious. Here is a sample: ',                       paste(head(score_dt$pos),collapse=' '))
if(score_dt[, any(grepl('[^ATCGNatcgn]', ref))     ]) message('Warning: score file\'s reference allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ', paste(head(score_dt$ref),collapse=' '))
if(score_dt[, any(grepl('[^ATCGNatcgn]', alt))     ]) message('Warning: score file\'s alternate allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ', paste(head(score_dt$alt),collapse=' '))
if(score_dt[, any(grepl('[^ATCGNatcgn]', ea ))     ]) message('Warning: score file\'s effect allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ',    paste(head(score_dt$ea ),collapse=' '))
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
if(geno_files_type=='bgen') { writeLines(score_dt[,paste0(chr,':', pos,'-', pos)], filter_ranges_fnm)
} else                      { writeLines(score_dt[,paste0(chr,'\t',pos,'\t',pos)], filter_ranges_fnm) }

var_extraction <- \(geno_file, filter_ranges_fnm, output_fnm) {
if(!file.exists(                                  output_fnm)) {
  cmd <-
    if(geno_files_type=='bgen') paste(
      args$bgenix_exe, '-list -g', geno_file, # -list Only list variant info, no genotype data
      '-incl-range', filter_ranges_fnm,
      '| awk \'{print $3,$4,$1,$6,$7}\'', # chr,pos,id,ref,alt columns
      '| sed -e \'/#/d\'', # Remove extra comment line # TODO eventually remove UNIX commands for native Windows compat? Soon data.table::fread() will have comment.char: https://github.com/Rdatatable/data.table/issues/856
      '>', output_fnm
    )
    else if(geno_files_type %in% c('vcf','bcf')) paste(
      args$bcftools_exe, 'query', geno_file, # -G means only list variant info, no genotype data nor header
      '-R', filter_ranges_fnm,
      '-f"%CHROM\t%POS\t%ID\t%REF\t%ALT{0}\n"',
      '>', output_fnm
    )
    else if(geno_files_type=='plink1') paste(
      args$plink2_exe, '--bfile', geno_file,
      '--rm-dup force-first', # If dups, PLINK errors unless this is here.
      '--extract bed1', filter_ranges_fnm,
      '--make-just-pvar cols= --out', sub('.pvar$','',output_fnm) # PLINK automatically adds '.pvar' to the name
    )
    else if(geno_files_type=='plink2') paste(
      args$plink2_exe, '--pfile', geno_file,
      '--rm-dup force-first',
      '--extract bed1', filter_ranges_fnm,
      '--make-just-pvar cols= --out', sub('.pvar$','',output_fnm)
    )
    else if(geno_files_type=='gds') {
      ranges <- fread(filter_ranges_fnm, col.names=c('chr','start','end'))
      gds <- seqOpen(geno_file)
      seqSetFilterPos(gds, ranges$chr, ranges$start)
      fwrite(seqGetData(gds, c('chromosome','position','annotation/id','$ref','$alt')), output_fnm)
      seqClose(gds)
      ' ' # Dummy cmd
    }
  system(cmd, ignore.stderr=T)
}}
vars_found_fnms <- paste0(args$scratch_folder,'/',basename(args$score_file),'-',basename(args$geno_files),'.pvar')

message('Checking for variants in the score_file missing from the geno_files...')
mcmapply(args$geno_files, filter_ranges_fnm, vars_found_fnms, FUN=var_extraction, mc.cores=args$threads) |> invisible()

vars_found <- do.call(rbind, lapply(vars_found_fnms, fread, col.names=c('chr','pos','id','ref','alt'))) |> suppressWarnings()
if(nrow(vars_found)>0 && all(vars_found$chr %in% score_dt$chr)) { # PLINK & seqSetFilterPos may return variants with integer-encoded chrs even if "chr#" strings were provided. Hence the additional check on chrs.
  vars_not_found <- fsetdiff(score_dt[,.(chr,pos,ref,alt)], vars_found[,.(chr,pos,ref,alt)]) # Will count as not found if ref/alt don't match
} else { vars_not_found <- score_dt } # TODO messy way of doing things

if(nrow(vars_not_found)==nrow(score_dt)) {
  unlink(vars_found_fnms)
  stop('No IDs from score_file were found in geno_files! This might be caused by:
        1. Chromosome naming mismatch between your score file and the geno type data, e.g. "1" vs. "chr1".
        2. Your ref/alt allele columns could be flipped.
        3. Your score file / genotype data positions could be of different reference genomes.')
  # TODO print out table of example chr Ids from each file.
}

vars_not_found <- merge(vars_not_found,score_dt) # We'll need the effect allele (ea) column in vars_not_found later, to determine which proxy allele corresponds to the original effect allele.
vars_not_found[, chrposid := paste0('chr',chr,':',pos) |> sub(pattern='chrchr',replacement='chr')] |> invisible() # LDproxy only takes rsIDs or chr:pos IDs

message('\x1b[33m',nrow(vars_not_found),'/',nrow(score_dt),'\x1b[m variant IDs in score file not found in the geno_files.', if(is.null(args$ldlink_token)) ' \x1b[31mWhich will be omitted because no --ldlink_token was provided for getting proxies!!\x1b[m')


# Find proxies for score file variants mising from the genotype data
if(!is.null(args$ldlink_token) & nrow(vars_not_found)>0) {
  if(nrow(vars_not_found)>100 && is.null(args$ldlink_yes_really)) stop('That\'s a lot of variants to find proxies for! Consider limiting the number of variants used to calculate your PRS by setting a p-value threshold.\nIf you really need all these proxies, please contact LDlink support at \x1b[34mNCILDlinkWebAdmin@mail.nih.gov\x1b[m to let them know you\'re planning to make a large number of API calls. If they give you the O.K. you can run this again with "--ldlink_yes_really" to skip this error message.')

  ## Don't want to ask the busy LDproxy server the same thing twice, so give the output a unique filename based on the inputs and check if it already exists before calling LDproxy.
                ldproxy_result_fnm <- paste0(args$scratch_folder,'/combined_query_snp_list_',args$ldlink_genome,'.txt') # LDproxy_batch() always writes to a file of this name
                  proxy_output_fnm <- paste0(args$scratch_folder,'/proxy-',basename(args$score_file),'-',args$ldlink_pop,'-',digest::digest(c(vars_not_found$chrposid,args$ldlink_pop,args$ldlink_r2,args$ldlink_winsize)),'.txt') # Better more specific filename that won't get overwritten in future runs on different data.
  if(!file.exists(proxy_output_fnm)) {
    old_wd <- setwd(args$scratch_folder) # LDproxy can only output in working directory, but we want it in the scratch dir
    LDproxy_batch(unique(vars_not_found$chrposid), pop=args$ldlink_pop, token=args$ldlink_token, genome_build=args$ldlink_genome, append=T, win_size=args$ldlink_winsize)
    old_wd |> setwd()
    if(file.exists(ldproxy_result_fnm)) file.rename(ldproxy_result_fnm, proxy_output_fnm) # Renaming so that the file doesn't get overwritten in future runs
  }

  if(file.exists(proxy_output_fnm)) {
    proxy_dt <- fread(proxy_output_fnm) |> suppressWarnings()
  
    proxy_dt <- proxy_dt[
      ][R2 > args$ldlink_r2
      ][!grepl('-',Alleles) # TODO ignoring these 1-bp indels b/c idk how to handle them

      ][, chr := sub(':.*','',Coord)
      ][, pos := sub('.*:','',Coord) |> as.integer()
      ][, ref := sub('\\(','',tstrsplit(Alleles,'/')[[1]])
      ][, alt := sub('\\)','',tstrsplit(Alleles,'/')[[2]])
      ][, chr_n := sub('^chr','',chr)
    ]

    proxy_dt <- merge(proxy_dt,
      vars_not_found[, .(chr_og=chr, pos_og=pos, ref_og=ref, alt_og=alt, ea_og=ea, chrposid)],
      all.x=T, allow.cartesian=T,                          by.x='query_snp', by.y='chrposid',
    ) #        allow.cartesian=T in case of multiallelic variants in the score file.

    # LDproxy only takes chr:pos, no alleles.
    # We need to ensure our alleles match what LDproxy thinks is the ref/alt for that position, and keep only those rows.
    # If we don't do this, consider a score file with two variants having the same position, but different alleles.
    #   The same proxies would be found for both, or one would be found as the proxy for the other,
    #     which would be bad because then there would be perfectly duplicate rows in the score file,
    #     meaning the variant would unfairly contribute to the score multiple times.
    #   This is a reasonable (albeit rare) case. Different alleles in a multiallelic site could all have different effects.
    #     Or, one could  simply be lift over genomic coords and discover that a site is now multiallic in the new reference genome, but you don't want to have to choose just one of the multi-alleles.
    #   Relatedly, the same proxy could be found for two variants at different position, but presumably this is okay.
    #     If so, this is more of an issue with the score_file (the variants should maybe be LD-pruned).
    proxy_dt <- proxy_dt[
      ][ mapply(Correlated_Alleles, ref_og, alt_og, FUN= \(as, ref_og, alt_og) {
           as <- strsplit(as, ',|=')[[1]] # 'A=C,T=G' -> list(A,C,T,G). A&T would be query variant's alleles, C&G the proxy's.
           as[1]==ref_og & as[3]==alt_og |
           as[3]==ref_og & as[1]==alt_og
         })

      # Figure out which allele of the proxy corresponds to the effect allele of the original variant
      ][, ea := mapply(Correlated_Alleles, ea_og, FUN = \(as, ea_og) {
                  as <- strsplit(as, ',|=')[[1]]
                  if(as[1]==ea_og) as[2] else as[4]
                })
    ]
  } else { # !file.exists(proxy_output_fnm)
    proxy_dt <- data.table() # empty
  }

  if(nrow(proxy_dt)==0) { # May happen if LDproxy found proxies for the wrong variant, i.e. for a variant that doesn't match our score_file variant's alleles.
    proxies_not_found <- vars_not_found
  } else {
    # Instead of detecting the format of the user input chromosomes in the score file, just throw both 'chr#' and plain '#' formats at the wall and one will work.
    if(geno_files_type=='bgen') { writeLines(proxy_dt[,paste0(c(chr,chr_n),':', pos,'-', pos)], filter_ranges_fnm)
    } else                      { writeLines(proxy_dt[,paste0(c(chr,chr_n),'\t',pos,'\t',pos)], filter_ranges_fnm) }

    proxies_found_fnms <- paste0(args$scratch_folder,'/',basename(proxy_output_fnm),'-',basename(args$geno_files),'.pvar')

    message('\nChecking for ~proxy~ variants in the geno_files...')
    if(!all(file.exists(proxies_found_fnms))) mcmapply(args$geno_files, filter_ranges_fnm, proxies_found_fnms, FUN=var_extraction, mc.cores=args$threads)

    proxies_found <- do.call(rbind, proxies_found_fnms |> lapply(fread, col.names=c('chr','pos','id','ref','alt'))) |> suppressWarnings()
    proxies_found[, chr_in_geno_format := chr][, chr := paste0('chr',chr) |> sub(pattern='chrchr',replacement='chr')] |> invisible() # proxy_dt always has chr prefix, but genotype_data might not

    if(nrow(proxies_found)==0) { # May happen if none of the proxies LDproxy found are present in the geno_files
      proxies_not_found <- vars_not_found
    } else {
      proxies_found <- merge(proxies_found,proxy_dt) # var_extraction() may have picked up extraneous variants whose chr:pos matches but not ref/alt (since we only filter on ranges). Merging w/ all=F (default) eliminates those.
      proxies_best_idx <- proxies_found[, .I[which.max(R2)], by=.(chr_og,pos_og,ref_og,alt_og)]$V1
      proxies_best <- proxies_found[proxies_best_idx][, chr := chr_in_geno_format]

      score_dt2 <- proxies_best[
                     ][score_dt, on=c(chr_og='chr', pos_og='pos', ref_og='ref', alt_og='alt') # merge proxies to their original variants
                     ][is.na(chr), `:=`(chr=chr_og, pos=pos_og, ref=ref_og, alt=alt_og, ea_og=ea)
                   ]

      proxies_not_found <- fsetdiff(vars_not_found[,.(chr,       pos,       ref,       alt       )],
                                     proxies_found[,.(chr=chr_og,pos=pos_og,ref=ref_og,alt=alt_og)])
    }
  }

  if(nrow(proxies_not_found)>0) message('\x1b[31mDid not find proxies in the geno_files for \x1b[33m', nrow(proxies_not_found),'/',nrow(vars_not_found),'\x1b[31m of the score_file variants which needed proxies:\x1b[m\n', paste(capture.output(proxies_not_found),collapse='\n'), '\n\x1b[31mWarning: these variants will be omitted from the PRS calculations!\x1b[m\nIf you do not want these variants to be omitted, consider decreasing --ldlink_r2 or increasing --ldlink_winsize.\n')

  # TODO if I ever feel like writing a more detailed message about multiallelic variants
  #n_multiallelics <- proxies_not_found[, mapply(chr,pos, FUN=\(c,p) nrow(score_dt[c==chr & p==pos]))] > 1
  #if(n_multiallelics>0) { }
} else { # --ldlink_token not provided
  proxies_not_found <- vars_not_found
  if(nrow(proxies_not_found)>0) message('\x1b[31mDid not find proxies in the geno_files for \x1b[33m', nrow(proxies_not_found),'/',nrow(vars_not_found),'\x1b[31m of the score_file variants which needed proxies:\x1b[m\n', paste(capture.output(proxies_not_found),collapse='\n'), '\n\x1b[31mWarning: these variants will be omitted from the PRS calculations!\x1b[m.\n')
}


# Extract full genotype data for the variants in the --score_file.  
geno_subset_file_paths <- paste0(args$scratch_folder,'/',basename(args$score_file),'-',basename(args$geno_files))

if(geno_files_type=='bgen') { writeLines(score_dt[,paste0(chr,':', pos,'-', pos)], filter_ranges_fnm)
} else                      { writeLines(score_dt[,paste0(chr,'\t',pos,'\t',pos)], filter_ranges_fnm) }

geno_extraction <- \(geno_file, filter_ranges_fnm, geno_subset_file) {
if(!file.exists(                                   geno_subset_file)) {
  cmd <-
    if(geno_files_type=='bgen') paste(
      args$bgenix_exe, '-g', geno_file,
      '-incl-range', filter_ranges_fnm,
      '>', geno_subset_file
    )
    else if(geno_files_type %in% c('vcf','bcf')) paste(
      args$bcftools_exe, 'view', geno_file,
      '-R', filter_ranges_fnm,
      '-Ob >', geno_subset_file # output compressed BCF format
    )
    else if(geno_files_type=='plink1') paste(
      args$plink2, '--bfile', geno_file,
      '--extract bed1', filter_ranges_fnm,
      '--make-pgen --out', geno_subset_file # Might as well conver to PLINK2 .pgen now
    )
    else if(geno_files_type=='plink2') paste(
      args4plink2, '--pfile', geno_file,
      '--extract bed1', filter_ranges_fnm,
      '--make-pgen --out', geno_subset_file
    )
    else if(geno_files_type=='gds') {
      gds <- seqOpen(geno_file)
      seqSetFilterPos(gds, score_dt$chr, score_dt$pos, score_dt$ref, score_dt$alt, verbose=F)
      seqGDS2BED(gds, geno_subset_file, verbose=F) # PLINK2 cannot handle GDS files, so have to convert
      seqClose(gds)
      ' '
    }
  system(cmd, ignore.stderr=T)
}}

message('Extracting genotype data...')
mcmapply(args$geno_files, filter_ranges_fnm, geno_subset_file_paths, FUN = geno_extraction, mc.cores=args$threads) |> invisible()

if(geno_files_type == 'gds') geno_files_type <<- 'plink1' # GDS were converted to PLINK1 format.

unlink(paste0(basename(args$geno_files),'.tbi')) # Clean up index files that bcftools pulls in
                                                                                                        # v TODO v: is there any way to be sure which is the reference allele in a user-provided .bgen? Answer: no. So make it is user flag. Or implement some allele-checking algo but that'd be miserable.
if(geno_files_type  ==       'bgen'   ) plink_input_flags <- paste0('--bgen  ',geno_subset_file_paths,' ref-first --sample ',args$sample_file)
if(geno_files_type %in% c('bcf','vcf')) plink_input_flags <- paste0('--bcf   ',geno_subset_file_paths)
if(geno_files_type  ==      'plink1'  ) plink_input_flags <- paste0('--bfile ',geno_subset_file_paths)
if(geno_files_type  ==      'plink2'  ) plink_input_flags <- paste0('--pfile ',geno_subset_file_paths)
stopifnot(!is.null(plink_input_flags))
plink_output_flags <- paste0('--allow-misleading-out-arg --out ', geno_subset_file_paths)
# PRS result files will be '<scratch_folder>/<score_filename>-<geno_filename>.sscore'

score_dt_simple <- score_dt[, cbind(cpraid=paste0(sub('chr','',chr),':',pos,':',ref,':',alt), ea, .SD), .SDcols=patterns('weights[0-9]+')]
setnames(score_dt_simple, grep('weights[0-9]+',names(score_dt_simple), value=T), score_colnms[grepl('weights',names(score_colnms))]) # TODO very messy way of restoring weight column nms
score_dt_simple_fnm <- file.path(args$scratch_folder,'score_file-formatted_for_plink.csv') # TODO better run-specific nm
fwrite(score_dt_simple, score_dt_simple_fnm, sep=' ')


# Finally calculate PRSes!!
plink_prs_cmds <- paste(
  args$plink2_exe, plink_input_flags,
  '--set-all-var-ids @:#:\\$r:\\$a', # chr:pos:ref:alt
  '--score', score_dt_simple_fnm,
    'ignore-dup-ids',
    'header-read', # Use row as names for clusters
    'cols=scoresums', # Output plain sum(dosages*weights), without averaging, so that we can sum scores across chromosomes. Then we can take the average.
    'list-variants',
  '--score-col-nums', paste0('3-',ncol(score_dt_simple)), # Skip ID & effect allele columns
  #'--rm-dup force-first',
  plink_output_flags
)

message('Calculating PRSes...')
err_codes <- mclapply(plink_prs_cmds,system, ignore.stdout=T,ignore.stderr=T, mc.cores=args$threads)

if(any(err_codes!=0)) message('\nPLINK PRS calculation errors happened for some files. Check these logs:\n    ', paste0(geno_subset_file_paths[err_codes!=0],'.log',collapse='\n    '), '\n(A common "error" is that no variants were in the input file, but that may just be because you provided files for all chromosomes but your score file didn\'t have a variant in every chromosome. If so, this is not a concern.)')
# TODO: "No variants found" should not be considered an error, because it is normal that a user provides separate chr 1-22 files, but there is not a weighted variant in every chr. To avoid this, right after the LDlink step, geno_subset_files with no score file variants should be removed. Something like that.


# Sum scores from each run
prs_files      <- list.files(path=args$scratch_folder, pattern=paste0(basename(geno_subset_file_paths),'.sscore$',collapse='|')) # Why list.files and not just add .sscore? Because some .sscore files might not exist if a chr had 0 of the score file's variants in it. # TODO this will no longer be a concern if I remove geno_subset_files with no score file variants as mentioned above.
prs_file_paths <- file.path(args$scratch_folder, prs_files)
prs_sums <- rbindlist(lapply(prs_file_paths,fread))[, lapply(.SD, sum), by='#IID'] # TODO: dup variants in a geno_file will be ignored. But dup variants ~across~ geno_files will get ~added~. Add a [!duplicated(.(chr,pos,ref,alt))] step here? Could save a user if they accidentally input the same chr's geno_file twice.

fwrite(prs_sums, args$output_fnm)
message('Done! PRS results file: ', args$output_fnm)

# TODO: --q-score-range option useful for global PRS, see Kenny script for referene
#   Although that could be easily implemented here in R by allowing the user to input an optional --score_file_p_col and --p-threshold(s). Tricky part would be to do good caching of var_extraction() & geno_extraction(), so maybe leaving it to plink is better?
