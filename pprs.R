deps <- c('data.table','digest','parallel','XML') 
install.packages(deps[lengths(sapply(deps,find.package,quiet=T))==0], repos='https://cloud.r-project.org')

library(data.table)
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
args$proxy_geno_file <- ifelse(file.exists('examples'), 'examples/data/1kG_plink2/all_hg38', 'data/1kG_plink2/all_hg38')
args$proxy_winsize_kb <- '100'
args$proxy_r2_cutoff <- '0.8'
args$memory_mb <- 8000
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
recognized_args <- c('geno_files','sample_file','score_file','scratch_folder','proxy_geno_file','proxy_sample_pop','proxy_winsize_kb','proxy_r2_cutoff','memory_mb','threads','output_fnm','bcftools_exe','bgenix_exe','plink2_exe','allow_allele_flips', paste0('score_file_',c('chr','pos','ref','alt','ea'),'_col'), 'score_file_weight_cols')
if(any(names(args) %ni% recognized_args)) stop('Unrecognized argument(s):', paste0(' --',setdiff(names(args),recognized_args)))

## Required args provided? Not too many? They exist?
if(is.null(args$geno_files            )) stop('--geno_files is required!\nThis may be .vcf[.gz]/.bcf[.gz], .gds, .bgen, .bed+.bim+.fam, or .pgen+.pvar+.psam file(s).')
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
if(all(grepl('\\.bgen$|\\.bgi$',args$geno_files))) {
  if(is.null(args$sample_file)   ) stop('--sample_file is required since your --geno_files are .bgen format!')
  if( length(args$sample_file)==0) stop('Provided --sample_file does not exist!')
  # TODO would be good to detect index files and warn if not found
  args$geno_files <- args$geno_files[!grepl('bgi$',args$geno_files)] # Don't include index files
  geno_files_type <-'bgen' 
} else if(all(grepl('\\.bcf(.gz)?$|\\.vcf(.gz)?$|\\.tbi$|\\.csi$',args$geno_files))) {
  # TODO would be good to detect index files and warn if not found
  if(all(grepl('bcf(.gz)?$|tbi$|csi$',args$geno_files))) geno_files_type <- 'bcf'
  if(all(grepl('vcf(.gz)?$|tbi$|csi$',args$geno_files))) geno_files_type <- 'vcf'
  args$geno_files <- args$geno_files[!grepl('\\.tbi$|\\.csi$',args$geno_files)] # Don't include index files
} else if(all(grepl('bed$|bim$|fam$',args$geno_files))) {
  args$geno_files <- unique(sub('\\.bed$|\\.bim$|\\.fam$','',args$geno_files)) # PLINK only takes the file prefix
  geno_files_type <- 'plink1'
} else if(all(grepl('\\.pgen$|\\.pvar$|\\.psam$|\\.bim$|\\.fam$',args$geno_files))) {
  args$geno_files <- unique(sub('\\.pgen$|\\.pvar$|\\.psam$|\\.bim$|\\.fam$','',args$geno_files)) # PLINK only takes the file prefix
  geno_files_type <- 'plink2'
} else if(all(grepl('\\.gds$',args$geno_files))) {
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
score_dt[, chr_n := as.numeric(sub('chr','',chr))] |> invisible()

if(score_dt[, !is.integer(chr) & !is.character(chr)]) message('Warning: score file\'s chromosome column is not character or integer type which is suspicious. Here is a sample: ',        paste(head(score_dt$chr),collapse=' '))
if(score_dt[, !is.integer(pos)                     ]) message('Warning: score file\'s position column is not integer type which is suspicious. Here is a sample: ',                       paste(head(score_dt$pos),collapse=' '))
if(score_dt[, any(grepl('[^ATCGNatcgn]', ref))     ]) message('Warning: score file\'s reference allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ', paste(head(score_dt$ref),collapse=' '))
if(score_dt[, any(grepl('[^ATCGNatcgn]', alt))     ]) message('Warning: score file\'s alternate allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ', paste(head(score_dt$alt),collapse=' '))
if(score_dt[, any(grepl('[^ATCGNatcgn]', ea ))     ]) message('Warning: score file\'s effect allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ',    paste(head(score_dt$ea ),collapse=' '))
if(score_dt[,!all(sapply(.SD,is.numeric)), .SDcols=patterns('weights')]) stop('Score file weights columns do not contain all-numeric values! Here is a sample:\n', paste(collapse='\n',capture.output(score_dt[,.SD,.SDcols=patterns('weights')])))

# TODO some sort of check if proxy_geno_file exists, and actually move the hardcoded 1kG download here as well.

# TODO (maybe): all this is hardcoded for the 1kG data from the PLINK2 resources page as the proxy_geno_file
plink_pop_flag <- ''
valid_spops  <- sort(c('EUR','EAS','AMR','SAS','AFR'))
valid_pops <- sort(c('GBR','FIN','CHS','PUR','CDX','CLM','IBS','PEL','PJL','KHV','ACB','GWD','ESN','BEB','MSL','STU','ITU','CEU','YRI','CHB','JPT','LWK','ASW','MXL','TSI','GIH'))
if(is.null(args$proxy_sample_pop)) {
  message('Not subsetting the proxy_geno_file samples (pooled sample population).')
} else {
  if(args$proxy_sample_pop %ni% union(valid_spops,valid_pops)) stop('Invalid sample population "',args$proxy_sample_pop,'" . Please choose one of:\nSuper-populations: ',paste(valid_spops,collapse=' '),'\nSub-populations: ',paste(valid_pops,collapse=' '))

  message('Subsetting the proxy_geno_file samples by population: \x1b[36m', args$proxy_sample_pop,'\x1b[m')
  if(args$proxy_sample_pop %in% valid_spops) plink_pop_flag <- paste('--keep-if SuperPop   "=="',args$proxy_sample_pop)
  if(args$proxy_sample_pop %in% valid_pops ) plink_pop_flag <- paste('--keep-if Population "=="',args$proxy_sample_pop)
}

if(!is.null(args$allow_allele_flips)) message('\x1b[36mAutomatic allele flipping enabled.\x1b[m')

# --- Finished input validation ---

dir.create(args$scratch_folder, recursive=T, showWarnings=F)

# Get list of score file variant IDs present in the geno_files (not extracting genotype data yet). Used for getting proxies later.
# Why not use PLINK to read all formats? Because it cannot make use of .bgi (bgen) or .csi (vcf/bcf) index files. Instead it would try to convert everything to plink2 format first, (very slow for huge datasets).
filter_ranges_fnm <- file.path(args$scratch_folder,'filter_ranges.txt')
if(geno_files_type=='bgen') { writeLines(score_dt[,paste0(chr,':', pos,'-', pos)], filter_ranges_fnm)
} else                      { writeLines(score_dt[,paste0(chr,'\t',pos,'\t',pos)], filter_ranges_fnm) }

var_extraction <- \(geno_file, filter_ranges_fnm, output_fnm, geno_files_type) {
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
      '--memory', args$memory_mb,
      '--rm-dup force-first', # If dups, PLINK errors unless this is here.
      '--extract bed1', filter_ranges_fnm,
      '--make-just-pvar cols= --out', sub('.pvar$','',output_fnm), # PLINK automatically adds '.pvar' to the name
      ' 2>&1 >/dev/null'
    )
    else if(geno_files_type=='plink2') paste(
      args$plink2_exe, '--pfile', geno_file,
      '--memory', args$memory_mb,
      '--rm-dup force-first',
      '--extract bed1', filter_ranges_fnm,
      '--make-just-pvar cols= --out', sub('.pvar$','',output_fnm),
      ' 2>&1 >/dev/null'
    )
    else if(geno_files_type=='gds') {
      ranges <- fread(filter_ranges_fnm, col.names=c('chr','start','end'))
      gds <- seqOpen(geno_file)
      seqSetFilterPos(gds, ranges$chr, ranges$start, verbose=F)
      fwrite(seqGetData(gds, c('chromosome','position','annotation/id','$ref','$alt')), output_fnm)
      seqClose(gds)
      ' ' # Dummy cmd
    }
  system(cmd)
  # TODO fread() the output_fnm and rewrite it back out w/ nicer colnames so I don't have to specify col.names= every time I fread it later
}
vars_found_fnms <- paste0(args$scratch_folder,'/',basename(args$score_file),'-',basename(args$geno_files),'.pvar')

message('Checking for variants in the score_file missing from the geno_files...')
mcmapply(args$geno_files, filter_ranges_fnm, vars_found_fnms, geno_files_type, FUN=var_extraction, mc.cores=args$threads) |> invisible()

vars_found <- do.call(rbind, lapply(vars_found_fnms, fread, col.names=c('chr','pos','id','ref','alt'))) |> suppressWarnings()
if(nrow(vars_found)>0 && all(vars_found$chr %in% score_dt$chr)) { # PLINK & seqSetFilterPos may return variants with integer-encoded chrs even if "chr#" strings were provided. Hence the additional check on chrs.
  vars_not_found         <- fsetdiff(score_dt[,.(chr,pos,ref,alt)], vars_found[,.(chr,pos,ref,    alt    )]) # Will count as not found if ref/alt don't match
  vars_not_found_flipped <- fsetdiff(score_dt[,.(chr,pos,ref,alt)], vars_found[,.(chr,pos,ref=alt,alt=ref)])
  if(!is.null(args$allow_allele_flips)) {
    score_dt[!vars_not_found_flipped, on=.(chr,pos,ref,alt), `:=`(ref=alt,alt=ref)]
    vars_not_found <- fintersect(vars_not_found, vars_not_found_flipped)
  }
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
proxies_not_found <- vars_not_found

message('\x1b[33m',nrow(vars_not_found),'/',nrow(score_dt),'\x1b[m variant IDs in score file not found in the geno_files.')


# Find proxies for score file variants mising from the genotype data
# TODO create 'eur_sample_ids.txt' 'afr_sample_ids.txt' etc. for later
if(nrow(vars_not_found)>0) {

  # If default 1kG proxy_geno_file and it's not already present, download it.
  if(args$proxy_geno_file %in% c('data/1kG_plink2/all_hg38','examples/data/1kG_plink2/all_hg38') & !file.exists(paste0(args$proxy_geno_file,'.pgen'))) {
      options(timeout = max(18000, getOption("timeout")))
      message('Downloading PLINK2 format 1000 Genomes data (https://www.cog-genomics.org/plink/2.0/resources) to use as a reference panel to find proxies...\nNote: this download will consume ~11GB of disk space.')
      d <- dirname(args$proxy_geno_file)
      d |> dir.create(recursive=T)
      download.file('https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1', file.path(d,'all_hg38.psam'    ))
      download.file('https://www.dropbox.com/s/ngbo2xm5ojw9koy/all_hg38_noannot.pvar.zst?dl=1',                                      file.path(d,'all_hg38.pvar.zst'))
      download.file('https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1',                                              file.path(d,'all_hg38.pgen.zst'))
      system(paste('plink2 --zst-decompress',file.path(d,'all_hg38.pgen.zst'),file.path(d,'all_hg38.pgen')))
      system(paste('plink2 --zst-decompress',file.path(d,'all_hg38.pvar.zst'),file.path(d,'all_hg38.pvar')))
      unlink(file.path(d,'all_hg38.pgen.zst'))
      unlink(file.path(d,'all_hg38.pvar.zst'))
  }

  writeLines(vars_not_found[,paste0(chr_n,'\t',pos,'\t',pos)], filter_ranges_fnm)
  vars_proxy_geno_file_has_fnm <- paste0(args$scratch_folder,'/vars_in_proxy_geno_file.pvar') # FIXME better unique name that won't confict if running multiple instances of the pipeline concurrently
  unlink(vars_proxy_geno_file_has_fnm) # Otherwise, if var_extraction finds nothing and thus doesn't return a file, file from a previous run could be mistakenly used. Could also be fixed by fixing the FIXME above.
  var_extraction(args$proxy_geno_file, filter_ranges_fnm, vars_proxy_geno_file_has_fnm, 'plink2')

  if(file.exists(vars_proxy_geno_file_has_fnm)) {
    vars_proxy_geno_file_has <- fread(vars_proxy_geno_file_has_fnm, col.names=c('chr_n','pos','id','ref','alt'))
    if(nrow(vars_proxy_geno_file_has)>0) {
      vars_proxy_geno_file_hasnt         <- fsetdiff(vars_not_found[,.(chr_n,pos,ref,alt)], vars_proxy_geno_file_has[,.(chr_n,pos,ref,    alt    )])
      vars_proxy_geno_file_hasnt_flipped <- fsetdiff(vars_not_found[,.(chr_n,pos,ref,alt)], vars_proxy_geno_file_has[,.(chr_n,pos,ref=alt,alt=ref)])
      if(!is.null(args$allow_allele_flips)) {
        score_dt[!vars_proxy_geno_file_hasnt_flipped, on=.(chr_n,pos,ref,alt), `:=`(ref=alt,alt=ref)]
        vars_proxy_geno_file_hasnt <- fintersect(vars_proxy_geno_file_hasnt, vars_proxy_geno_file_hasnt_flipped)
      }
    }
    vars_proxy_geno_file_has <- vars_proxy_geno_file_has |> merge(score_dt, by=c('chr_n','pos','ref','alt')) # Merge w/ all=F to eliminate multiallelics, since var_extraction only filters by chr:pos.
  } else { vars_proxy_geno_file_hasnt <- vars_not_found }
  # TODO merge the vars_proxy_geno_file_hasn't file w/ score_dt s.t. more info shows up in final recap messages.

  # TODO what happens if a variant is monoallelic in chosen population? I'm guessing PLINK handles it w/o fuss, but just to make sure
  message('\x1b[33m',nrow(vars_proxy_geno_file_hasnt),'/',nrow(vars_not_found),'\x1b[m were not found in the proxy_geno_file either.')
  if(nrow(vars_proxy_geno_file_hasnt) != nrow(vars_not_found)) {
    vars_to_find_proxies_for_fnm <- paste0(
      args$scratch_folder,'/vars_to_find_proxies_for',
      '-r2_',args$proxy_r2_cutoff,
      '-win_', args$proxy_winsize_kb,
      '-pop_', args$proxy_sample_pop,
      '-', digest::digest(vars_not_found, algo='md5'),
      # TODO args$proxy_geno_file too
      '.txt')
    writeLines(vars_proxy_geno_file_has[,paste0(sub('chr','',chr),':',pos,':',ref,':',alt)], vars_to_find_proxies_for_fnm) # Hardcoded no-chr-prefix naming for the PLINK2 1kG files.
    proxy_output_fnm <- paste0(vars_to_find_proxies_for_fnm,'.vcor')

    message('Finding possible proxies in proxy_geno_file, within ',args$proxy_winsize_kb,'kb windows and with r^2 > ',args$proxy_r2_cutoff,'...')
    if(!file.exists(proxy_output_fnm)) system(paste(
      args$plink2_exe, '--pfile', args$proxy_geno_file,
      plink_pop_flag,
      '--r2-unphased cols=+ref,+alt1,+maj,-id',
        '--ld-window-kb', args$proxy_winsize_kb,
        '--ld-window-r2', args$proxy_r2_cutoff,
        '--ld-snp-list', vars_to_find_proxies_for_fnm,
        '--memory', args$memory_mb,
        '--out', vars_to_find_proxies_for_fnm
    ), ignore.stdout=T)

    proxy_dt <- fread(proxy_output_fnm, col.names=c('chr_n','pos','ref','alt','maj','chr_n_proxy','pos_proxy','ref_proxy','alt_proxy','maj_proxy','r2'))

    if(nrow(proxy_dt)==0) {
      proxies_not_found <- vars_not_found
    } else {
      proxy_dt <- score_dt[,.(chr_n,chr,pos,ref,alt,ea)][proxy_dt,on=c('chr_n','pos','ref','alt')] # Merge to score_dt by chr_n, so that we can recover chr, whose names should agree with the geno_files.
      proxy_dt[, ea_proxy := fifelse((maj==ea & maj_proxy==ref_proxy), ref_proxy, alt_proxy)]

      if(nrow(proxy_dt)==0) { proxies_not_found <- vars_not_found } else {
        if(geno_files_type=='bgen') { writeLines(proxy_dt[,paste0(chr,':', pos_proxy,'-', pos_proxy)], filter_ranges_fnm)
        } else                      { writeLines(proxy_dt[,paste0(chr,'\t',pos_proxy,'\t',pos_proxy)], filter_ranges_fnm) }

        # FIXME name with a hash in it because as of now cache not invalidated if winsize or r2 cutoff changed
        proxies_found_fnms <- paste0(args$scratch_folder,'/proxies_found',basename(args$geno_files),'.pvar')

        message('Checking which of those proxies from the proxy_geno_file also exist in the geno_files...')
        mcmapply(args$geno_files, filter_ranges_fnm, proxies_found_fnms, geno_files_type, FUN=var_extraction, mc.cores=args$threads)

        proxies_found <- do.call(rbind, proxies_found_fnms |> lapply(fread, col.names=c('chr_proxy','pos_proxy','id_proxy','ref_proxy','alt_proxy'))) |> suppressWarnings()
        proxies_found[, chr_n_proxy := as.numeric(sub('chr','',chr_proxy))] |> invisible() # proxy_dt always has no chr prefix because it's always from 1kG file, but genotype_data might have it

        if(nrow(proxies_found)==0) { # May happen if none of the proxies are present in the geno_files
          proxies_not_found <- vars_not_found
        } else {
          proxies_found <- proxy_dt[proxies_found, on=c('chr_n_proxy','pos_proxy','ref_proxy','alt_proxy')] # var_extraction() may have picked up extraneous variants whose chr:pos matches but not ref/alt (since we only filter on ranges). Merging w/ all=F (default) eliminates those.
          proxies_best_idx <- proxies_found[, .I[which.max(r2)], by=.(chr,pos,ref,alt)]$V1
          proxies_best <- proxies_found[proxies_best_idx]

          score_dt <- proxies_best[
            ][score_dt, on=c('chr','pos','ref','alt') # merge proxies to their original variants
            ][!is.na(chr_n_proxy), `:=`(chr=chr_proxy, pos=pos_proxy, ref=ref_proxy, alt=alt_proxy, ea=ea_proxy) # Replace original variants that needed proxies by their proxy
          ]

          proxies_not_found <- vars_not_found[!proxies_found, on=c('chr','pos','ref','alt')]
        }
      }
    }
  }
}

message('\x1b[33m',nrow(proxies_not_found),'/',nrow(score_dt),"\x1b[m score_file variants had no candidate proxies found in the geno_files.")

# Extract full genotype data for the variants in the --score_file.  
geno_subset_file_paths <- paste0(args$scratch_folder,'/',basename(args$score_file),'-',basename(args$geno_files))

if(geno_files_type=='bgen') { writeLines(score_dt[,paste0(chr,':', pos,'-', pos)], filter_ranges_fnm)
} else                      { writeLines(score_dt[,paste0(chr,'\t',pos,'\t',pos)], filter_ranges_fnm) }

geno_extraction <- \(geno_file, filter_ranges_fnm, geno_subset_file) {
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
      args$plink2_exe, '--bfile', geno_file,
      '--extract bed1', filter_ranges_fnm,
      '--make-bed --out', geno_subset_file,
      '--memory', args$memory_mb,
      ' 2>&1 >/dev/null'
    )
    else if(geno_files_type=='plink2') paste(
      args$plink2_exe, '--pfile', geno_file,
      '--extract bed1', filter_ranges_fnm,
      '--make-pgen --out', geno_subset_file,
      '--memory', args$memory_mb,
      ' 2>&1 >/dev/null'
    )
    else if(geno_files_type=='gds') {
      gds <- seqOpen(geno_file)
      seqSetFilterPos(gds, score_dt$chr, score_dt$pos, score_dt$ref, score_dt$alt, verbose=F)
      seqGDS2BED(gds, geno_subset_file, verbose=F) # PLINK2 cannot handle GDS files, so have to convert
      seqClose(gds)
      ' '
    }
  system(cmd, ignore.stderr=T)
}

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
  '--memory', args$memory_mb,
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
message('\n\x1b[32mDone!\x1b[m PRS results file: ', args$output_fnm)

if(nrow(vars_not_found>0)) {
  message('\n\x1b[31mSUMMARY OF OMITTED VARIANTS: \x1b[33m',nrow(proxies_not_found),'/',nrow(score_dt),'\x1b[m')
  proxies_otherwise_not_found <- proxies_not_found[!vars_proxy_geno_file_hasnt, on=c('chr_n','pos','ref','alt')][,chr_n:=NULL]

  if(nrow(vars_proxy_geno_file_hasnt )>0) message(
    "\x1b[33m",nrow(vars_proxy_geno_file_hasnt ),"/",nrow(vars_not_found),
    "\x1b[31m score_file variants were not found in the geno_files nor proxy_geno_file\x1b[m:\n",
    paste(capture.output(vars_proxy_geno_file_hasnt[,chr_n:=NULL] ),collapse='\n'),
    "\n\x1b[36mTo recover these, you'll need to provide a different --proxy_geno_file.\x1b[m\n"
  )

  if(nrow(proxies_otherwise_not_found)>0) message(
    "\x1b[33m",nrow(proxies_otherwise_not_found),"/",nrow(vars_not_found),
    "\x1b[31m score_file variants were not found besides\x1b[m:\n",
    paste(capture.output(proxies_otherwise_not_found),collapse='\n'),
    "\n\x1b[36mTo recover these, you'll need to provide either:",
    "\n  1. Bigger --proxy_winsize_kb (but it will take longer to calcualte)",
    "\n  2. Lower  --proxy_r2_cufoff  (but your proxies may be less accurate)",
    "\n  3. Other  --proxy_sample_pop to subset the reference panel to a population which more closely matches the population of your geno_files",
    "\n  4. A new  --proxy_geno_file  entirely.\x1b[m\n"
  )

}
