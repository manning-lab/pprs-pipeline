library(data.table)
library(LDlinkR)
library(parallel)

# TMP: Test inputs
commandArgs <- \(unused) { "--geno_files data/ukbb/*.bgen$ --score_file score_files/alzheimers.csv --sample_file data/ukbb/*.sample$" }
commandArgs <- \(unused) { "--geno_files data/mxbb/*       --score_file score_files/udler2018.csv$" }

args <- list() # Parse command-line args
pieces <- tstrsplit(split="--", paste(commandArgs(T), collapse=" "))[-1]
for(piece in pieces) { parts <- unlist(tstrsplit(piece," ")); args[[parts[[1]]]] <- parts[-1] }

args$geno_files  <- list.files(dirname(args$geno_files ), pattern=basename(args$geno_files ), full.names=T)
args$sample_file <- list.files(dirname(args$sample_file), pattern=basename(args$sample_file), full.names=T)
args$score_file  <- list.files(dirname(args$score_file), pattern=basename(args$score_file  ), full.names=T)
if(is.null(args$scratch_dir)) { args$scratch_dir <- "pprs_scratch" }; dir.create(args$scratch_dir, showWarnings=F)

# --- Input validation ---
## Required files provided? Not too many? They exist?
if(is.null(args$geno_files))   stop("--geno_file  is required!")
if(is.null(args$score_file))   stop("--score_file is required!")
if(length(args$score_file) >1) stop("Detected multiple score   files. Please provide only one:", paste(args$score_file, collapse=' '))
if(length(args$sample_file)>1) stop("Detected multiple .sample files. Please provide only one:", paste(args$sample_file,collapse=' '))
if(length(args$geno_files)==0) stop("Provided --geno_files do not exist!  \n(If you're submitting a job on a compute cluster, you may need to specify absolute file paths.)")
if(length(args$score_file)==0) stop("Provided --score_file does not exist!\n(If you're submitting a job on a compute cluster, you may need to specify absolute file paths.)")

## Detect --geno_file type (and --sample_file)
if(all(grepl("bgen$",args$geno_files))) {
  if(is.null(args$sample_file)   ) stop("--sample_file is required since your --geno_files are .bgen format!")
  if( length(args$sample_file)==0) stop("Provided --sample_file does not exist!")
  # TODO would be good to detect index files and warn if not found
  geno_files_type <-"bgen" 
} else if(all(grepl("bcf(.gz)?$|vcf(.gz)?$", args$geno_files))) {
  message(paste(c("Detected .vcf/.bcf files:",args$geno_files),collapse='\n'))
  # TODO would be good to detect index files and warn if not found
  geno_files_type <- "vcf/bcf"
} else if(all(grepl("bed$|bim$|fam$"))) {
  message(paste(c("Detected PLINK1 files:",args$geno_files),collapse='\n'))
  geno_files_type <- "plink1"
} else if(all(grepl("pgen$|pvar$|psam$|bim$|fam$"))) {
  message(paste(c("Detected PLINK1 files:",args$geno_files),collapse='\n'))
  geno_files_type <- "plink2"
} else { stop(paste(c("Input geno_files format is unsupported, or is a mix of different formats. Below is the list of detected geno_files:",args$geno_files),collapse='\n')) }
message("Detected ",geno_files_type," files:\n", paste(args$geno_files,collapse='\n'))

## Make sure --score_file's format is correct.
args$score_file <- setnames(fread(args$score_file), old=1:6, new=c("chr","pos","ref","alt","id","ea"))
anyDuplicated(args$score_file$id)
###TODO

## LDlink inputs
###TODO: If neither token nor pop provided, warn that no proxies will be gotten if some variants are missing from geno_files. Complain if not both token & pop provided. If pop not provided, LDlinkR::list_pop(), but warn that ALL or multiple pops may be slower. 

# --- Finished input validation ---


# TODO: Proxies. See what variants are missing from the genotype files, grab proxies, check if any of _those_ are missing, repeat. Don't forget to update score file df with the proxies.

# Extract genotype data for the variants in the --score_file.  
if(geno_files_type=="bgen") {
  ranges <- file.path(args$scratch_dir,"variant_ranges.txt")
  ids    <- file.path(args$scratch_dir,"variant_ids.txt"   )
  args$score_file[, writeLines(paste0(chr,":",pos,"-",pos), ranges) ]
  args$score_file[, writeLines(id,                            ids ) ]
  mclapply(geno_files, \(gf) system(paste0(
    "bgenix -g ",gf,
    " -incl-ranges ",ranges,
    " -incl-rsids ",ids,
    " > ",args$scratch_dir,"/",gf,"-subset"
  )))
} else if(geno_files_type="vcf/bcf") {
} else if(geno_files_type="plink1") {
} else if(geno_files_type="plink2") {
} else stopifnot(F)

system(paste("bgenix -list -g",files[[16]],"-incl-range","chr16:5378598-5378599"))
system(paste("bgenix -list -g",files[[18]],"-incl-range","chr18:60217517-60217517"))

# TODO: Take the list of variants found, compare against original score file to see if any missing
#   separate fxns for bgenix and bcftools and plink extract, they should both output a df harmonzied to same format
# Then if any missing put them all into LDproxy_batch
#   TODO need to see what kind of format LDproxy_batch returns

ldprox_test <- LDproxy_batch(c("chr16:53785981","chr18:60217517"), token="1d5a0d43036d", genome_build="grch38")

tmp <- fread("chr16:53785981_grch38.txt")

# PLINK2
bgenix \
	-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
	-incl-rsids ${prs_dir}/${tag}_prs_snps.txt \
	> ${prs_dir}/${tag}_prs_subset_chr${chr}.bgen

system(paste("plink2",
	"--bgen", f, # will be looping through all the files
	"--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample",
	"--rm-dup force-first",
  "--score", args$score_file, "#? #? --score-col-nums",paste0("#?-",ncol(args$score_file)),
  "header-read ignore-dup-ids cols=scoresums"
))



	"--make-pgen",
	"--out ${prs_dir}/${tag}_prs_subset_chr${chr}"
	"--memory 5000"


# Calculate a series of P&T-based scores
plink2
  --score f
  --q-score-range ${prs_dir}/${tag}_pt_range_list.txt ${input_file} 1 8 header \
	--out ${prs_dir}/${tag}_chr${chr} \


# -------
