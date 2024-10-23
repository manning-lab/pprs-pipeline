library(data.table)
library(LDlinkR)
library(parallel)

# TMP: Test inputs
commandArgs <- \(unused=T) { "--geno_files data/ukbb/*.bgen$ --score_file score_files/alzheimers.csv --sample_file data/ukbb/*.sample$" }
#commandArgs <- \(unused=T) { "--geno_files data/mxbb/*       --score_file score_files/udler2018.csv$ --score_file_pos_col 5" }

# Default arguments
args <- list()
args$scratch_dir <- "pprs_scratch"
args$score_file_chr_col <- 1
args$score_file_pos_col <- 2
args$score_file_id_col  <- 3
args$score_file_ref_col <- 4
args$score_file_alt_col <- 5
args$score_file_ea_col  <- 6

# Parse command-line args
pieces <- tstrsplit(split="--", paste(commandArgs(T), collapse=" "))[-1]
for(piece in pieces) { parts <- unlist(tstrsplit(piece," ")); args[[parts[[1]]]] <- parts[-1] }

# Resolve patterns to actual files (e.g. *.bgen -> chr1.bgen, chr2.bgen)
args$geno_files  <- list.files(dirname(args$geno_files ), pattern=basename(args$geno_files  ), full.names=T)
args$score_file  <- list.files(dirname(args$score_file ), pattern=basename(args$score_file  ), full.names=T)
if(!is.null(args$sample_file)) {
args$sample_file <- list.files(dirname(args$sample_file), pattern=basename(args$sample_file), full.names=T)
}

# Misc.
dir.create(args$scratch_dir, recursive=T, showWarnings=F)
args$col_idxs <- sapply(c(chr=args$score_file_chr_col, pos=args$score_file_pos_col, id=args$score_file_id_col, ref=args$score_file_ref_col, alt=args$score_file_alt_col, ea=args$score_file_ea_col),
                        as.integer)

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
} else if(all(grepl("bcf(.gz)?$|vcf(.gz)?$",args$geno_files))) {
  # TODO would be good to detect index files and warn if not found
  geno_files_type <- "vcf/bcf"
} else if(all(grepl("bed$|bim$|fam$",args$geno_files))) {
  args$geno_files <- unique(sub(".bed$|.bim$|.fam$","",args$geno_files)) # PLINK only takes the file prefix
  geno_files_type <- "plink1"
} else if(all(grepl("pgen$|pvar$|psam$|bim$|fam$",args$geno_files))) {
  geno_files_type <- "plink2"
} else { stop(paste(c("Input geno_files format is unsupported, or is a mix of different formats. Below is the list of detected geno_files:",args$geno_files),collapse='\n')) }
message("Detected ",geno_files_type," files:\n    ", paste(args$geno_files,collapse='\n    '))

## Make sure --score_file's format is correct.
if(anyNA(args$col_idxs)) stop(paste0("Command-line argument --score_file_",names(args$col_idxs)[is.na(args$col_idxs)],"_col"), " was invalid. Please provide integers, not column names.")
message(paste0("Using score file's \"",names(fread(args$score_file,nrow=0))[args$col_idxs],"\" column as ",names(args$col_idxs),"\n"))

args$score_file <- fread(args$score_file) |> setnames(old=args$col_idxs, new=names(args$col_idxs))

weights_cols <- args$score_file[, max(args$col_idxs+1):ncol(args$score_file)]
message("Using score file weight columns:\n    \"",paste(names(weights_cols),collapse="\"\n    \""),"\"")

anyDuplicated(args$score_file$id)
###TODO
#### Message "using names(args$score_file)[args$chr_col] as chromosome", etc.
#### Is now the time to check if the id format is the same between the geno_files and the score_file? Or wait until later?

## LDlink inputs
###TODO: If neither token nor pop provided, warn that no proxies will be gotten if some variants are missing from geno_files. Complain if not both token & pop provided. If pop not provided, LDlinkR::list_pop(), but warn that ALL or multiple pops may be slower. 

# --- Finished input validation ---



# TODO: Proxies. See what variants are missing from the genotype files, grab proxies, check if any of _those_ are missing, repeat. Don't forget to update score file df with the proxies.
# Why not use PLINK to read all formats? Because it cannot make use of .bgi (bgen) or .csi (vcf/bcf) index files. Instead it would try to convert everything to plink2 format first, (very slow for huge datasets).
ldprox_test <- LDproxy_batch(c("chr16:53785981","chr18:60217517"), token="1d5a0d43036d", genome_build="grch38")
tmp <- fread("chr16:53785981_grch38.txt")

# Extract full genotype data for the variants in the --score_file.  
# Why not use PLINK to read all formats? Because it cannot make use of .bgi (bgen) or .csi (vcf/bcf) index files. Instead it would try to convert everything to plink2 format first, (very slow for huge datasets).
ranges <- file.path(args$scratch_dir,"variant_ranges.txt")
ids    <- file.path(args$scratch_dir,"variant_ids.txt"   )
if(geno_files_type=="bgen") {
  args$score_file[, writeLines(paste0(chr,":",pos,"-",pos), ranges) ]
  args$score_file[, writeLines(id,                            ids ) ]
  mclapply(args$geno_files, \(gf) system(paste0(
    "bgenix -g ",gf,
    " -incl-range ",ranges,
    " -incl-rsids ",ids,
    " > ",args$scratch_dir,"/",basename(gf),"-subset"
  )))
} else if(geno_files_type=="vcf/bcf") {
  args$score_file[, writeLines(paste0(chr,'\t',pos), ranges) ]
  args$score_file[, writeLines(id,                     ids ) ]
  mclapply(args$geno_files, \(gf) system(paste0(
    "bcftools view ",gf,
    " -R ",ranges,
    " -i'ID=@",ids,"'",
    " > ",args$scratch_dir,"/",basename(df),"-subset"
  )))
} else if(geno_files_type=="plink1") {
  args$score_file[, writeLines(paste0(chr,":",pos,"-",pos), ranges) ]
  mclapply(args$geno_files, \(gf) system(paste0(
    "plink2 --bfile",gf,
    " --extract bed0 ",ranges,
    " --make-pgen --out ", args$scratch_dir,"/",basename(gf),"-subset" # Might as well conver to pgen now
  )))
} else if(geno_files_type=="plink2") {
  args$score_file[, writeLines(paste0(chr,":",pos,"-",pos), ranges) ]
  mclapply(args$geno_files, \(gf) system(paste0(
    "plink2 --pfile",gf,
    " --extract bed0 ",ranges,
    " --make-pgen --out ", args$scratch_dir,"/",basename(gf),"-subset"
  )))
} else { stopifnot(F) }





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
