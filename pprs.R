library(data.table)
library(LDlinkR)
library(parallel)
"%ni%" <- Negate("%in%")

# TMP: Test inputs
#commandArgs <- \(unused=T) { "--geno_files data/ukbb/*bgen$ --score_file score_files/alzheimers.csv --sample_file data/ukbb/*.sample$" }
#commandArgs <- \(unused=T) { "--geno_files data/ukbb/*bgen$ --score_file score_files/alzheimers.csv --sample_file data/ukbb/*.sample$ --score_file_id_col 6 --score_file_ref_col 3 --score_file_alt_col 4 --score_file_ea_col 7" }
commandArgs <- \(unused=T) { "--geno_files data/ukbb/*bgen$ --score_file score_files/alzheimers.csv --sample_file data/ukbb/*.sample$ --score_file_id_col cpra_id --score_file_ref_col Ref --score_file_alt_col Alt --score_file_ea_col Effect_Allele --ldlink_token 1d5a0d43036d" }
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
args$ldlink_genome <- "grch38_high_coverage" # TODO input validation, must be exactly grch38 or grch37 or grch38_high_coverage

# Parse command-line args
pieces <- tstrsplit(split="--", paste(commandArgs(T), collapse=" "))[-1]
for(piece in pieces) { parts <- unlist(tstrsplit(piece," |=")); args[[parts[[1]]]] <- parts[-1] }

# --- Input validation ---
recognized_args <- c("geno_files","sample_file","score_file","ldlink_token","ldlink_pop","ldlink_genome","scratch_dir", paste0("score_file_",c("chr","pos","id","ref","alt","ea"),"_col"))
if(any(names(args) %ni% recognized_args)) stop("Unrecognized argument(s):", paste0(" --",setdiff(names(args),recognized_args)))

## Required args provided? Not too many? They exist?
if(is.null(args$geno_files)) stop("--geno_files is required!")
if(is.null(args$score_file)) stop("--score_file is required!")

### Resolve patterns to actual files (e.g. *.bgen -> chr1.bgen, chr2.bgen)
args$geno_files  <- list.files(dirname(args$geno_files ), pattern=paste(collapse='|',basename(args$geno_files )), full.names=T)
args$score_file  <- list.files(dirname(args$score_file ), pattern=paste(collapse='|',basename(args$score_file )), full.names=T)
if(!is.null(args$sample_file)) {
args$sample_file <- list.files(dirname(args$sample_file), pattern=paste(collapse='|',basename(args$sample_file)), full.names=T)
}

if(any(lengths(args[names(args)!="geno_files"]) > 1)) stop("Only one parameter is allowed for each argument, except for --geno_files.")
if(length(args$geno_files)==0) stop("Provided --geno_files do not exist!  \n(If you're submitting a job on a compute cluster, you may need to specify absolute file paths.)")
if(length(args$score_file)==0) stop("Provided --score_file does not exist!\n(If you're submitting a job on a compute cluster, you may need to specify absolute file paths.)")

## Detect --geno_file type (and --sample_file)
if(all(grepl("bgen$|bgi$",args$geno_files))) {
  if(is.null(args$sample_file)   ) stop("--sample_file is required since your --geno_files are .bgen format!")
  if( length(args$sample_file)==0) stop("Provided --sample_file does not exist!")
  # TODO would be good to detect index files and warn if not found
  args$geno_files <- args$geno_files[!grepl("bgi$",args$geno_files)] # Don't include index files
  geno_files_type <-"bgen" 
} else if(all(grepl("bcf(.gz)?$|vcf(.gz)?$",args$geno_files))) {
  # TODO would be good to detect index files and warn if not found
  if(all(grepl("bcf(.gz)?$|tbi$|csi$"))) geno_files_type <- "bcf"
  if(all(grepl("vcf(.gz)?$|tbi$|csi$"))) geno_files_type <- "vcf"
  args$geno_files <- args$geno_files[!grepl("tbi$|csi$",args$geno_files)] # Don't include index files
} else if(all(grepl("bed$|bim$|fam$",args$geno_files))) {
  args$geno_files <- unique(sub(".bed$|.bim$|.fam$","",args$geno_files)) # PLINK only takes the file prefix
  geno_files_type <- "plink1"
} else if(all(grepl("pgen$|pvar$|psam$|bim$|fam$",args$geno_files))) {
  args$geno_files <- unique(sub("pgen$|pvar$|psam$|bim$|fam$","",args$geno_files)) # PLINK only takes the file prefix
  geno_files_type <- "plink2"
} else { stop(paste(c("Input geno_files format is unsupported, or is a mix of different formats. Below is the list of detected geno_files:",args$geno_files),collapse='\n')) }
message("Detected ",geno_files_type," files:\n    ", paste(args$geno_files,collapse="\n    "), "\n")

## Make sure --score_file's format is correct.
score_dt <- fread(args$score_file)
score_cols <- list(chr=args$score_file_chr_col,
                   pos=args$score_file_pos_col,
                   id =args$score_file_id_col,
                   ref=args$score_file_ref_col,
                   alt=args$score_file_alt_col,
                   ea =args$score_file_ea_col)

score_col_idxs <- sapply( score_cols, \(col) if(is.numeric(col)) col else which(names(score_dt)==col) )

score_col_farthest <- max(unlist(score_col_idxs))
score_cols_not_found <- score_cols[lengths(score_col_idxs )==0]
if(score_col_farthest  > ncol(score_dt)) stop("There are fewer than ",score_col_farthest," columns in the score file!")
if(length(score_cols_not_found)>0      ) stop("Some columns were not found in the score file! These columns were: ", paste(collapse=' ',score_cols_not_found))

message(paste0("Using score file's \"",names(score_dt)[score_col_idxs],"\" column as ",names(score_cols),"\n"))

score_dt <- fread(args$score_file) |> setnames(old=score_col_idxs, new=names(score_cols))

if(score_dt[, !is.integer(chr) & !is.character(chr)]) message("Warning: score file's chromosome column is not character or integer type which is suspicious. Here is a sample: ",        paste(head(score_dt$chr),collapse=' '))
if(score_dt[, !is.integer(pos)                     ]) message("Warning: score file's position column is not integer type which is suspicious. Here is a sample: ",                       paste(head(score_dt$pos),collapse=' '))
if(score_dt[, any(grepl("[^ATCGNatcgn]", ref))     ]) message("Warning: score file's reference allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ", paste(head(score_dt$ref),collapse=' '))
if(score_dt[, any(grepl("[^ATCGNatcgn]", alt))     ]) message("Warning: score file's alternate allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ", paste(head(score_dt$alt),collapse=' '))
if(score_dt[, any(grepl("[^ATCGNatcgn]", ea ))     ]) message("Warning: score file's effect allele column has non-nucleotide letters in it, which is suspicious. Here is a sample: ",    paste(head(score_dt$ea ),collapse=' '))
if(score_dt[,!any(grepl("[0-9]",         id ))     ]) message("Warning: score file's id column has no numbers in it, which is suspicious, here is a sample: ",                           paste(head(score_dt$id ),collapse=' ')) 
if(score_dt[,              anyDuplicated(id)>0     ]) stop("Found duplicated IDs in the score file! Please remove or merge the duplicates.") # TODO suggest data.table by= or tibble group_by code maybe

weights_cols <- score_dt[, (score_col_farthest+1):ncol(score_dt) ]
message('Using score file weight columns:\n    "',paste(names(weights_cols),collapse='"\n    "'),'"')

## LDlink inputs
###TODO: If neither token nor pop provided, warn that no proxies will be gotten if some variants are missing from geno_files. Complain if not both token & pop provided. If pop not provided, LDlinkR::list_pop(), but warn that ALL or multiple pops may be slower. 
if(is.null(args$ldlink_token)) message("\x1b[31m No LDlink API token provided\x1b[m, so no proxies will be gotten for missing variants.")
if(any(args$ldlink_pop %ni% list_pop()$pop_code)) stop("Invalid LDlink population code(s):") # TODO finish this, and print out list_pop() for the user

# --- Finished input validation ---

dir.create(args$scratch_dir, recursive=T, showWarnings=F)

# Get list of score file variant IDs present in the geno_files (not extracting genotype data yet). Used for getting proxies later.
# Why not use PLINK to read all formats? Because it cannot make use of .bgi (bgen) or .csi (vcf/bcf) index files. Instead it would try to convert everything to plink2 format first, (very slow for huge datasets).
ranges_file    <- file.path(args$scratch_dir,"score_variant_ranges.txt"   )
score_ids_file <- file.path(args$scratch_dir,"score_variant_ids.txt"      )
writeLines(score_dt$id, score_ids_file)
writeLines(score_dt[,paste0(chr,"\t",pos)], ranges_file)

found_id_files <- paste0(args$scratch_dir,"/",basename(args$score_file),"-",basename(args$geno_files),".snplist")

id_extraction_cmds <- mapply(args$geno_files, found_id_files, FUN = \(gf,idf) {
  if(geno_files_type=="bgen") paste0(
    "bgenix -list -g ", gf, # -list Only list variant info, no genotype data
    #" -incl-range ", ranges_file, # Purposefully commented out. bgenix extracts the _union_ of -incl-range and -incl-rsids. This can lead to variants we don't want which are at the same position as the desired variants. Doesn't help speed because bgen files are already indexed by ID, not only chr+pos. 
    " -incl-rsids ", score_ids_file,
    " | cut -f1", # Keep only ID column
    " | sed '/# bgenix/d' | sed '/alternate_ids/d'", # Remove comment and colname lines
    " > ", idf
  )
  else if(geno_files_type %in% c("vcf","bcf")) paste0(
    "bcftools view -GH", gf, # -GH means only list variant info, no genotype data nor header
    " -R ", ranges_file, # Need to filter by ranges otherwise will be very slow, because VCF/BCF files are only indexed by chr+pos, not ID. Unlike bgenix, -R combined with a -i filter will take the intersection of variants.
    " -i'ID=@",score_ids_file,"'",
    " -Ob", # Output compressed BCF format
    " | cut -f3",
    " > ", idf
  )
  else if(geno_files_type=="plink1") paste0(
    "plink2 --bfile", gf,
    " --rm-dup force-first", # If dups, PLINK errors unless this is here.
    " --extract ", score_ids_file,
    " --write-snplist --out ", sub(".snplist$","",idf) # PLINK automatically adds ".snplist" to the name
  )
  else if(geno_files_type=="plink2") paste0(
    "plink2 --pfile", gf,
    " --rm-dup force-first",
    " --extract ", score_ids_file,
    " --write-snplist --out ", sub(".snplist$","",idf)
  )
  else stopifnot(F)
})

if(!all(file.exists(found_id_files))) { # Don't rerun this if the found id files alredy exist, bottleneck.
  invisible(mclapply(id_extraction_cmds, system))
}

ids_found <- sapply(found_id_files, scan, what=character()) |> unlist()
ids_not_found <- setdiff(score_dt$id, ids_found)

if(length(ids_found)==0) stop("No IDs from score_file were found in geno_files! This is probably an ID format mismatch issue, e.g. rsID vs. chr:pos ID, etc..") # TODO print out a couple Ids from each file.

# TODO: Proxies. See what variants are missing from the genotype files, grab proxies, check if any of _those_ are missing, repeat. Don't forget to update score file df with the proxies.
ids_not_found_chrpos <- score_dt[id %in% ids_not_found, paste0("chr",chr,":",pos)] |> sub(pattern="chrchr",replacement="chr")

# TODO: rename file or s/t and detect when already exists s.t. don't cal LDlink's servers more than you need to
old_wd <- setwd(args$scratch_dir) # LDproxy outputs a file in working directory
LDproxy_batch(ids_not_found_chrpos, pop=args$ldlink_pop, token=args$ldlink_token, genome_build=args$ldlink_genome, append=T, win_size=100000)
setwd(old_wd)
proxy_dt <- fread(paste0(args$scratch_dir,"/combined_query_snp_list_",args$ldlink_genome,".txt"))
# TODO: take best proxy for each missing variant, overwrite the "score_ids_file", then check for variants that are STILL missing... put this in a while loop or recurse


# Write a new score file updated with proxies, and might as well reduce the columns to only those needed for PLINK --score input to simplify the code later.
## TODO maybe don't make a "simple" score file version, instead use score_col_farthest to tell PLINK what score cols. That way in the scratch folder, the user can find their score file updated w/ proxies if they want, w/ all their other additional metadata still there as well.
score_dt_simple <- score_dt[, .SD, .SDcols=c(args$col_idxs["id"],args$col_idxs["ea"], max(args$col_idxs+1):ncol(score_dt)) ]
score_file_path_simple <- file.path(args$scratch_dir,"score_file-formatted_for_plink-with_proxies.csv")
fwrite(score_dt_simple, score_file_path_simple, sep=' ')

# Extract full genotype data for the variants in the --score_file.  
geno_subset_file_paths <- paste0(args$scratch_dir,"/",basename(args$score_file),"-",basename(args$geno_files))

geno_extraction_cmds <- mapply(args$geno_files, geno_subset_file_paths, FUN=\(gf,sgf) {
  if(geno_files_type=="bgen") paste0(
    "bgenix -g ", gf,
    #" -incl-range ", ranges_file, # Purposefully commented out. bgenix extracts the _union_ of -incl-range and -incl-rsids. This can lead to variants we don't want which are at the same position as the desired variants. Doesn't help speed because bgen files are already indexed by ID, not only chr+pos. 
    " -incl-rsids ", ids_file,
    " > ", sgf
  )
  else if(geno_files_type %in% c("vcf","bcf")) paste0(
    "bcftools view ", gf,
    " -R ", ranges_file, # Need to filter by ranges otherwise will be very slow, because VCF/BCF files are only indexed by chr+pos, not ID. Unlike bgenix, -R combined with a -i filter will take the intersection of variants.
    " -i'ID=@",ids_file,"'",
    " -Ob", # Output compressed BCF format
    " > ", sgf
  )
  else if(geno_files_type=="plink1") paste0(
    "plink2 --bfile", gf,
    " --extract ", ids_file,
    " --make-pgen --out ", sgf # Might as well conver to PLINK2 .pgen now
  )
  else if(geno_files_type=="plink2") paste0(
    "plink2 --pfile", gf,
    " --extract ", ids_file,
    " --make-pgen --out ", sgf
  )
  else stopifnot(F)
})
message("\nRunning the following commands to extract genotype data...:\n    ", paste(geno_extraction_cmds,collapse="\n    "))
  
if(!all(file.exists(geno_subset_file_paths))) { # Don't rerun this if the subsetted files alredy exist, big bottleneck.
  invisible(mclapply(geno_extraction_cmds,system,ignore.stderr=T))
}
                                                                                              # v TODO v
if(geno_files_type=="bgen")   plink_input_flags <- paste0("--bgen  ",geno_subset_file_paths," ref-unknown --sample ",args$sample_file)
if(geno_files_type=="bcf")    plink_input_flags <- paste0("--bcf   ",geno_subset_file_paths)
if(geno_files_type=="vcf")    plink_input_flags <- paste0("--vcf   ",geno_subset_file_paths)
if(geno_files_type=="plink1") plink_input_flags <- paste0("--bfile ",geno_subset_file_paths)
if(geno_files_type=="plink2") plink_input_flags <- paste0("--pfile ",geno_subset_file_paths)
stopifnot(!is.null(plink_input_flags))
plink_output_flags <- paste0("--allow-misleading-out-arg --out ", geno_subset_file_paths)
# PRS result files will be "<scratch_dir>/<score_filename>-<geno_filename>.sscore

# Finally calculate PRSes!!
plink_prs_cmds <- paste(
  "plink2", plink_input_flags,
  "--score", score_file_path_simple,
    "ignore-dup-ids",
    "header-read", # Use row as names for clusters
    "cols=scoresums", # Output plain sum(dosages*weights), without averaging, so that we can sum scores across chromosomes. Then we can take the average.
    "list-variants",
  "--score-col-nums", paste0("3-",ncol(score_dt_simple)), # Skip ID & effect allele columns
  #"--rm-dup force-first",
  plink_output_flags
)
message("\nRunning the following commands to calculate PRSes...:\n    ", paste(plink_prs_cmds,collapse="\n    "))

err_codes <- mclapply(plink_prs_cmds,system, ignore.stdout=T,ignore.stderr=T)

if(any(err_codes!=0)) message("\nPLINK PRS calculation errors happened for some files. Check these logs:\n    ", paste0(geno_subset_file_paths[err_codes!=0],".log",collapse="\n    "))
# TODO: No variants found should not be considered an error, because it is normal that a user provides separate chr 1-22 files, but there is not a weighted variant in every chr. To avoid this, right after the LDlink step, geno_subset_files with no score file variants should be removed. Something like that.


# Sum scores from each run
prs_files      <- list.files(path=args$scratch_dir, pattern=paste0(basename(geno_subset_file_paths),".sscore",collapse='|')) # Why list.files and not just add .sscore? Because some .sscore files might not exist if a chr had 0 of the score file's variants in it.
prs_file_paths <- file.path(args$scratch_dir, prs_files)
prs_sums <- rbindlist(lapply(prs_file_paths,fread))[, lapply(.SD, sum), by="#IID"]
#prs_sums[, 3:ncol(prs_sums) := lapply(.SD,'/',ALLELE_CT), .SDcols = 3:ncol(prs_sums)] # [.SD := lapply(.SD,'/',ALLELE_CT), .SDcols=-c("#IID","ALLELE_CT")] should work but doesn't I think because old-ish version of data.table on Broad server
prs_sums

fwrite(prs_sums, "my_results.txt") #TMP


# TODO: --q-score-range option useful for global PRS, see Kenny script for referene
