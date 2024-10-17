library(data.table)

args <- list() # Parse command-line args
pieces <- tstrsplit(split="--", paste(commandArgs(T), collapse=" "))[-1]
for(piece in pieces) { parts <- unlist(tstrsplit(piece," ")); args[[parts[[1]]]] <- parts[-1] }

if(is.null(args$geno_file)) stop("--geno_file is required!")
if(length(list.files(pattern=args$geno_file))==0) stop("Provided --geno_file does not exist! (If you're submitting a job on a compute cluster, you might need to specify absolute file paths.)")

# TODO, tricky case: what if .txt w/ lots of bgen files but no supply sample file?

if(grepl("bgen$",geno_file)) { # BGEN format
  if(   is.null  (args$sample_file)) stop("--sample_file is required since your --geno_file is bgen format!")
  if(!file.exists(args$sample_file)) stop("Provided --sample_file does not exist!")
  # bgenix cmd here
} else if(grepl("bcf$|vcf$", geno_file)) { # VCF format
  # bcftools cmd here
} else if(grepl("txt$", geno_file)) { # File specifying Multiple input genotype files
  readLines(geno_file)
} else { # PLINK1 or PLINK2 format
  # plink cmd here
}

# TODO: instead, make 
