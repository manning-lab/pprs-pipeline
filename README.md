Mac/Linux only.
```
Rscript pprs.R "
  # (Required arguments)
  --geno_files <f1.bgen f2.bcf data/*.vcf>

  --score_file <my_cluster_weights.csv>

  --score_file_chr_col     <colname    or index    >
  --score_file_pos_col     <colname    or index    >
  --score_file_id_col      <colname    or index    >
  --score_file_ref_col     <colname    or index    >
  --score_file_alt_col     <colname    or index    >
  --score_file_ea_col      <colname    or index    >
  --score_file_weight_cols <colname(s) or index(es)>

  # (Optional arguments)
  --sample_file <my_bgen_samples.sample>

  --scratch_folder (default scratch/)

  --output_fnm (default "my_results.txt")

  --ldlink-token <e.g. 7ad621fa1bd2>
  --ldlink-pop <population code like GBR>
  --ldlink-r2 <between 0-1, proxy minimum R^2> (default 0.8)
  --ldlink-winsize <distance to search for proxies> (default 100000)
  --ldlink-genome <grch37/grch38/grch38_high_coverage> (default: grch38_high_coverage)

  --fill_vcf_ids_with <"chr:pos:ref:alt", or a vcf file with IDs to annotate with>

  --bcftools_exe <path/to/bcftools> (default "bcftools" or "tools/bcftools/bcftools")
  --bgenix_exe   <path/to/bgenix>   (default "bgenix"   or "tools/bgen/build/apps/bgenix")
  --plink2_exe   <path/to/plink2>   (default "plink2"   or "tools/plink2")

  --threads <#> (default 1)
"
```

+ **`--geno_files`**: [`.vcf[.gz]`/`.bcf[.gz]`](https://www.cog-genomics.org/plink/2.0/formats#vcf), [`gds`](https://www.bioconductor.org/packages/release/bioc/html/SeqArray.html), [`.bgen`](https://www.cog-genomics.org/plink/2.0/formats#bgen), [`.bed`](https://www.cog-genomics.org/plink/2.0/formats#bed)+[`.bim`](https://www.cog-genomics.org/plink/2.0/formats#bim)+[`.fam`](https://www.cog-genomics.org/plink/2.0/formats#fam), or [`.pgen`](https://www.cog-genomics.org/plink/2.0/formats#pgen)+[`.pvar`](https://www.cog-genomics.org/plink/2.0/formats#pvar)+[`.psam`](https://www.cog-genomics.org/plink/2.0/formats#psam) file formats.
  - Can specify (multiple) files, or patterns e.g. `data/*.bgen`.
  - VCF files can be accessed directly through URLs. See examples.
  - `--fill_vcf_ids_with`: VCF geno_files with missing IDs (such as 1000G VCFs) can be filled with "chr:pos:ref:alt" IDs, or populated with IDs from a different VCF files (probably [a VCF from dbSNP](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/)). See Examples.
  - `--sample_file`: an accompanying [`.sample`](https://www.cog-genomics.org/plink/2.0/formats#sample) file is required if using `.bgen` file(s).
  - If there are variants with duplicate chr,pos,ref,alt in your genotype data, the first will be chosen.
+ **`--score_file`** Must contain columns for chromosome, position, variant IDs, reference allele, alternate allele, effect allele, and at least one column of weights.
  - No duplicate IDs allowed.
+ `--ldlink-token`: This pipeline can use [LDproxy](https://ldlink.nih.gov/?tab=ldproxy) via the `LDlinkR` package to find variants that are highly correlated with the variants in your score file, so that if a variant is missing from your genotype file, a highly correlated proxy can be used instead. A secret token is required to access the LDproxy servers, to limit bad actors from overloading the servers. You can easily [obtain a token by registering here](https://ldlink.nih.gov/?tab=apiaccess).
  - If no token is provided, the pipeline will proceed without finding proxies, dropping missing variants.
  - `--ldlink-pop`: One or more population codes. See below. LD will be calculated from the pool of these population groups. Note that the more populations you select, the longer it will take to find proxies.
  - `--ldlink-winsize`: The larger the window size, the longer it will take to calculate proxies.

# Dependencies
+ **[R](https://cloud.r-project.org/) (>=4.1)**
  - Packages: `install.packages(c("data.table","LDlinkR","parallel","XML")`
  - (If using `.gds` files) `BiocManager::install("SeqArray")`
+ **[`plink2`](https://www.cog-genomics.org/plink/2.0/)**
+ (If using `.vcf`/`.bcf` files) [`bcftools`](http://samtools.github.io/bcftools/howtos/install.html)
+ (If using `.bgen` files) [`bgenix`](https://enkre.net/cgi-bin/code/bgen/dir?ci=tip)

`tool/get_tools.sh` can automatically install bgenix, bcftools, and plink2.

# LDproxy population codes
For use in `--ldlink-pop`. You can specify more than one or choose a super-population, but LD calculation will be slower.
| pop_code | super_pop_code |                                  pop_name |
| -------- | -------------- | ----------------------------------------- |
|      ALL |            ALL |                           ALL POPULATIONS |
|      AFR |            AFR |                                   AFRICAN |
|      YRI |            AFR |                  Yoruba in Ibadan, Nigera |
|      LWK |            AFR |                    Luhya in Webuye, Kenya |
|      GWD |            AFR |                 Gambian in Western Gambia |
|      MSL |            AFR |                     Mende in Sierra Leone |
|      ESN |            AFR |                            Esan in Nigera |
|      ASW |            AFR |   Americans of African Ancestry in SW USA |
|      ACB |            AFR |           African Carribbeans in Barbados |
|      AMR |            AMR |                         AD MIXED AMERICAN |
|      MXL |            AMR |    Mexican Ancestry from Los Angeles, USA |
|      PUR |            AMR |            Puerto Ricans from Puerto Rico |
|      CLM |            AMR |        Colombians from Medellin, Colombia |
|      PEL |            AMR |                 Peruvians from Lima, Peru |
|      EAS |            EAS |                                EAST ASIAN |
|      CHB |            EAS |              Han Chinese in Bejing, China |
|      JPT |            EAS |                  Japanese in Tokyo, Japan |
|      CHS |            EAS |                      Southern Han Chinese |
|      CDX |            EAS |       Chinese Dai in Xishuangbanna, China |
|      KHV |            EAS |         Kinh in Ho Chi Minh City, Vietnam |
|      EUR |            EUR |                                  EUROPEAN |
|      CEU |            EUR | Utah Residents from North and West Europe |
|      TSI |            EUR |                         Toscani in Italia |
|      FIN |            EUR |                        Finnish in Finland |
|      GBR |            EUR |           British in England and Scotland |
|      IBS |            EUR |               Iberian population in Spain |
|      SAS |            SAS |                               SOUTH ASIAN |
|      GIH |            SAS |  Gujarati Indian from Houston, Texas, USA |
|      PJL |            SAS |             Punjabi from Lahore, Pakistan |
|      BEB |            SAS |                   Bengali from Bangladesh |
|      STU |            SAS |              Sri Lankan Tamil from the UK |
|      ITU |            SAS |                 Indian Telugu from the UK |

# Tips
+ If you run this as a job on a computer cluster, you may need to specify absolute file paths in order for the compute node you dispatched to find your files.
+ If you get errors using VCF files over URL like "Failed to read" or "Failed to seek" or "Could not retrieve index", try running again. It is probably just be a network issue.

# Citations
* LDlink: Machiela MJ, Chanock SJ. [LDlink: a web-based application for exploring population-specific haplotype structure and linking correlated alleles of possible functional variants.](http://www.ncbi.nlm.nih.gov/pubmed/?term=26139635) Bioinformatics. 2015 Jul 2.
* PLINK2: Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) [Second-generation PLINK: rising to the challenge of larger and richer datasets.](https://doi.org/10.1186/s13742-015-0047-8) GigaScience, 4.
* BCFtools: Danecek P, Bonfield JK, et al. [Twelve years of SAMtools and BCFtools.](https://doi.org/10.1093/gigascience/giab008) Gigascience (2021) 10(2):giab008
* BGEN: Band, G. and Marchini, J., [BGEN: a binary file format for imputed genotype and haplotype data](https://doi.org/10.1101/308296) bioArxiv 308296
* SeqArray: Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir B, Laurie C, Levine D (2017). [SeqArray – A storage-efficient high-performance data format for WGS variant calls.](https://doi:10.1093/bioinformatics/btx145) Bioinformatics, 33(15), 2251-2257.

# TODO
+ Small examples for all file types
+ Remove requirement for variant IDs, just need chr,pos,ref,alt
