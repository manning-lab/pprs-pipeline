Rscript ../pprs.R "
  --geno_files data/1kG_gds/*.gds

  --score_file score_files/alzheimers.csv

  --score_file_chr_col Chr_n
  --score_file_pos_col Pos
  --score_file_ref_col Ref
  --score_file_alt_col Alt
  --score_file_ea_col Effect_Allele
  --score_file_weight_cols Apolipoprotein_B_Pos CRP_Neg Diffusivity_inferior_longitudinal_fasciculus_Neg ICVF_Cingulum_CingulateR Proinsulin_Neg T2stat_PutamenL_Neg Total_GRS

  --proxy_sample_pop GBR

  --output_fnm alzheimers-1000G.txt
"
