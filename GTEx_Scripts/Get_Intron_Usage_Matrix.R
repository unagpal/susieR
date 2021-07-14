#Creating 11049x208 SuSiE-Ann output matrix Y
liver_leafcutter_phenotypes <- read.table(gzfile("/gpfs/commons/datasets/controlled/GTEx/portal/data/v8/GTEx_Analysis_v8_sQTL_phenotype_matrices/Liver.v8.leafcutter_phenotypes.bed.gz"))
liver_leafcutter_intron_usage <- liver_leafcutter_phenotypes[,5:dim(liver_leafcutter_phenotypes)[2]]
sqtl_leafcutter_indices <- read.table(file = "/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/sqtl_leafcutter_indices.txt", sep = "\t")
liver_sgene_intron_usage <- liver_leafcutter_intron_usage[sqtl_leafcutter_indices[,1],] 
write.table(liver_sgene_intron_usage, file = "/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Intron_Usage_Matrix.txt", sep = "\t")
