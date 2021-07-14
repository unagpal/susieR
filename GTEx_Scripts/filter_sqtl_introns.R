#Obtain chr, start, end, gene of all introns/junctions
liver_leafcutter_phenotypes <- read.table(gzfile("/gpfs/commons/datasets/controlled/GTEx/portal/data/v8/GTEx_Analysis_v8_sQTL_phenotype_matrices/Liver.v8.leafcutter_phenotypes.bed.gz"))
liver_leafcutter_introns <- sub("\n", "", liver_leafcutter_phenotypes[, 4])
liver_leafcutter_introns <- strsplit(liver_leafcutter_introns, ":") #69303 introns/junctions

#Obtain chr, start, end, gene of all introns/junctions with significant sQTLs
liver_sgenes <- read.table(gzfile("/gpfs/commons/datasets/controlled/GTEx/portal/data/v8/GTEx_Analysis_v8_sQTL/Liver.v8.sgenes.txt.gz"), header=TRUE)
liver_sqtl_introns <- sub("\n", "", liver_sgenes[,1])
liver_sqtl_introns <- strsplit(liver_sqtl_introns, ":") #11049 introns/junctions

#Identify correspondence between introns/junctions with significant sQTLs
#and all Leafcutter introns/junctions
sqtl_leafcutter_indices <- rep(0, length(liver_sqtl_introns))
for (intr_ind in 1:length(liver_sqtl_introns)){
  print(intr_ind)
  sqtl_intron <- liver_sqtl_introns[[intr_ind]]
  intron_match_found <- FALSE
  leafcutter_intron_ctr <- 1
  while (intron_match_found==FALSE && leafcutter_intron_ctr<= length(liver_leafcutter_introns)){
    #Checking if chr, start, and end of introns match
    if (liver_leafcutter_introns[[leafcutter_intron_ctr]][1]==sqtl_intron[1]){
      if (liver_leafcutter_introns[[leafcutter_intron_ctr]][2]==sqtl_intron[2]){
        if (liver_leafcutter_introns[[leafcutter_intron_ctr]][3]==sqtl_intron[3]){
          sqtl_leafcutter_indices[intr_ind] = leafcutter_intron_ctr
          intron_match_found <- TRUE
        }
      }
    }
    if (intron_match_found==FALSE)
      leafcutter_intron_ctr = leafcutter_intron_ctr + 1
  }
}

write.table(sqtl_leafcutter_indices, file = "/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/sqtl_leafcutter_indices.txt", sep = "\t")
