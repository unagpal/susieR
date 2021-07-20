source("~/SuSiE-Ann/susieR/GTEx_Scripts/get_cis_geno.R")
source("~/SuSiE-Ann/susieR/GTEx_Scripts/easy_impute.R")

#Obtain Donor IDs of HepG2 Individuals (some code from https://github.com/broadinstitute/gtex-tutorials/blob/master/GTEx_data_tutorial_R.ipynb)
sample.df <- read.delim("/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt",
           as.is=TRUE, header=TRUE, row.names=1)
liver.sample.df = sample.df[sample.df['SMTSD']=='Liver', ]
liver_sample_ids <- rownames(liver.sample.df)
liver_donor_ids <- unlist(lapply(strsplit(rownames(liver.sample.df), '-'), FUN=function(x){paste(x[1],x[2],sep="-")})) #226 donor IDs

#Obtain chr, start, end, gene of all introns/junctions with significant sQTLs
liver_sgenes <- read.table(gzfile("/gpfs/commons/datasets/controlled/GTEx/portal/data/v8/GTEx_Analysis_v8_sQTL/Liver.v8.sgenes.txt.gz"), header=TRUE)
liver_sqtl_introns <- sub("\n", "", liver_sgenes[,1])
liver_sqtl_introns <- strsplit(liver_sqtl_introns, ":") #11049 introns/junctions

#Obtain genotype for each intron/junction with significant sQTLs
tab = TabixFile("/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz",
                "/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz.tbi")
window_dist <- 10e5
for (locus_ind in 1:length(liver_sqtl_introns)){
  sqtl_intron <- liver_sqtl_introns[[locus_ind]]
  locus_chr <- sqtl_intron[1]
  start_pos <- max(strtoi(sqtl_intron[2])-window_dist, 0)
  end_pos <- strtoi(sqtl_intron[3])+window_dist
  imputed_geno_matrix <- easy_impute(t(get_cis_geno(tab, locus_chr, start_pos, end_pos)))
  genotype_donor_ids <- rownames(imputed_geno_matrix)
  genotype_liver_donor_ids <- intersect(liver_donor_ids, genotype_donor_ids) #208 samples
  liver_id_genotype_indices <- match(genotype_liver_donor_ids, genotype_donor_ids)
  liver_geno_matrix <- imputed_geno_matrix[liver_id_genotype_indices, ]
  filename <- paste("/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Genotype_Matrices/X", locus_ind, sep="_")
  filename <- paste(filename, ".txt", sep="")
  #dim (liver_geno_matrix) = [# liver samples = 208] x [# SNPs]
  print(locus_ind)
  print(dim(liver_geno_matrix))
  write.table(liver_geno_matrix, file = filename, sep = "\t")
}
