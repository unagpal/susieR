#Obtaining RBP peaks from get_peak_list.R
peak_fname <- "/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Annotation_Matrices/"
all_peak_chrom_lst <- read.table(file =  paste(peak_fname, "peak_chrom_lst.txt", sep=""), sep = "\t")
all_peak_start_lst <- read.table(file =  paste(peak_fname, "peak_start_lst.txt", sep=""), sep = "\t")
all_peak_stop_lst <- read.table(file =  paste(peak_fname, "peak_stop_lst.txt", sep=""), sep = "\t")

print(dim(all_peak_chrom_lst))

num_loci <- 1

#Obtaining annotations of A) whether SNP is in an exon vs. intron vs. intergenic: 
#Annotation 1: 1[SNP within gene boundaries]
#Annotation 2: 1[SNP within exon]
#Subsequent annotations: 1[SNP in peak for each HepG2 RBP]
gene_boundaries <-  read.table(gzfile("/gpfs/commons/groups/knowles_lab/index/hg38/genes.tsv.gz"), header=TRUE)[,1:3]
gene_boundary_chr <- sub("\n", "", gene_boundaries[,1])
gene_boundary_start <- sub("\n", "", gene_boundaries[,2])
gene_boundary_end <- sub("\n", "", gene_boundaries[,3])
exon_boundaries <- read.table(gzfile("/gpfs/commons/groups/knowles_lab/index/hg38/gencode.v30.exons.txt.gz"), header=TRUE)[,1:3]
exon_boundary_chr <- sub("\n", "", exon_boundaries[,1])
exon_boundary_start <-  sub("\n", "", exon_boundaries[,2])
exon_boundary_end <-  sub("\n", "", exon_boundaries[,3])
#for (loc_num in 1:num_loci){
for (loc_num in 1:1){
  print("Beginning new locus")
  filename <- paste("/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Genotype_Matrices/X", loc_num, sep="_")
  filename <- paste(filename, ".txt", sep="")
  X_g <- read.table(filename, sep = "\t")
  locus_SNP_loc <- strsplit(colnames(X_g), "_")
  locus_annotations <- matrix(0, 2+length(all_peak_chrom_lst), length(locus_SNP_loc))
  for (SNP_ind in (1:length(locus_SNP_loc))){
    SNP_loc <- locus_SNP_loc[[SNP_ind]]
    SNP_pos <- strtoi(SNP_loc[2])
    for (gene_ind in 1:length(gene_boundary_chr)){
      if (gene_boundary_chr[gene_ind] == SNP_loc[1]){
        if (SNP_pos>=gene_boundary_start[gene_ind] && SNP_pos<=gene_boundary_end[gene_ind]){
          locus_annotations[1,SNP_ind] = 1
        }
      }
    }
    for (exon_ind in 1:length(exon_boundary_chr)){
      if (exon_boundary_chr[exon_ind]==SNP_loc[1]){
        if (SNP_pos>=exon_boundary_start[exon_ind] && SNP_pos<=exon_boundary_end[exon_ind]){
          locus_annotations[2,SNP_ind] = 1
        }
      }
    }
    for (snp_annot_ind in 3:2+length(all_peak_chrom_lst)){
      rbp_peak_chrom_lst <- all_peak_chrom_lst[[snp_annot_ind-2]]
      rbp_peak_start_lst <- all_peak_start_lst[[snp_annot_ind-2]]
      rbp_peak_stop_lst <- all_peak_stop_lst[[snp_annot_ind-2]]
      for (file_num in 1:2){
        for (peak_ctr in 1:length(rbp_peak_chrom_lst[[file_num]])){
          if (rbp_peak_chrom_lst[[file_num]][peak_ctr]==SNP_loc[1]){
            if (SNP_pos >= strtoi(rbp_peak_start_lst[[file_num]][peak_ctr]) && SNP_pos <= strtoi(rbp_peak_stop_lst[[file_num]][peak_ctr])){
              locus_annotations[snp_annot_ind, SNP_ind] = 1
            }
          }
        }
      }
    }
  }
  locus_annot_fname <- paste("/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Annotation_Matrices/A_", loc_num, sep="")
  locus_annot_fname = paste(locus_annot_fname, ".txt", sep="")
  write.table(locus_annotations, file = locus_annot_fname, sep = "\t")
}

