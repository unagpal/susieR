num_loci <- 11049
peak_fname <- "/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Annotation_Matrices/"
unique_rbps <- sub("\n", "",unlist(read.table(file =  paste(peak_fname, "unique_rbps.txt", sep=""), sep = "\t")))
#Obtaining max (peak end pos - peak start pos) for each RBP and chromosome
#This information expedites identifying if SNPs are in RBP peaks
max_peak_len_lst <- list()
for (rbp in unique_rbps){
  max_peak_len_lst[[rbp]] <- list()
  rbp_fname_prefix <- paste(peak_fname, rbp, sep="")
  peak_start_pos_by_chrom <- readRDS(paste(rbp_fname_prefix,"_peak_sorted_start_pos_by_chrom.txt", sep=""))
  peak_stop_pos_by_chrom <- readRDS(paste(rbp_fname_prefix,"_peak_sorted_stop_pos_by_chrom.txt", sep=""))
  for (chr in names(peak_start_pos_by_chrom)){
    max_peak_len_lst[[rbp]][[chr]] <- max(peak_stop_pos_by_chrom[[chr]] - peak_start_pos_by_chrom[[chr]])
  }
}

#Obtain chr, start, end, gene of all introns/junctions with significant sQTLs
#for calculating SNP distance from junction annotations
liver_sgenes <- read.table(gzfile("/gpfs/commons/datasets/controlled/GTEx/portal/data/v8/GTEx_Analysis_v8_sQTL/Liver.v8.sgenes.txt.gz"), header=TRUE)
liver_sqtl_introns <- sub("\n", "", liver_sgenes[,1])
liver_sqtl_introns <- strsplit(liver_sqtl_introns, ":") #11049 introns/junctions

#Loading sorted (by start position) locations of gene/exon boundaries by chromosome (see get_sorted_boundary_list.R)
#To expedite obtaining exon vs. intron vs. intergenic annotations below
boundary_fname <- "/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Annotation_Matrices/"
sorted_gene_boundary_start <- readRDS(paste(boundary_fname,"sorted_gene_boundary_start_pos.txt", sep=""))
sorted_gene_boundary_end <-  readRDS(paste(boundary_fname,"sorted_gene_boundary_end_pos.txt", sep=""))
sorted_exon_boundary_start <- readRDS(paste(boundary_fname,"sorted_exon_boundary_start_pos.txt", sep=""))
sorted_exon_boundary_end <- readRDS(paste(boundary_fname,"sorted_exon_boundary_end_pos.txt", sep=""))
#Calculating lengths of longest genes and exons by chromosome to expedite exon vs. intron vs. intergenic annotations 
max_gene_len_lst <- list()
max_exon_len_lst <- list()
for (chr in names(sorted_gene_boundary_start)){
  max_gene_len_lst[[chr]] <- max(sorted_gene_boundary_end[[chr]]-sorted_gene_boundary_start[[chr]])
  max_exon_len_lst[[chr]] <- max(sorted_exon_boundary_end[[chr]]-sorted_exon_boundary_start[[chr]])
}
#Obtaining and saving annotation matrices
for (loc_num in 1:2000){
  print(loc_num)
  geno_filename <- paste("/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Genotype_Matrices/X", loc_num, sep="_")
  geno_filename <- paste(geno_filename, ".txt", sep="")
  X_g <- read.table(geno_filename, sep = "\t")
  locus_SNP_loc <- strsplit(colnames(X_g), "_")
  locus_annotations <- matrix(0, 4+length(unique_rbps), length(locus_SNP_loc))
  locus_chr <- locus_SNP_loc[[1]][1]
  max_chr_gene_len <- max_gene_len_lst[[locus_chr]]
  max_chr_exon_len <- max_exon_len_lst[[locus_chr]]
  chr_sorted_gene_boundary_start <- sorted_gene_boundary_start[[locus_chr]]
  chr_sorted_gene_boundary_end <- sorted_gene_boundary_end[[locus_chr]]
  chr_sorted_exon_boundary_start <- sorted_exon_boundary_start[[locus_chr]]
  chr_sorted_exon_boundary_end <- sorted_exon_boundary_end[[locus_chr]]
  sqtl_intron <- liver_sqtl_introns[[loc_num]]
  junction_start_pos <- strtoi(sqtl_intron[2])
  junction_end_pos <- strtoi(sqtl_intron[3])
  #Obtaining peak annotations: for each RBP, annotation=1[SNP in peak for each HepG2 RBP]
  rbp_index <- 1
  print("Starting RBP annotations")
  for (rbp in unique_rbps){
    rbp_fname_prefix <- paste(peak_fname, rbp, sep="")
    #Obtaining RBP peaks from get_peak_list.R
    peak_start_pos_by_chrom <- readRDS(paste(rbp_fname_prefix,"_peak_sorted_start_pos_by_chrom.txt", sep=""))[[locus_chr]]
    peak_stop_pos_by_chrom <- readRDS(paste(rbp_fname_prefix,"_peak_sorted_stop_pos_by_chrom.txt", sep=""))[[locus_chr]]
    max_peak_len  <- max_peak_len_lst[[rbp]][[locus_chr]]
    num_rbp_chr_peaks <- length(peak_start_pos_by_chrom)
    first_pot_peak_ctr <- 1
    last_pot_peak_ctr <- 1
    for (SNP_ind in (1:length(locus_SNP_loc))){
      SNP_loc <- locus_SNP_loc[[SNP_ind]]
      SNP_pos <- strtoi(SNP_loc[2])
      min_peak_start_pos <- SNP_pos - max_peak_len
      while (first_pot_peak_ctr+1<=num_rbp_chr_peaks && peak_start_pos_by_chrom[first_pot_peak_ctr+1] < min_peak_start_pos){
        first_pot_peak_ctr = first_pot_peak_ctr + 1
      }
      while (last_pot_peak_ctr+1<=num_rbp_chr_peaks && peak_start_pos_by_chrom[last_pot_peak_ctr+1] <= SNP_pos){
        last_pot_peak_ctr = last_pot_peak_ctr + 1
      }
      #Checking if SNP_pos is before the end of any peak: if so, SNP is in a RBP peak
      #since for all considered peaks, SNP_pos >= peak start pos
      pot_peak_stop_pos <- peak_stop_pos_by_chrom[first_pot_peak_ctr:last_pot_peak_ctr]
      if (SNP_pos <= max(pot_peak_stop_pos)){
        locus_annotations[rbp_index, SNP_ind] = 1
      }
    }
    rbp_index = rbp_index + 1
  }
  # Obtaining two annotations specifying whether SNP is in an exon vs. intron vs. intergenic: 
  # Annotation 1: 1[SNP within gene boundaries]
  # Annotation 2: 1[SNP within exon]
  # Additionally, obtaining two annotations specifying SNP distance from nearest junction
  # Annotation 3: min(|junction start pos - SNP pos|, |junction end pos - SNP pos|)/1000
  # Annotation 4: ln(min(|junction start pos - SNP pos|, |junction end pos - SNP pos|))
  print("Starting exon vs. intron vs. intergenic")
  first_pot_gene_ctr <- 1
  last_pot_gene_ctr <- 1
  first_pot_exon_ctr <- 1
  last_pot_exon_ctr <- 1
  num_chr_genes <- length(chr_sorted_gene_boundary_start)
  num_chr_exons <- length(chr_sorted_exon_boundary_start)
  for (SNP_ind in (1:length(locus_SNP_loc))){
    SNP_loc <- locus_SNP_loc[[SNP_ind]]
    SNP_pos <- strtoi(SNP_loc[2])
    min_gene_start_pos <- SNP_pos - max_chr_gene_len
    min_exon_start_pos <- SNP_pos - max_chr_exon_len
    #Updating first and last potentially overlapping gene & exon
    while (first_pot_gene_ctr+1<=num_chr_genes && chr_sorted_gene_boundary_start[first_pot_gene_ctr+1] < min_gene_start_pos){
      first_pot_gene_ctr = first_pot_gene_ctr + 1
    }
    while (last_pot_gene_ctr+1<=num_chr_genes && chr_sorted_gene_boundary_start[last_pot_gene_ctr+1] <= SNP_pos){
      last_pot_gene_ctr = last_pot_gene_ctr + 1
    }    
    #Checking if SNP is within gene boundaries
    pot_gene_end_pos <- chr_sorted_gene_boundary_end[first_pot_gene_ctr:last_pot_gene_ctr]
    if (SNP_pos <= max(pot_gene_end_pos)){
      locus_annotations[length(unique_rbps)+1,SNP_ind] = 1
    }
    while (first_pot_exon_ctr+1<=num_chr_exons && chr_sorted_exon_boundary_start[first_pot_exon_ctr+1] < min_exon_start_pos){
      first_pot_exon_ctr = first_pot_exon_ctr + 1
    }
    while (last_pot_exon_ctr+1<=num_chr_exons && chr_sorted_exon_boundary_start[last_pot_exon_ctr+1] <= SNP_pos){
      last_pot_exon_ctr = last_pot_exon_ctr + 1
    }
    #Checking if SNP is within exon boundaries
    pot_exon_end_pos <- chr_sorted_exon_boundary_end[first_pot_exon_ctr:last_pot_exon_ctr]
    if (SNP_pos <= max(pot_exon_end_pos)){
      locus_annotations[length(unique_rbps)+2,SNP_ind] = 1
    }
    #Calculating SNP junction distance annotations
    min_junction_dist <- min(abs(junction_start_pos-SNP_pos),abs(junction_end_pos-SNP_pos))
    ln_min_junction_dist <- log(min_junction_dist)
    locus_annotations[length(unique_rbps)+3,SNP_ind] = min_junction_dist/1000
    locus_annotations[length(unique_rbps)+4,SNP_ind] = ln_min_junction_dist
  }
  locus_annot_fname <- paste("/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Annotation_Matrices/A_", loc_num, sep="")
  locus_annot_fname = paste(locus_annot_fname, ".txt", sep="")
  write.table(locus_annotations, file = locus_annot_fname, sep = "\t")
}

