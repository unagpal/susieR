num_loci <- 1
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

# Loading files for exon vs. intron vs. intergenic annotations
gene_boundaries <-  read.table(gzfile("/gpfs/commons/groups/knowles_lab/index/hg38/genes.tsv.gz"), header=TRUE)[,1:3]
gene_boundary_chr <- sub("\n", "", gene_boundaries[,1])
gene_boundary_start <- sub("\n", "", gene_boundaries[,2])
gene_boundary_end <- sub("\n", "", gene_boundaries[,3])
exon_boundaries <- read.table(gzfile("/gpfs/commons/groups/knowles_lab/index/hg38/gencode.v30.exons.txt.gz"), header=TRUE)[,1:3]
exon_boundary_chr <- sub("\n", "", exon_boundaries[,1])
exon_boundary_start <-  sub("\n", "", exon_boundaries[,2])
exon_boundary_end <-  sub("\n", "", exon_boundaries[,3])

#Obtain chr, start, end, gene of all introns/junctions with significant sQTLs
#for calculating SNP distance from junction annotations
liver_sgenes <- read.table(gzfile("/gpfs/commons/datasets/controlled/GTEx/portal/data/v8/GTEx_Analysis_v8_sQTL/Liver.v8.sgenes.txt.gz"), header=TRUE)
liver_sqtl_introns <- sub("\n", "", liver_sgenes[,1])
liver_sqtl_introns <- strsplit(liver_sqtl_introns, ":") #11049 introns/junctions

#Storing index of first & last gene/exon boundaries corresponding to each chromosome
#To expedite obtaining exon vs. intron vs. intergenic annotations below
gene_bdry_ind_by_chr <- list()
exon_bdry_ind_by_chr <- list()
unique_gene_bdry_chr <- unique(gene_boundary_chr)
unique_exon_bdry_chr <- unique(exon_boundary_chr)
rev_gene_boundary_chr <- rev(gene_boundary_chr)
for (gene_bdry_chr in unique_gene_bdry_chr){
  gene_bdry_ind_by_chr[[gene_bdry_chr]] <- c(match(gene_bdry_chr,gene_boundary_chr), length(gene_boundary_chr)+1-match(gene_bdry_chr,rev_gene_boundary_chr))
}
rev_exon_boundary_chr <- rev(exon_boundary_chr)
for (exon_bdry_chr in unique_exon_bdry_chr){
  exon_bdry_ind_by_chr[[exon_bdry_chr]] <- c(match(exon_bdry_chr,exon_boundary_chr),length(exon_boundary_chr)+1-match(exon_bdry_chr,rev_exon_boundary_chr))
}

#Obtaining and saving annotation matrices
for (loc_num in 1:num_loci){
  geno_filename <- paste("/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Genotype_Matrices/X", loc_num, sep="_")
  geno_filename <- paste(geno_filename, ".txt", sep="")
  X_g <- read.table(geno_filename, sep = "\t")
  locus_SNP_loc <- strsplit(colnames(X_g), "_")
  locus_annotations <- matrix(0, 4+length(unique_rbps), length(locus_SNP_loc))
  locus_chr <- locus_SNP_loc[[1]][1]
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
      #Updating first_pot_peak_ctr and last_pot_peak_ctr
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
  #Annotation 3: min(|junction start pos - SNP pos|, |junction end pos - SNP pos|)/1000
  #Annotation 4: ln(min(|junction start pos - SNP pos|, |junction end pos - SNP pos|))
  print("Starting exon vs. intron vs. intergenic")
  chr_gene_bdry_range <- gene_bdry_ind_by_chr[[locus_chr]]
  chr_exon_bdry_range <- exon_bdry_ind_by_chr[[locus_chr]]
  for (SNP_ind in (1:length(locus_SNP_loc))){
    SNP_loc <- locus_SNP_loc[[SNP_ind]]
    SNP_pos <- strtoi(SNP_loc[2])
    for (gene_bdry_ind in chr_gene_bdry_range[1]:chr_gene_bdry_range[2]){
      if (SNP_pos >= gene_boundary_start[gene_bdry_ind] && SNP_pos <= gene_boundary_end[gene_bdry_ind]){
        locus_annotations[length(unique_rbps)+1,SNP_ind] = 1
      }
    }
    for (exon_bdry_ind in chr_exon_bdry_range[1]:chr_exon_bdry_range[2]){
      if (SNP_pos >= exon_boundary_start[exon_bdry_ind] && SNP_pos <= exon_boundary_end[exon_bdry_ind]){
        locus_annotations[length(unique_rbps)+2,SNP_ind] = 1
      }
    }
    min_junction_dist <- min(abs(junction_start_pos-SNP_pos),abs(junction_end_pos-SNP_pos))
    ln_min_junction_dist <- log(min_junction_dist)
    locus_annotations[length(unique_rbps)+3,SNP_ind] = min_junction_dist/1000
    locus_annotations[length(unique_rbps)+4,SNP_ind] = ln_min_junction_dist
  }
  locus_annot_fname <- paste("/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Annotation_Matrices/A_", loc_num, sep="")
  locus_annot_fname = paste(locus_annot_fname, ".txt", sep="")
  write.table(locus_annotations, file = locus_annot_fname, sep = "\t")
}

