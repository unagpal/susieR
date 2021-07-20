output_fpath <- "/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Annotation_Matrices/"

# Loading files for exon vs. intron vs. intergenic annotations
gene_boundaries <-  read.table(gzfile("/gpfs/commons/groups/knowles_lab/index/hg38/genes.tsv.gz"), header=TRUE)[,1:3]
gene_boundary_chr <- sub("\n", "", gene_boundaries[,1])
gene_boundary_start <- strtoi(sub("\n", "", gene_boundaries[,2]))
gene_boundary_end <- strtoi(sub("\n", "", gene_boundaries[,3]))
exon_boundaries <- read.table(gzfile("/gpfs/commons/groups/knowles_lab/index/hg38/gencode.v30.exons.txt.gz"), header=TRUE)[,1:3]
exon_boundary_chr <- sub("\n", "", exon_boundaries[,1])
exon_boundary_start <-  strtoi(sub("\n", "", exon_boundaries[,2]))
exon_boundary_end <-  strtoi(sub("\n", "", exon_boundaries[,3]))

#Storing index of first & last gene/exon boundaries corresponding to each chromosome
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

#Sorting gene boundaries and exon boundaries by start position within each chromosome
#Output lists, indexed by each chromosome, contain start and end positions
sorted_gene_boundary_start <- list()
sorted_gene_boundary_end <- list()
sorted_exon_boundary_start <- list()
sorted_exon_boundary_end <- list()
for (chr_name in unique_gene_bdry_chr){
  gene_bdry_ind <- gene_bdry_ind_by_chr[[chr_name]]
  chr_gene_bdry_start <- gene_boundary_start[gene_bdry_ind[1]:gene_bdry_ind[2]]
  chr_gene_bdry_end <- gene_boundary_end[gene_bdry_ind[1]:gene_bdry_ind[2]]
  ordered_gene_start_ind <- order(chr_gene_bdry_start)
  sorted_gene_boundary_start[[chr_name]] <- chr_gene_bdry_start[ordered_gene_start_ind]
  sorted_gene_boundary_end[[chr_name]] <- chr_gene_bdry_end[ordered_gene_start_ind]
  exon_bdry_ind <- exon_bdry_ind_by_chr[[chr_name]]
  chr_exon_bdry_start <- exon_boundary_start[exon_bdry_ind[1]:exon_bdry_ind[2]]
  chr_exon_bdry_end <- exon_boundary_end[exon_bdry_ind[1]:exon_bdry_ind[2]]
  ordered_exon_start_ind <- order(chr_exon_bdry_start)
  sorted_exon_boundary_start[[chr_name]] <-chr_exon_bdry_start[ordered_exon_start_ind]
  sorted_exon_boundary_end[[chr_name]] <- chr_exon_bdry_end[ordered_exon_start_ind]
}
saveRDS(sorted_gene_boundary_start, file=paste(output_fpath, "sorted_gene_boundary_start_pos.txt", sep=""))
saveRDS(sorted_gene_boundary_end, file=paste(output_fpath, "sorted_gene_boundary_end_pos.txt", sep=""))
saveRDS(sorted_exon_boundary_start, file=paste(output_fpath, "sorted_exon_boundary_start_pos.txt", sep=""))
saveRDS(sorted_exon_boundary_end, file=paste(output_fpath, "sorted_exon_boundary_end_pos.txt", sep=""))
