
#Obtaining list of HepG2 RBPs and corresponding bed file names (before ".bed.gz")
peak_fname <- "/gpfs/commons/home/unagpal/SuSiE-Ann/Real_Data_Exp/Data_Processing/Annotation_Matrices/"
rbp_file_metadata <- read.delim("/gpfs/commons/groups/knowles_lab/data/ENCODE_eCLIP/HepG2_replicates.txt", header = TRUE,sep = "\t")
rbp_list_metadata_unparsed <- sub("\n", "", rbp_file_metadata[,6])
rbp_list_metadata <- strsplit(rbp_list_metadata_unparsed, "-human")
rbp_bed_files <- sub("\n", "",rbp_file_metadata[,1])
unique_rbps <- unique(unlist(rbp_list_metadata))
write.table(unique_rbps, file =  paste(peak_fname, "unique_rbps.txt",sep=""), sep = "\t")
num_unique_rbps <- length(unique_rbps)
#Compiling all chromosome, start pos, and end pos of RBP peaks
#Three lists below indexed by: 1) RBP index (numerically) 2) chromosome name (via $chr_name)
for (rbp_ind in 1:length(unique_rbps)){
  print(rbp_ind)
  rbp <- unique_rbps[rbp_ind]
  start_pos_by_chrom <- list()
  stop_pos_by_chrom <- list()
  for (file_num in 1:2){
    peak_filename <- paste("/gpfs/commons/groups/knowles_lab/data/ENCODE_eCLIP/bed/", rbp_bed_files[rbp_ind*2-file_num+1], sep="")
    peak_filename <- paste(peak_filename, ".bed.gz", sep="")
    rbp_peak_data <- read.table(gzfile(peak_filename))
    rbp_peak_chrom <- sub("\n", "",rbp_peak_data[,1])
    #Resolving different chromosome formats, e.g. "chr14_GL000194v1_random" vs. "chr7"
    for (peak_ind in 1:length(rbp_peak_chrom)){
      rbp_peak_chrom[[peak_ind]] <- unlist(strsplit(rbp_peak_chrom[peak_ind], "_"))[1]
    }
    if (file_num == 1){
      rbp_peak_chrom_lst <- unlist(rbp_peak_chrom)
      rbp_peak_start_lst <- sub("\n", "",rbp_peak_data[,2])
      rbp_peak_stop_lst <- sub("\n", "",rbp_peak_data[,3])
    }
    else{
      rbp_peak_chrom_lst <- c(rbp_peak_chrom_lst, unlist(rbp_peak_chrom))
      rbp_peak_start_lst <- strtoi(c(rbp_peak_start_lst,sub("\n", "",rbp_peak_data[,2])))
      rbp_peak_stop_lst <- strtoi(c(rbp_peak_stop_lst, sub("\n", "",rbp_peak_data[,3])))
    }
  }
  unique_chr_lst <- unique(rbp_peak_chrom_lst)
  for (chr_name in unique_chr_lst){
    start_pos_by_chrom[[chr_name]] <- list()
    stop_pos_by_chrom[[chr_name]] <- list()
  }
  for (peak_ind in 1:length(rbp_peak_chrom_lst)){
    peak_chr_name = rbp_peak_chrom_lst[peak_ind]
    start_pos_by_chrom[[peak_chr_name]]=c(start_pos_by_chrom[[peak_chr_name]],rbp_peak_start_lst[peak_ind])
    stop_pos_by_chrom[[peak_chr_name]]=c(stop_pos_by_chrom[[peak_chr_name]], rbp_peak_stop_lst[peak_ind])
  }
  for (chr_name in unique_chr_lst){
    sorted_start_pos_ind <- order(unlist(start_pos_by_chrom[[chr_name]]))
    start_pos_by_chrom[[chr_name]] = unlist(start_pos_by_chrom[[chr_name]][sorted_start_pos_ind])
    stop_pos_by_chrom[[chr_name]] = unlist((stop_pos_by_chrom[[chr_name]])[sorted_start_pos_ind])
  }
  rbp_peak_fname <- paste(peak_fname, rbp, sep="")
  env_start_pos <- as.environment(start_pos_by_chrom)
  env_stop_pos <- as.environment(stop_pos_by_chrom)
  save(list = ls(env_start_pos), file = paste(rbp_peak_fname, "_peak_sorted_start_pos_by_chrom.txt", sep=""), envir = env_start_pos)
  save(list = ls(env_stop_pos), file = paste(rbp_peak_fname, "_peak_sorted_stop_pos_by_chrom.txt", sep=""), envir = env_stop_pos)
}
