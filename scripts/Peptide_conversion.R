id_conversion <- function(table, input_tissue, output_tissue) {
  
  con <- read.csv(table, sep = "\t")
  con <- con[,-3]
  names(con) <- c("Chromosome", "NCBI_ACC")
  
  tissue <- read.csv(input_tissue, sep = "\t")
  tissue.m <- merge(tissue, con, by = "Chromosome")
  tissue.m <- tissue.m[,-1]
  tissue.m <- tissue.m[,c(11,1:10)]

  
  write.table(adipose.m, file = output_tissue, quote = F, sep = "\t", col.names = T, row.names = F)
  
}

for (i in list.files(".", pattern = ".tsv")) {
  print(i)
  new <- substr(basename(i), 1, nchar(basename(i)) - 4)
  new <- paste0(new, "_converted.tsv")
  print(new)
  id_conversion(table = "ncbi_ens_chrs_mapping.out", input_tissue = i, output_tissue = new)
}
