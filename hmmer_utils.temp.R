hmmer_utils <- function (protein_fasta,hmm) {
  require('dplyr', quietly =TRUE)
  require('stringr', quietly =TRUE)
  hmm_dir <- '~/msz/bioinformatics/rippinr-master/hmm_db'
  setwd(hmm_dir)
  hmm_folder <- list.files(hmm_dir, pattern = hmm)
  setwd(hmm_folder)
  write.table(protein_fasta, 
              file="proseq.temp.txt",sep = '\t',quote = FALSE,row.names = FALSE, col.names = FALSE)
  output_name <- paste(hmm, '_profile_output.temp.tab', sep = '')
  hmm_file <- paste(hmm, '.hmm',sep = '')
  cmd <- paste('hmmscan --domtblout', output_name, hmm_file,'proseq.temp.txt') %>% 
    system(intern = TRUE)
  hmm_connect <-file(output_name)
  line <- readLines(hmm_connect,n=4,encoding='URF-8')
  t<-str_split(line, pattern = ' ')
  t[[4]] <- t[[4]][which(t[[4]]!='')]
  close(hmm_connect)
  hmm_result <- as.vector(c(str_extract_all(t[[4]][2],pattern='^PF[:digit:]+'),t[[4]][1],t[[4]][7]),mode ='character')
  file.remove(output_name, 'proseq.temp.txt')
  dirname <- '~/msz/bioinformatics/rippinr-master'
  setwd(dirname)
  return(hmm_result)
}