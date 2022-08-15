#=====================================================================#
#                                                                     #
#     RIPPINR: A TOOL SERIES FOR RIPP NATURAL PRODUCT DISCOVERY       #
#                                                                     #
#     DEVELOPER: SUZE MA   DATE: JUNE 2020                            #
#                                                                     #
#     DEPARTMENT OF CHEMISTRY, FUDAN UNIVERSITY                       #
#                                                                     #
#=====================================================================#
rentrez_args <- commandArgs(trailingOnly = T)
if (length(rentrez_args)!=0){
  message("
 #============================================================#
 #                                                            #
 #      You are currently running rentrez_utils.R script      #          
 #                                                            #
 #                  from RIPPINR tool series                  #
 #                                                            #
 #     Automating data-gathering process for genome mining    #
 #                                                            #
 #============================================================#
         ") 
  package<-c('rentrez','stringr','protr','seqinr','dplyr')
  package.check <- lapply(package,
                          FUN = function(x) {
                            if (!require(x, character.only = TRUE, warn.conflicts=FALSE,quietly=TRUE)) {
                              install.packages(x, dependencies = TRUE)
                              require(package,character.only = TRUE,warn.conflicts=FALSE,quietly=TRUE)
                            }
                          }
  )
  dirname <- '~/msz/bioinformatics/rippinr-master' %>% 
    setwd()
  source('hmmer_utils.temp.R')
  amino_acid<-read.csv("dataset/amino_acid.csv")
  aa<-amino_acid$Symbol
  scale_value<-amino_acid$Hydropathicity
  aa_GRAVY <- data.frame (aa, scale_value)
  listfile <- list.files(dirname, pattern=rentrez_args[1], full.names=TRUE) 
  if (length (listfile)==0){
    stop (paste(Sys.time(),"-", rentrez_args[1]," is not found in current working directory -"))
  } else {
    message(
      paste("-",Sys.time(),"- info - Parsing", rentrez_args[1],"-")
    ) 
  }
  acc_list <- lapply(listfile, 
                     read.file <- function(File) { 
                       read.table(File, quote="\"", comment.char="")
                     }
  )[[1]] %>% 
    unlist()
  output_name<-strsplit(rentrez_args[1],"\\.txt$")
  writeDir <- paste(output_name,"rentrez_utils")
  writetxt <- paste(output_name,"rentrez_utils_output.txt")
  writecsv <- paste(output_name,"rentrez_utils_output.csv")
  if (file.exists(writeDir)){
    message(
      paste("-", Sys.time(),"- warning - Overwriting output file for", rentrez_args[1], "-")
    )
  }else{
    dir.create(writeDir)
  }
}
rentrez_utils<-function(acc_list){
  protein_accession <- unique(acc_list) %>% as.vector()
  redundant_num <- length(acc_list)-length(protein_accession)
  if(redundant_num>0){
    message(
      paste("-",Sys.time(),"- warning - Redundancy report -", redundant_num, "redundant accession number(s) found in input list-")
    )
    setwd(writeDir)
    writeacclist <- paste(output_name,"non_redundant_accesion_list.csv")
    write.table(protein_accession,file = writeacclist,sep = '\t',quote = FALSE,row.names = FALSE, col.names = FALSE) 
    setwd(dirname)
  }
  protein_retrieve <- vector(mode="character",length=0)
  info_retrieve <- matrix(nrow = length(protein_accession),ncol = 20)
  lineage <- c('superkingdom', 'phylum','class','order','family','genus','species')
  hmm_info <- c('pfam_id','pfam_name','e-value') 
  colnames(info_retrieve) <- c('accession_number','producing_organism','gi','nuc_accession','txid',lineage,hmm_info,'aa_seq','length','molecular_weight','isoelectric_point(Pi)','hydropathicity(GRAVY)')
  num_interval<-50*c(1:1000)
  for (protein_num in 1:length(protein_accession)) { 
    if(protein_num%in%num_interval){
      interval<-60
      message(
        paste("-",Sys.time(),"- info - Suspend execution for",interval/60,"min -")
      )
      Sys.sleep(interval)
    }
    message(
      paste("-",Sys.time(),"- info - No.",protein_num, "- Fetching protein sequence for", protein_accession[protein_num],"-")
    )
    protein_retrieve[protein_num] <- entrez_fetch(db="protein", id=protein_accession[protein_num],rettype="fasta")
    link_result <- entrez_link (dbfrom='protein', id=protein_accession[protein_num], db='nuccore')
    if (length(link_result$links$protein_nuccore)!=0) {
      gi <- link_result$links$protein_nuccore[1] 
      summary_result <- entrez_summary (db='nuccore', id=gi)
      nuccore_accession <- summary_result$accessionversion
    } else{
      gi <- NA
      nuccore_accession <- NA
    }
    message(
      paste("-",Sys.time(),"- info - No.",protein_num, "- Bioinformatically scanning", protein_accession[protein_num],"-")
    )
    hmm <- hmmer_utils (protein_retrieve[protein_num], 'Pfam-A')
    seq <- str_replace_all(str_remove_all(protein_retrieve[protein_num],pattern = "^>.*\n"),pattern = "\n",'')
    length <- strsplit(seq,'') %>% unlist() %>% length()
    molweight <- strsplit(seq,'') %>% unlist() %>% pmw()
    Pi <- strsplit(seq,'') %>% unlist() %>% computePI
    gravy<-aa_GRAVY$scale_value[match(unlist(strsplit(seq,'')),aa_GRAVY$aa)]
    GRAVY<-0
    for (i in 1:length(gravy)){
      if (!is.na(gravy[i])){
        GRAVY<-GRAVY+gravy[i]
      }
    }
    message(
      paste("-",Sys.time(),"- info - No.",protein_num, "- Parsing taxonomy information for", protein_accession[protein_num],"-")
    )
    orgn_name<-unlist(str_extract_all(protein_retrieve[protein_num],"\\[.*\\]"))
    if(length(orgn_name)!=0){
      orgn<-paste(str_replace_all(orgn_name,'\\[|\\]',''),'[ORGN]')
      taxon <- entrez_search(db="taxonomy", term=orgn)
      if(length(unlist(taxon$ids))==0){
        info_retrieve[protein_num,]<-c(protein_accession[protein_num],str_replace_all(orgn_name,'\\[|\\]',''),gi,nuccore_accession,c(rep(NA,8)),hmm,seq,length,molweight,Pi,GRAVY/length)
        next
      }
      taxon_fetch <- entrez_fetch(db="taxonomy", id=taxon$ids[1], rettype="xml", parsed=TRUE)
      rank_name <- XML::xpathSApply(taxon_fetch, "//LineageEx/Taxon/ScientificName", XML::xmlValue)
      rank <- XML::xpathSApply(taxon_fetch, "//LineageEx/Taxon/Rank", XML::xmlValue)
      info_retrieve[protein_num,]<-c(protein_accession[protein_num],str_replace_all(orgn_name,'\\[|\\]',''),gi,nuccore_accession,taxon$ids[1],rank_name[match(lineage,rank)],hmm,seq,length,molweight,Pi,GRAVY/length)
      rm(taxon_fetch);rm(rank_name);rm(rank)
    }else{
      info_retrieve[protein_num,]<-c(protein_accession[protein_num],str_replace_all(orgn_name,'\\[|\\]',''),gi,nuccore_accession,c(rep(NA,8)),hmm,seq,length,molweight,Pi,GRAVY/length) 
    }
    setwd(writeDir)
    write.table(protein_retrieve,file = writetxt,sep = '\t',quote = FALSE,row.names = FALSE, col.names = FALSE) 
    write.csv(info_retrieve,file = writecsv,quote = FALSE,row.names = FALSE)
  }
  message(
    paste("-",Sys.time(),"- info - Fetching completed -", length(protein_retrieve), "non-redundant sequence processed -")
  )
}   
rentrez_utils(acc_list)



