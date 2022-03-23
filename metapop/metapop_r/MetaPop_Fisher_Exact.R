options <- commandArgs(trailingOnly = T)

input = options[1]
library_location <- options[2]

if(library_location == ""){
  library_location <- .libPaths()
}

suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))

linkage_test <- function(all_ref, second_is_ref, first_is_ref, all_snp){
  
  contingency_table <- matrix(c(all_ref, second_is_ref, first_is_ref, all_snp), nrow = 2, ncol = 2)
  
  if(any(contingency_table < 0)){
    return(list(NA, NA))
  }
  
  
  p <- fisher.test(contingency_table)$p
  top <- (all_ref*all_snp)-(second_is_ref*first_is_ref)
  botA <- sqrt(all_ref+first_is_ref)
  botB <- sqrt(all_ref+second_is_ref)
  botC <- sqrt(all_snp+first_is_ref)
  botD <- sqrt(all_snp+second_is_ref)
  botDiv <- botA*botB*botC*botD
  
  phi <- top/botDiv
  
  return(list(p, phi))
  
}

linked_data = fread(input, sep = "\t")


linkage_results <- mapply(linkage_test, linked_data$ref_count, linked_data$ref_second, linked_data$ref_first, linked_data$snp_count)

linked_data$fisher_p <- unlist(linkage_results[c(T,F)])
linked_data$phi_coef <- unlist(linkage_results[c(F,T)])

linked_data <- linked_data[!is.na(linked_data$fisher_p),]

fwrite(linked_data, input, sep = "\t", col.names = T)
