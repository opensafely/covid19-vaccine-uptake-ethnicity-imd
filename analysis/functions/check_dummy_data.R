check_dummy_data <- function(
    studydef,
    custom
) {
  
  not_in_studydef <- names(custom)[!( names(custom) %in% names(studydef) )]
  not_in_custom  <- names(studydef)[!( names(studydef) %in% names(custom) )]
  
  
  if(length(not_in_custom)!=0) stop(
    paste(
      "These variables are in studydef but not in custom: ",
      paste(not_in_custom, collapse = ", ")
    )
  )
  
  if(length(not_in_studydef)!=0) stop(
    paste(
      "These variables are in custom but not in studydef: ",
      paste(not_in_studydef, collapse=", ")
    )
  )
  
  # reorder columns
  studydef <- studydef[,names(custom)]
  
  unmatched_types <- cbind(
    map_chr(studydef, ~paste(class(.), collapse=", ")),
    map_chr(custom, ~paste(class(.), collapse=", "))
  )[ (map_chr(studydef, ~paste(class(.), collapse=", ")) != map_chr(custom, ~paste(class(.), collapse=", ")) ), ] %>%
    as.data.frame() %>% rownames_to_column()
  
  if(nrow(unmatched_types)>0) stop(
    #unmatched_types
    "inconsistent typing in studydef : dummy dataset\n",
    apply(unmatched_types, 1, function(row) paste(paste(row, collapse=" : "), "\n"))
  )
  
  print("studydef and custom dummy data match!")
  
}