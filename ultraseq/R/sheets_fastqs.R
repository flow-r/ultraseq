##
## Author: Sahil Seth
## Date:   2014/09/18
## Edited: 2015/07/15
## major edits
## sseth@mdanderson.org

detect_fq_format2 <- function(x){
  
  if(grepl(".*_([ATGC]*|NoIndex)_L00([0-9]*)_R([0-9]*)_([0-9]*).fastq.gz", basename(x[1]))){
    
    flog.info("Using CASAVA 1.8 naming format")
    format <- "{{samplename}}_{{index}}_L00{{lane}}_R{{read}}_{{num}}.fastq.gz"
    
  # if any of them as have S1 in their names
  # need a better detection system
  }else if(grepl(".*_S[0-9]*_L00([0-9]*)_R([0-9]*)_([0-9]*).fastq.gz", basename(x[1]))){ 
    # miseq output
    flog.info("Using MiSeq/bcl2fastq 2.0 naming format")
    format <- "{{samplename}}_S[0-9]*_L00{{lane}}_[RI]{{read}}_{{num}}.fastq.gz"

  }else if(grepl(".*_S[0-9]*_L([0-9]*)_R([0-9]*)_([0-9]*).fastq.gz", basename(x[1]))){ 
    # "1732FL-04-03_S0_L3_R2_001.fastq.gz"
    flog.info("Using bcl2fastq 3.0?? naming format (L1, instead of L001)")
    format <- "{{samplename}}_S[0-9]*_L{{lane}}_[RI]{{read}}_{{num}}.fastq.gz"
  
    # "SRR8518122_pass_1.fastq.gz"
   }else if( grepl("SRR[0-9]*_pass_[1-2]{1}.fastq.gz", basename(x[1])) ){ 
     flog.info("Using SRA PE format")
     format <- "{{samplename}}_pass_{{read}}.fastq.gz"

  }else{
    stop(c("Looks like we could not understand pattern in names of fastq files\n",
           paste(head(basename(x)), collapse = "\n")))
  }
}


#' split_names_fastq2
#'
#' @param x a fastq file
#' @param format naming format for the file
#' @param strict_format_checking TRUE (error)/FALSE(warning)
#' 
#' @export
#'
split_names_fastq2 <- function(x, 
                               format = "{{samplename}}_{{index}}_L00{{lane}}_R{{read}}_{{num}}.fastq.gz", 
                               # regex pattern for each piece
                               lst_patterns = list(
                                 samplename = "(.*)",
                                 index = "([ATGC]*|NoIndex)",
                                 lane = "([0-9]*)",
                                 read = "([1-3]*)", ## ideally would be 1 OR 2 (rarely 3)
                                 num = "([0-9]*)"),
                               strict_format_checking = FALSE, 
                               verbose = T){
  
  
  # replace and get final pattern !
  tmp = whisker_render(format, lst_patterns)
  fmt = tmp$out
  vars = tmp$vars
  flog.debug(fmt)
  
  ## get the parse matrix
  repl <- paste("\\",1:length(vars), sep="",collapse=",")
  mat <- gsub(fmt, repl, basename(x))
  
  # longer, but more robust
  mat = lapply(1:length(mat), function(i){
    out = strsplit(mat[i], ",")[[1]]
    add_vars = length(vars) - length(out)
    if(add_vars > 0){
      out = c(out, rep(NA, add_vars))
      warning("there was a issue parsing this filename: ", basename(x[i]))
    }
    return(out)
  })
  warnings()
  
  mat <- do.call(rbind, mat)
  mat <- data.frame(mat, x, stringsAsFactors = FALSE)
  
  # test if none of the file have the correct pattern 
  # if all files were not parsed properly
  # show extra message, to help debug
  #
  # This does not catch, when some file are parsed improperly.
  if(ncol(mat) != length(c(vars, "file"))){
    message("Final evaluated regular expression: ", fmt, " ---> ", repl)
    message("Format was not able to split fastq names properly.")
    print(head(mat))
    stop("Exiting create_fq_sheet...")
  }
  
  # ------- test if even on file fails pattern ----------
  if(strict_format_checking & mean(complete.cases(mat)) < 1){
    stop("Some file names do not have the correct format, please check.",
         "\nRefer to the warning message provided below.")
  }
  
  # if all is well
  colnames(mat) <- c(vars, "file")
  
  return(mat)
}


#' check_fastq_sheet
#'
#' Checks the nuber of samples, and that number of files in each sample should be
#' same for each read (if paired end)
#'
#' @param x a data.frame with details regarding fastq files
#' @param id_column All rows, where this column is NA will be removed. 
#' @param file_column column name for fastq files
#' @param read_column column name for read
#'
#' @export
check_fastq_sheet <- function(x,
                              id_column = "samplename",
                              file_column = "file",
                              read_column = "read"){
  mat = x
  message("There are ", length(unique(mat[, id_column])), " samples in this dataset")
  dat_list <- split.data.frame(mat, mat[, id_column])
  if(length(mat$files) > 0){
    tmp <- sapply(1:length(dat_list), function(i){
      tmp <- tapply(dat_list[[i]][, file_column], dat_list[[i]][, read_column], length)
      ## check in case of paired end
      if(length(tmp) > 1){
        if(diff(tmp) != 0) stop("Number of fastq files are not the same in",
                                dat_list[[i]][, id_column][1], "sample")
      }
      return(0)
    })}
  return(x)
}

# https://stat.ethz.ch/R-manual/R-devel/library/base/html/regex.html

#' Create a table with details on names of fastq files
#' @param x path to a fastq folder
#' @param ext extensions to look for. A regex to search for fastq files
#' @param format auto detect. Pattern specify pattern acceptable to split_fastq_names. If missing will detect hiseq and miseq
#' @param strict_format_checking If some file names do not follow the format properly die with error (TRUE), continue with warning (FALSE)
#' @param sample_prefix A prefix to add to all sample names, run, project etc.....
#' @param lst_patterns a list of variables to be extracted, along with their regex patterns [typically the default should work for most cases]
#' 
#' @export
#' @import futile.logger
create_fq_sheet <- function(x,
                            ext = opts_flow$get("fastq_extension"),
                            format = opts_flow$get("fastq_format"),
                            lst_patterns = list(
                              samplename = "(.*)",
                              index = "([ATGC]*|NoIndex)",
                              lane = "([0-9]*)",
                              read = "([0-9]*)", 
                              num = "([0-9]*)"),
                            strict_format_checking = FALSE,
                            sample_prefix = ""){
  
  # get the extension to use
  if(is.null(ext)) ext = "fastq.gz$|fastq$"
  if(!grepl("$", ext, fixed = TRUE))
    ext = paste0(ext, "$") ## add a dollar
  
  message("> fetch the files to parse (", ext, ")")
  fqs <- unlist(lapply(x, list.files, pattern = ext, 
                       full.names=TRUE, recursive=TRUE))
  if(length(fqs) == 0)
    stop("no fastq files found. Please check the folder provided.")
  
  # get the format to use
  if(missing(format))
    format <- detect_fq_format2(fqs)
  
  if(is.null(format))
    format <- detect_fq_format2(fqs)
  
  fqmat <- split_names_fastq2(fqs, format, 
                              strict_format_checking = strict_format_checking, 
                              lst_patterns = lst_patterns)
  
  # if one needs to change the samplename to include anything in this
  fqmat$samplename = paste0(sample_prefix, fqmat$samplename)

  invisible(fqmat)
}

detect_fq_format <- function(x){
  .Deprecated("detect_fq_format2")
  if(grepl("_S1.*fastq.gz", x[1])){ ## miseq output
    ## ------ casava output
    message("\nUsing MiSeq naming format")
    format <- "$samplename$_S[0-9]*_L00$lane$_R$read$_$num$.fastq.gz"
  }else if(grepl(".*_([ATGC]*|NoIndex).*L00([0-9]*)_R([0-9]*)_([0-9]*).fastq.gz",basename(x[1]))){
    ## ------ casava output
    message("\nUsing CASAVA naming format")
    format <- "$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"
  }else{
    stop(c("Looks like we could not understand pattern in names of fastq files",
           print(head(x))))
  }
}





#' @title split.names.fastq
#' @description Given a format split the files provided. Tried and tested on fastq files
#'
#' @param files a character vector of filenames, to be parsed
#' @param format the regex type format to be used to extract information from filenames. See default format as an example.
#'
#' @export
split_names_fastq <- function(files, 
                              format="$samplename$_$index$_L00$lane$_R$read$_$num$.fastq.gz"){
  ## process format:
  .Deprecated("split_names_fastq2")
  fmt <- gsub("\\$samplename\\$","(.*)",format)
  fmt <- gsub("\\$index\\$","([ATGC]*|NoIndex)", fmt)
  fmt <- gsub("\\$lane\\$","([0-9]*)", fmt)
  fmt <- gsub("\\$read\\$","([0-9]*)", fmt)
  fmt <- gsub("\\$num\\$","([0-9]*)", fmt)
  ## -- column names
  out = strsplit(format,"\\$")[[1]]
  ## every second one would be a variable
  cols = out[seq(2,length(out), by=2)]
  repl <- paste("\\",1:length(cols),sep="",collapse=",")
  mat <- gsub(fmt, repl,basename(files))
  mat <- data.frame(cbind(do.call(rbind,strsplit(mat,",")),files))
  colnames(mat) <- c(cols,"files")
  return(mat)
}
split.names.fastq=split_names_fastq






if(FALSE){
  x = "/rsrch1/genetics/htep/hiseq/150707_ST-J00106_0022_AH35LGBBXX/Unaligned/YLD-NS-KO_RNA109"
  fqmat = create_fq_sheet(x)
  fqmat = check_fastq_sheet(fqmat)
  
}

########################## ------------------- functions below are not used
#      are here for archival purposes only





if(FALSE){
  
  path = "/rsrch1/genetics/htep/hiseq/150707_ST-J00106_0022_AH35LGBBXX/Unaligned/YLD-NS-KO_RNA109"
  runid= "140626_M01692_0046_000000000-A99M8"
  runid="140627_M01692_0047_000000000-A9J0V";
  #debug(create_sample_mat)
  
  path=sprintf("/scratch/iacs/gcc/leveli/%s",runid)
  fq_mat <- create_sample_mat(path,project="ANDY",subproject="Futreal-AS", runid=runid,
                              outpath="~/flows/ANDY-Futreal-AS")
}

