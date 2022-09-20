
library(ultraseq)


# get path to this script
# addional functionality to R
#script_path = funr::get_script_path()
#opts_flow$load(file.path(script_path, "ultraseq.conf"))


# this is a common help arg, and will be moved to ultraseq
help_text_flowr_run = "
Optional (used by flowr::run):
def             flow definition file
platform        platform, to over-ride flow definition file
conf            tab delimited configuration file
wd              base directory used for this run
execute         logical, controlling submission of this flow to the cluster
                OR running on the local platform

use 'flowr run -h' for more details and additional optional arguments.

"

# roxygen2:::parse_file("../pipelines/bam_mutect.R", as.environment("package:ultraseq"))


#' Title
#'
#' @param tumor_bam
#' @param normal_bam
#'
#' @return
#' @export
#'
#' @examples
bam_mutect <- function(tumor_bam, 
                       normal_bam, 
                       tumor_name,
                       normal_name,# usually tumor___normal
                       out_prefix
){
  
  help_text = "

Usage:
flowr run [args for bam_mutect] [args for flowr::run]
flowr run x=bam_mutect tumor_bam=tumor.bam normal_bam=normal.bam tumor_name=TCGA-00-ABCD-01 normal_name=TCGA-00-ABCD-10 out_prefix=tumor_normal platform=local execute=FALSE


Required:
tumor_bam       tumor bam file
normal_bam      normal bam file
out_prefix      output paired name
tumor_bam       tumor sample ID
normal_bam      normal sample ID

"
  
  help_text = paste0(help_text, help_text_flowr_run)
  
  
  if(missing(tumor_bam) | missing(normal_bam))
    stop(help_text)
  
  check_args()
  
  tumor_bam = sapply(tumor_bam, tools::file_path_as_absolute)
  normal_bam = sapply(normal_bam, tools::file_path_as_absolute)
  
  # bam file names depend on the sample names we supply
  #fq_pre = fetch_pipes("fastq_preprocess")
  paired_name = paste0(tumor_name, "___", normal_name)
  #debug(preprocess)
  preproc_tum = preprocess(bam = tumor_bam, 
                           samplename = tumor_name,
                           split_by_chr = TRUE)
  
  preproc_ref = preprocess(bam = normal_bam, 
                           samplename = normal_name,
                           split_by_chr = TRUE)
  
  # add suffix to jobnme
  preproc_tum$flowmat$jobname = paste0(preproc_tum$flowmat$jobname, "_t")
  preproc_ref$flowmat$jobname = paste0(preproc_ref$flowmat$jobname, "_n")
  
  # -- run mutect on each of the smaller chrs
  #print(preproc_tum$outfiles$recalibed_bam)
  
  mut = mutect(tumor_bam = preproc_tum$outfiles,
               normal_bam = preproc_ref$outfiles,
               out_prefix = out_prefix,
               is_merged = FALSE, # input files are by chr
               split_by_chr = TRUE,
               samplename = paired_name)
  
  ## generate a bigger flowdef with all the elements
  flowmat = rbind(preproc_tum$flowmat, preproc_ref$flowmat, mut$flowmat)
  
  ## since the three modules are to be submitted as a single pipeline
  ## we will put the same sample name throughout
  flowmat$samplename = paired_name
  
  #def = to_flowdef(flowmat)
  #write_sheet(def, "analysis/fastq_mutect.def")
  
  return(list(flowmat = flowmat, outfiles = mut$outfiles))
  
  
}


#debug(bam_mutect)



# --- testing:

if(FALSE){
  
  source('../pipelines/bam_mutect.R')
  devtools::load_all()
  #debug(bam_mutect)
  #debug(mutect);
  #debug(paired_bam_set)
  #debug(preprocess)
  out = bam_mutect(
    tumor_bam="inst/testdata/bams/tumor.bam", 
    normal_bam="inst/testdata/bams/normal.bam",
    paired_name = "tumor_normal", out_prefix="tumor_normal")
  
  run("../pipelines/bam_mutect", 
      tumor_bam="tumor.bam", 
      normal_bam="normal.bam",
      paired_name = "tumor_normal",
      execute = FALSE,
      out_prefix="tumor_normal")
  
  
  
}
