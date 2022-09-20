
#' Flow following Broad's best practices for variant calling, starting from sorted bam
#'
#' @title Pre-process bam files following Broad's best practices for variant calling, starting from aligned BAM file
#' @description This function provides a wrapper around the best practices described on \href{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}{GATK's website}.
#' If the link is broken google 'GATK best practices'
#'
#' This aims to perform the following steps ( for DNA ):
#'
#' \itemize{
#' \item mark duplicates
#' \item realign indels
#' \item recalibrate bases
#' \item current version: \emph{3.4-46}
#' }
#'
#' For RNA GATK recommends a additional step of split n trim, which is not currently supported (contributions welcome !).
#'
#' \strong{NOTE}:
#'
#' Some GATK tools use \href{https://www.broadinstitute.org/gatk/guide/article?id=1975}{CPU threads while others use data threads},
#' flowr tries to use efficiently make the best use of both/either depending on tool's compatibility.
#'
#' @param bam bam file path
#' @param samplename name of the sample
#' @param outfile output file name
#' 
#' @param java_exe path to java
#' @param java_tmp path to java tmp, can leave blank
#' 
#' @param split_by_chr split processing by chromosomr where ever possible
#' 
#' @param gatk_jar_path path to gatk jar file
#' @param picard_jar_path path to picard jar file
#' @param samtools_exe path to samtools
#' 
#' @param ref_fasta_path reference fasta file
#' 
#' @param picard_markdup_opts a character vector of options for picard mark duplication step
#' @param gatk_target_opts a character vector of options for gatk target step
#' @param gatk_realign_opts a character vector of options for gatk realign step
#' @param gatk_baserecalib_opts a character vector of options for gatk baserecalib step
#' @param gatk_printreads_opts a character vector of options for gatk printreads step
#' 
#' 
#' @param mem_markdup memory used by java, example -Xmx1g
#' @param mem_target memory used by java, example -Xmx1g
#' @param mem_realign memory used by java, example -Xmx1g
#' @param mem_baserecalib memory used by java, example -Xmx1g
#' @param mem_printreads memory used by java, example -Xmx1g
#' 
#' 
#' @param cpu_markdup not used.
#' @param cpu_target number of threads used for GATK target creation step
#' @param cpu_realign number of cpu used
#' @param cpu_baserecalib number of cpu used
#' @param cpu_printreads number of cpu used
#' 
#' @param execute_cmds run commands, after creation. Useful for testing/debugging and running on local platforms.
#'
#' @export
#'
#' @examples \dontrun{
#' ## load options, including paths to tools and other parameters
#' opts_flow$load(flowr::fetch_conf("ultraseq.conf"), check = FALSE)
#' out = preprocess("my_wex.bam", samplename = "samp", split_by_chr = TRUE)
#'
#' }
preprocess <- function(bam,
                       outfile,
                       samplename = opts_flow$get("samplename"),
                       split_by_chr = opts_flow$get("split_by_chr"),

                       java_exe = opts_flow$get("java_exe"),
                       java_tmp = opts_flow$get("java_tmp"),

                       gatk_jar_path = opts_flow$get('gatk_jar_path'),
                       picard_jar_path = opts_flow$get('picard_jar_path'),
                       samtools_exe = opts_flow$get('samtools_exe'),

                       cpu_markdup = 1,
                       mem_markdup = opts_flow$get("mem_markdup"),
                       
                       cpu_target = opts_flow$get("cpu_target"),  ## not used
                       mem_target = opts_flow$get("mem_target"),

                       cpu_realign = opts_flow$get("cpu_realign"),
                       mem_realign= opts_flow$get("mem_realign"),

                       ## scatter 8 per node nct=8
                       cpu_baserecalib = opts_flow$get("cpu_baserecalib"),
                       mem_baserecalib = opts_flow$get("mem_baserecalib"),

                       ## scatter 8 per node nct=8
                       cpu_printreads = opts_flow$get("cpu_printreads"),
                       mem_printreads = opts_flow$get("mem_printreads"),

                       ref_fasta_path = opts_flow$get('ref_fasta_path'),

                       picard_markdup_opts = opts_flow$get('picard_markdup_opts'),
                       gatk_target_opts = opts_flow$get('gatk_target_opts'),
                       gatk_realign_opts = opts_flow$get('gatk_realign_opts'),
                       gatk_baserecalib_opts = opts_flow$get('gatk_baserecalib_opts'),
                       gatk_printreads_opts = opts_flow$get('gatk_printreads_opts'), 

                       execute_cmds = FALSE
                       ){

  check_args(ignore = "outfile")
  
  bamset = bam_set(bam = bam, outfile = outfile, ref_fasta_path = ref_fasta_path, split_by_chr = split_by_chr)

  # get the name of the function
  pipename = match.call()[[1]]
  message("Generating a ", pipename, " flowmat for sample: ", samplename)

  # ------------ dedup; SINGLE FILE
  dedupbam <- paste0(bamset$out_prefix, ".marked.bam")
  metricsfile <- paste0(bamset$out_prefix, ".marked.metrics")
  cmd_markdup <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s %s",
                         java_exe, mem_markdup, java_tmp, picard_jar_path, bamset$bam,
                         dedupbam, metricsfile, 
                         picard_markdup_opts)
  cmd_markdup

  ## ------------ realign; SINGLE FILE
  intervalsfiles <- paste0(bamset$out_prefix, ".realign.intervals")
  ## ------------ do this for all chrs
  cmd_target <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s %s",
                        java_exe, mem_target, java_tmp, gatk_jar_path, ref_fasta_path, dedupbam,
                        intervalsfiles, cpu_target, gatk_target_opts)

  realignedbams <- paste0(bamset$out_prefix_chr ,".realigned.bam")
  cmd_realign <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s %s %s",
                         java_exe, mem_realign, java_tmp, gatk_jar_path, ref_fasta_path, dedupbam,
                         intervalsfiles, realignedbams, gatk_realign_opts, bamset$gatk_intervals)

  ## ------------ base recalibration
  ## explicity define intervals here as well, though not required.
  recalibbams <- paste0(bamset$out_prefix_chr, ".recalibed.bam")
  recalibtabfile <- paste0(bamset$out_prefix_chr, ".recalib.tab")
  cmd_baserecalib <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T BaseRecalibrator -R %s -I %s -o %s -nct %s %s %s",
                             java_exe, mem_baserecalib, java_tmp, gatk_jar_path, ref_fasta_path,
                             realignedbams, recalibtabfile, cpu_baserecalib,
                             gatk_baserecalib_opts, bamset$gatk_intervals)
  
  cmd_printreads1 <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T PrintReads -R %s -I %s -BQSR %s -o %s -nct %s %s %s",
                            java_exe, mem_printreads, java_tmp, gatk_jar_path, ref_fasta_path, realignedbams,
                            recalibtabfile, recalibbams, cpu_printreads,
                            gatk_printreads_opts, bamset$gatk_intervals)
  cmd_printreads2 <- sprintf("%s index %s", samtools_exe, recalibbams)
  cmd_printreads = sprintf("%s;%s", cmd_printreads1, cmd_printreads2)

  cmds <- list(markdup = cmd_markdup,
               target = cmd_target, realign = cmd_realign,
               baserecalib = cmd_baserecalib, printreads = cmd_printreads)
  sapply(cmds, length)
  

  if(execute_cmds)
    sapply(cmds, system)
  
  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfiles = recalibbams))

}




