% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{preprocess}
\alias{preprocess}
\title{Pre-process bam files following Broad's best practices for variant calling, starting from aligned BAM file}
\usage{
preprocess(bam, outfile, samplename = opts_flow$get("samplename"),
  split_by_chr = opts_flow$get("split_by_chr"),
  java_exe = opts_flow$get("java_exe"),
  java_tmp = opts_flow$get("java_tmp"),
  gatk_jar_path = opts_flow$get("gatk_jar_path"),
  picard_jar_path = opts_flow$get("picard_jar_path"),
  samtools_exe = opts_flow$get("samtools_exe"), cpu_markdup = 1,
  mem_markdup = opts_flow$get("mem_markdup"),
  cpu_target = opts_flow$get("cpu_target"),
  mem_target = opts_flow$get("mem_target"),
  cpu_realign = opts_flow$get("cpu_realign"),
  mem_realign = opts_flow$get("mem_realign"),
  cpu_baserecalib = opts_flow$get("cpu_baserecalib"),
  mem_baserecalib = opts_flow$get("mem_baserecalib"),
  cpu_printreads = opts_flow$get("cpu_printreads"),
  mem_printreads = opts_flow$get("mem_printreads"),
  ref_fasta_path = opts_flow$get("ref_fasta_path"),
  picard_markdup_opts = opts_flow$get("picard_markdup_opts"),
  gatk_target_opts = opts_flow$get("gatk_target_opts"),
  gatk_realign_opts = opts_flow$get("gatk_realign_opts"),
  gatk_baserecalib_opts = opts_flow$get("gatk_baserecalib_opts"),
  gatk_printreads_opts = opts_flow$get("gatk_printreads_opts"),
  execute_cmds = FALSE)
}
\arguments{
\item{bam}{bam file path}

\item{outfile}{output file name}

\item{samplename}{name of the sample}

\item{split_by_chr}{split processing by chromosomr where ever possible}

\item{java_exe}{path to java}

\item{java_tmp}{path to java tmp, can leave blank}

\item{gatk_jar_path}{path to gatk jar file}

\item{picard_jar_path}{path to picard jar file}

\item{samtools_exe}{path to samtools}

\item{cpu_markdup}{not used.}

\item{mem_markdup}{memory used by java, example -Xmx1g}

\item{cpu_target}{number of threads used for GATK target creation step}

\item{mem_target}{memory used by java, example -Xmx1g}

\item{cpu_realign}{number of cpu used}

\item{mem_realign}{memory used by java, example -Xmx1g}

\item{cpu_baserecalib}{number of cpu used}

\item{mem_baserecalib}{memory used by java, example -Xmx1g}

\item{cpu_printreads}{number of cpu used}

\item{mem_printreads}{memory used by java, example -Xmx1g}

\item{ref_fasta_path}{reference fasta file}

\item{picard_markdup_opts}{a character vector of options for picard mark duplication step}

\item{gatk_target_opts}{a character vector of options for gatk target step}

\item{gatk_realign_opts}{a character vector of options for gatk realign step}

\item{gatk_baserecalib_opts}{a character vector of options for gatk baserecalib step}

\item{gatk_printreads_opts}{a character vector of options for gatk printreads step}

\item{execute_cmds}{run commands, after creation. Useful for testing/debugging and running on local platforms.}
}
\description{
This function provides a wrapper around the best practices described on \href{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}{GATK's website}.
If the link is broken google 'GATK best practices'

This aims to perform the following steps ( for DNA ):

\itemize{
\item mark duplicates
\item realign indels
\item recalibrate bases
\item current version: \emph{3.4-46}
}

For RNA GATK recommends a additional step of split n trim, which is not currently supported (contributions welcome !).

\strong{NOTE}:

Some GATK tools use \href{https://www.broadinstitute.org/gatk/guide/article?id=1975}{CPU threads while others use data threads},
flowr tries to use efficiently make the best use of both/either depending on tool's compatibility.
}
\details{
Flow following Broad's best practices for variant calling, starting from sorted bam
}
\examples{
\dontrun{
## load options, including paths to tools and other parameters
opts_flow$load(flowr::fetch_conf("ultraseq.conf"), check = FALSE)
out = preprocess("my_wex.bam", samplename = "samp", split_by_chr = TRUE)

}
}

