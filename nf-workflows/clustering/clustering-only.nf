#!/usr/bin/env nextflow

/**
 * Input parameters
 */
params.raw_dir = "${baseDir}/test"

/**
 * Search engine + clustering parameters
 */
// precursor tolerance can only be specified in ppm
params.prec_tol = 10
// fragment tolerance can only be specified in Th
params.frag_tol = 0.5
// missed cleavages
params.mc = 1

/**
 * Id transfer parameters
 */
params.min_ident = 2
params.min_ratio = 0.7

// number of threads per search engine
threads = 1

/**
 * Create a channel for all MGF files
 **/
(mgf_files, mgf_files_2) = Channel.fromPath("${params.raw_dir}/*.mgf").into(2)

process runClustering {
	container 'biocontainers/spectra-cluster-cli:vv1.1.2_cv2'
	publishDir "result"

	input:
	file mgf_file from mgf_files

	output:
	file "*.clustering" into spectracluster_results

	script:
	"""
	spectra-cluster-cli -major_peak_jobs ${threads} -threshold_start 1 -threshold_end 0.99 -rounds 5 -precursor_tolerance ${params.prec_tol} -precursor_tolerance_unit ppm -fragment_tolerance ${params.frag_tol} -filter mz_150 -output_path ${mgf_file}.clustering ${mgf_file}
	"""
}

process runMaRaCluster(){
   container 'ypriverol/maracluster:1.0'

   publishDir "result"

   input:
   file mgf_file from mgf_files_2

   output:
   file "maracluster_output/*.tsv" into maracluster_results

   script:
   """
   echo ${mgf_file} > bash_files.txt
   maracluster batch -b bash_files.txt -t -10
   """



}

