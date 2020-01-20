#!/usr/bin/env nextflow

/*
 *  This pipeline is designed to perform peptide identifications of mzML files.
 *
 *   Params:
 *
 *   --mzml_folder The folder containing the mzML folder
 *   --fasta      Fasta database
 *   --id_config  MSGF+ Identification step configuration file (see example configs/msgf.ini)
 *   --index_config Peptide indexer config file.
 *   --fdr_config The FDR config provides the parameters to filter peptides/proteins
 *   --idfilter_config The Identification filter config
 *   --result_folder Export to the following folder the results.
 *
 */

(mz_files, mz_result_files) = Channel.fromPath("${params.mzml_folder}/*.mzML").into(2)
fasta_file   = file(params.fasta)

/*
 * Config files for each processing step
 */

id_config    = file(params.id_config)
index_config = file(params.index_config)
fdr_config   = file(params.fdr_config)
idfilter_config = file(params.idfilter_config)

params.result_folder ='results'
result_folder = params.result_folder

/**
 * Identification step using MSGF+
 */
process peptideIdentification {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   memory { 10.GB * task.attempt }

   input:
   file mz_ml from mz_files
   file "database.fasta" from fasta_file
   file id_config

   output:
   file "*.idXML" into id_xmls

   script:
   """
   MSGFPlusAdapter -ini "${id_config}" -database database.fasta -in ${mz_ml} -out ${mz_ml.baseName}.idXML -executable /opt/thirdparty/MSGFPlus/MSGFPlus.jar
   """
}

/**
 * This function is needed to compute the PSM features
 */
process psm_features{

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   input:
   file id_xml from id_xmls

   output:
   file "*.idXML" into psm_features_xmls

   script:
   """
   PSMFeatureExtractor -in ${id_xml} -out ${id_xml.baseName}-features.idXML
   """
}

/**
 * PeptideIndexer is an step that is used to map the identified peptides to protein ids.
 */
process peptideIndexer {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   input:
   file id_xml from psm_features_xmls
   file "database.fasta" from fasta_file
   file index_config

   output:
   file "*.idXML" into index_xmls

   script:
   """
   PeptideIndexer -ini "${index_config}" -fasta database.fasta -in ${id_xml} -out ${id_xml.baseName}-index.idXML
   """
}

/**
 * Boost peptide identification using percolator tool
 */
process percolator{

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   input:
   file id_xml from index_xmls

   output:
   file "*.idXML" into percolator_xmls

   script:
   """
   PercolatorAdapter -in ${id_xml} -out ${id_xml.baseName}-percolator.idXML -percolator_executable /opt/thirdparty/Percolator/percolator
   """
}

/**
 * FalseDiscoveryRate This step compute the FDR for the identified peptides.
 */
process peptideFDRCompute {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   input:
   file index_xml from percolator_xmls
   file fdr_config

   output:
   file "*.idXML" into fdr_xmls

   script:
   """
   FalseDiscoveryRate -ini "${fdr_config}" -in ${index_xml} -out ${index_xml.baseName}-fdr.idXML
   """
}


/**
 * IDFilter This step filter the peptides using the FDR computation
 */
process peptideFDRFilter {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   input:
   file fdr_xml from fdr_xmls
   file idfilter_config

   output:
   file "*.idXML" into peptide_xmls, peptide_convert_xmls

   script:
   """
   IDFilter -ini "${idfilter_config}" -in ${fdr_xml} -out ${fdr_xml.baseName}-filter.idXML
   """
}


process convertMZIdML{

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   input:
   file filter_file from peptide_convert_xmls


   output:
   file "*.mzid" into peptide_mzids

   script:
   """
   IDFileConverter -in ${filter_file} -out ${filter_file.baseName}.mzid
   """
}

(mzMLs, mzML_print) = mz_result_files.map { file -> tuple(file.name, file)}.into(2)
(peptide_collection, peptide_print) = peptide_xmls.map { file -> tuple(file.baseName.replaceFirst("-index-fdr-filter",""), file)}.into(2)
(combined_results, print_combined) = mzMLs.combine(peptide_collection, by: 0).into(2)


mzML_print.subscribe{ println "value: $it"}
peptide_print.subscribe{ println "value: $it"}
print_combined.subscribe{ println "value: $it"}

process idQualityControl{
   container 'mwalzer/openms-batteries-included:V2.3.0_pepxmlpatch'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   input:
   set val(mzML), file(mzML_file), file(idx_file) from combined_results

   output:
   file "*.qcML" into qc_MLs

   script:
   """
   QCCalculator -in ${mzML_file} -id ${idx_file} -out ${mzML}.qcML
   """
}

qc_tool_parameters = ['fracmass', 'auctic', 'charge_histogram', 'dppm', 'dppm_time', 'dppm_percentiles', 'esiinstability', 'gravy',
                                  'idmap', 'id_oversampling', 'lengthdistro', 'ms1peakcount', 'ms2peakcount', 'repeatid',
                                  'rt_events',  'tic', 'ticric', 'topn']

/**
 * 'ms1sn', 'ms2sn', 'sn'
 */

process plotQualityControl{

  container 'mwalzer/qc-plotter:latest'
  publishDir "${params.result_folder}", mode: 'copy', overwrite: true

  input:
  file qc_ML from qc_MLs
  each param from qc_tool_parameters


  output:
  file "*.png" into plots

  script:
  """
  qc_plot.sh -$param ${qc_ML}
  """

}