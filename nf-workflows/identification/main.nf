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
 *   --phospho (true | false) Perform Phospho Localization
 *   --phospho_config Provide a Phospho localization configuration file for LuciphorAdapter
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
phospho_config  = file(params.phospho_config)

params.phospho = false


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
 * Merge Identification files from multiple RAW files.
 */
process idMerger {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   memory { 16.GB * task.attempt }
   errorStrategy 'retry'

   input:
   file(input_files) from id_xmls.collect()

   output:
   file "*.idXML" into merged_xmls

   script:
   """
   IDMerger -in ${(input_files as List).join(" ")} -out merged-global-index.idXML
   """

}

/**
 * This function is needed to compute the PSM features
 */
process psm_features{

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   memory { 16.GB * task.attempt }
   errorStrategy 'retry'

   input:
   file id_xml from merged_xmls

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

   memory { 16.GB * task.attempt }
   errorStrategy 'retry'

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

   memory { 16.GB * task.attempt }
   errorStrategy 'retry'

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

   memory { 16.GB * task.attempt }
   errorStrategy 'retry'

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

   memory { 16.GB * task.attempt }
   errorStrategy 'retry'

   input:
   file fdr_xml from fdr_xmls
   file idfilter_config

   output:
   file "*.idXML" into peptide_xmls, peptide_convert_xmls, phospho_xmls

   script:
   """
   IDFilter -ini "${idfilter_config}" -in ${fdr_xml} -out ${fdr_xml.baseName}-filter.idXML
   """
}


/**
 * IDFilter This step filter the peptides using the FDR computation
 */
process phosphoLocalization {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   when:
    params.phospho

   memory { 16.GB * task.attempt }
   errorStrategy 'retry'

   input:
   file fdr_xml from phospho_xmls
   file phospho_config

   output:
   file "*.idXML" into convert_phospho_xmls

   script:
   """
   LuciphorAdapter -ini "${phospho_config}" -in ${fdr_xml} -out ${fdr_xml.baseName}-phospho.idXML
   """
}

merged_files = peptide_convert_xmls.mix(convert_phospho_xmls)

/**
 * Convert a set of idXML files to mzIdentML
 */
process convertMZIdML{

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq_2'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   memory { 16.GB * task.attempt }
   errorStrategy 'retry'

   input:
   file filter_file from merged_files

   output:
   file "*.mzid" into peptide_mzids

   script:
   """
   IDFileConverter -in ${filter_file} -out ${filter_file.baseName}.mzid
   """
}