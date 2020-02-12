#!/usr/bin/env nextflow

/**
 * Input parameters
 */
params.raw_dir = "${baseDir}/test"
params.mzml_dir  = "${baseDir}/test"

/** Create a channel for all RAW files **/
raw_files = Channel.fromPath("${params.raw_dir}/*.raw")

process msconvert {

    container 'quay.io/biocontainers/thermorawfileparser:1.2.0--0'

    memory { 10.GB * task.attempt }
    errorStrategy 'retry'
    publishDir "${params.mzml_dir}", mode:'copy', overwrite: true

    input:
    file rawFile from raw_files

    output:
    file '*.json' into metaResults mode flatten
    file '*.mzML' into spectraFiles

    script:
    """
    ThermoRawFileParser.sh -i=${rawFile} -m=0 -f=1 -o=./ --ignoreInstrumentErrors
    """
}