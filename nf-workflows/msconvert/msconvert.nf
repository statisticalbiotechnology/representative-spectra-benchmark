#!/usr/bin/env nextflow


/** Create a channel for all RAW files **/
raw_files = Channel.fromPath("${params.raw_folder}/*.raw")

process msconvert {

    container 'quay.io/biocontainers/thermorawfileparser:1.2.0--0'

    memory { 10.GB * task.attempt }
    errorStrategy 'retry'
    publishDir "${params.mzml_folder}", mode:'copy', overwrite: true

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