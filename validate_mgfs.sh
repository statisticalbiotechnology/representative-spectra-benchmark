#!/bin/bash

# uses OpenMS 2.5 pre-release from github but should work with 2.4 as well
# make sure that third-party executables are available (can be checkout out from OpenMS/THIRDPARTY)
export PATH="$PATH:/home/sachsenb/OMS/OpenMS-build/bin"
export PATH="$PATH:/home/sachsenb/OpenMS/THIRDPARTY/Linux/64bit/Comet"
export PATH="$PATH:/home/sachsenb/OpenMS/THIRDPARTY/Linux/64bit/XTandem"
export PATH="$PATH:/home/sachsenb/OpenMS/THIRDPARTY/All/MSGFPlus"
export PATH="$PATH:/home/sachsenb/OpenMS/THIRDPARTY/Linux/64bit/Percolator"
export PATH="$PATH:/home/sachsenb/OpenMS/THIRDPARTY/Linux/64bit/MaRaCluster"
export PATH="$PATH:/home/sachsenb/OpenMS/THIRDPARTY/All/ThermoRawFileParser"

FASTA="id/final_concatenated_target_decoy.fasta"

mgf=$1
fn=${mgf%.mgf} 

FileConverter -in ${mgf} -out ${fn}.mzML
CometAdapter  -ini id/comet.ini -in ${fn}.mzML -out id/comet_${fn}.idXML -database ${FASTA}
XTandemAdapter -in ${fn}.mzML -out id/xtandem_${fn}.idXML -database ${FASTA} -ini id/xtandem.ini
MSGFPlusAdapter -in ${fn}.mzML -out id/msgf_${fn}.idXML -database ${FASTA} -ini id/MSGFPlus.ini

PeptideIndexer -in id/comet_${fn}.idXML -out id/comet_${fn}.idXML -fasta ${FASTA} -decoy_string_position "suffix"  -decoy_string "REVERSED"
PeptideIndexer -in id/xtandem_${fn}.idXML -out id/xtandem_${fn}.idXML -fasta ${FASTA} -decoy_string_position "suffix"  -decoy_string "REVERSED"
PeptideIndexer -in id/msgf_${fn}.idXML  -out id/msgf_${fn}.idXML -fasta ${FASTA} -decoy_string_position "suffix"  -decoy_string "REVERSED"

FalseDiscoveryRate -in id/comet_${fn}.idXML -out results/comet_${fn}.idXML -FDR:PSM 0.01 -protein false 
FalseDiscoveryRate -in id/xtandem_${fn}.idXML -out results/xtandem_${fn}.idXML -FDR:PSM 0.01 -protein false
FalseDiscoveryRate -in id/msgf_${fn}.idXML -out results/msgf_${fn}.idXML -FDR:PSM 0.01 -protein false


