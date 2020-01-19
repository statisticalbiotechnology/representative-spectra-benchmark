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

for mgf in *.mgf
do
fn=${mgf%.mgf} 
FileConverter -in ${mgf} -out ${fn}.mzML
CometAdapter    -in ${fn}.mzML -out id/comet_${fn}.idXML -database ${FASTA} -ini id/comet.ini 
XTandemAdapter  -in ${fn}.mzML -out id/xtandem_${fn}.idXML -database ${FASTA} -ini id/xtandem.ini
MSGFPlusAdapter -in ${fn}.mzML -out id/msgf_${fn}.idXML -database ${FASTA} -ini id/MSGFPlus.ini

PSMFeatureExtractor -in id/msgf_${fn}.idXML -out id/msgf_${fn}.idXML
PSMFeatureExtractor -in id/xtandem_${fn}.idXML -out id/xtandem_${fn}.idXML
PSMFeatureExtractor -in id/comet_${fn}.idXML -out id/comet_${fn}.idXML

PeptideIndexer -in id/comet_${fn}.idXML -out id/comet_${fn}.idXML -fasta ${FASTA} -decoy_string_position "suffix"  -decoy_string "REVERSED"
PeptideIndexer -in id/xtandem_${fn}.idXML -out id/xtandem_${fn}.idXML -fasta ${FASTA} -decoy_string_position "suffix"  -decoy_string "REVERSED"
PeptideIndexer -in id/msgf_${fn}.idXML  -out id/msgf_${fn}.idXML -fasta ${FASTA} -decoy_string_position "suffix"  -decoy_string "REVERSED"

PercolatorAdapter -in id/comet_${fn}.idXML -out id/comet_${fn}.idXML  
PercolatorAdapter -in id/xtandem_${fn}.idXML -out id/xtandem_${fn}.idXML  
PercolatorAdapter -in id/msgf_${fn}.idXML -out id/msgf_${fn}.idXML  

FalseDiscoveryRate -in id/comet_${fn}.idXML -out results/comet_${fn}.idXML -FDR:PSM 0.01 -protein false 
FalseDiscoveryRate -in id/xtandem_${fn}.idXML -out results/xtandem_${fn}.idXML -FDR:PSM 0.01 -protein false
FalseDiscoveryRate -in id/msgf_${fn}.idXML -out results/msgf_${fn}.idXML -FDR:PSM 0.01 -protein false

TextExporter -in results/comet_${fn}.idXML -out results/comet_${fn}.csv -id:peptides_only -id:add_metavalues 0 -id:add_hit_metavalues 0 
TextExporter -in results/xtandem_${fn}.idXML -out results/xtandem_${fn}.csv -id:peptides_only -id:add_metavalues 0 -id:add_hit_metavalues 0 
TextExporter -in results/msgf_${fn}.idXML -out results/msgf_${fn}.csv -id:peptides_only -id:add_metavalues 0 -id:add_hit_metavalues 0 
done

wc -l results/*.csv
