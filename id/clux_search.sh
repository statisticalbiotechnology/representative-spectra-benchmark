#!/bin/bash
mkdir crux
cat peptides.txt | cut -f 1 | tail -n +2 | gawk '{print ">" $0; print $0}' > crux/pept.fa
cd crux
crux tide-index --mods-spec 3M+15.9949 pept.fa pept.idx
crux tide-search ../01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML pept.idx
crux percolator --overwrite T crux-output/tide-search.target.txt crux-output/tide-search.decoy.txt
