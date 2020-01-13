#!/bin/bash
conda create --name specpride python=3
conda install -c conda-forge -c bioconda -c defaults -y -q --name specpride --file requirements.txt
conda activate specpride

mkdir -p data
cd data
curl ftp://ftp.pride.ebi.ac.uk/pride/data/proteogenomics/projects/eubic-2020/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML
curl ftp://ftp.pride.ebi.ac.uk/pride/data/proteogenomics/projects/eubic-2020/msms.txt
curl ftp://ftp.pride.ebi.ac.uk/pride/data/proteogenomics/projects/eubic-2020/peptides.txt
cd ..
