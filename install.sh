#!/bin/bash

conda install -c conda-forge -c bioconda -c defaults -y -q --name specpride --file requirements.txt
conda activate specpride
