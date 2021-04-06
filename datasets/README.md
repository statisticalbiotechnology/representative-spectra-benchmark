# Benchmarking datasets

Here we will describe the datasets use to benchmark the different methods:

 - Synthetic peptides of [proteome tools]( http://www.proteometools.org/index.php?id=52) , i.e. [PXD004732](https://www.ebi.ac.uk/pride/archive/projects/PXD004732).

## Synthetic peptides (PXD004732)

We selected the run described by [TUM_first_pool_75_01_01_3xHCD-1h-R2-tryptic.zip](https://drive.google.com/open?id=1nDF2yOY2JU0UNml-MEMVGZHOoZD9uAB_), and its associated spectra [01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.raw](https://drive.google.com/open?id=1P0TWB9O0PzVCB1_1m3gxq93T2ZsYrul1).

 - [MzML](https://drive.google.com/open?id=1CeAFcPZNzjHd7tqAntXupwKt0n5E2fza)
 - [MGF](https://drive.google.com/open?id=1nRGllZeNmHupXnIaPoy27Oz1NLIEjw4F)

Some of the original parameters for the search are the following:

  - Oxidation M (Variable Modification)
  - Carbamidomethyl (C) (Fixed Modification)
  - Precursor mass Tolerance (10 ppm)
  - Fragment mass Tolerance (20.0 ppm)
  - Enzyme (Trypsin/P), fully , no - miss-cleavages

The Fasta database used to perform the search can be found [here](https://raw.githubusercontent.com/statisticalbiotechnology/specpride/master/datasets/PXD004732/final_concatenated_target_decoy.fasta). For convenience, we made the data available through a [google drive](https://drive.google.com/open?id=1UkI6Uvuo9AimRrGJGMjfLWSMZgFoDm9k).

### Running analysis in nextflow

The folder `configs` contains the configuration files for the identification workflow in nextflow.


In order to run the identification pipeline with this dataset please use the following command:

```bash
nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder data/PXD004732/ --fasta data/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/ -resume
```

Note, that you need to be in the `specpride` main folder and the mzML of the MSrun should be in the following folder `data/PXD004732/` and the fasta file should be in `data/PXD004732/final_concatenated_target_decoy.fasta`.



~ Addition(XiYang-Luo 20210204)

### The following headlines indicate the process, and the following processes all take 01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2 in the PXD004732 project as an example

## 1. MSconvert:

```bash
#Convert original files into mzML files and mgf files
activating environment:	conda activate nf-core-proteomicstmt-1.0dev
nextflow run nf-workflows/msconvert/msconvert.nf --raw_folder datasets/PXD004732/
```

## 2. identification original spectra

PXD004732:

// adding profile:	LuciphorAdapter -write_ini phospho_config.ini

The reason for adding this configuration file is that I reported an error when running the command. And some missing parameters are also added in the main.nf file and the corresponding configuration file. I will place the modified main.nf file in the compressed package

```bash
activating environment:	conda activate nf-core-proteomicstmt-1.0dev

nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder datasets/PXD004732/ --fasta datasets/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/ -resume PXD004732 --phospho false --phospho_config datasets/PXD004732/configs/phospho_config.ini
```

## 3. Clustering

// Here I decompose maracluster and pride into two nextflow files. The reason is that processing PRIDE in the original clustering-only.nf file requires input of MGF format files, not mzML files. So I generated two mextflow files for the two methods.

### // It should also be noted here that after converting the original file (.RAW) to mgf, the spectra will be lost. In other words, the number of spectra in the mzML file converted from the .raw file is always greater than that in the mgf format file. There are more spectra (of course, I don't know why this is the case).

```bash
activating environment:	conda activate nf-core-proteomicstmt-1.0dev
#spectra-cluster-cli:
nextflow run nf-workflows/clustering/t1.nf --raw_dir datasets/PXD004732/
#maracluster
nextflow run nf-workflows/clustering/t2.nf --raw_dir datasets/PXD004732/
```

## 4. Spectra_add_cluster:

activate environment (When running the python script):

```bash
conda activate specpride

python spectra_add_cluster.py --spectra /mnt/h/xiyang/specpride-master/datasets/PXD004732/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML --cluster /mnt/h/xiyang/specpride-master/datasets/PXD004732/spectra-cluster_output/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mgf.clustering 'spectra-cluster' --out /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mgf/S/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addSpectra-cluster.mgf
```

```bash
conda activate specpride

python spectra_add_cluster.py --spectra /mnt/h/xiyang/specpride-master/datasets/PXD004732/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML --cluster /mnt/h/xiyang/specpride-master/datasets/PXD004732/maracluster_output/MaRaCluster.clusters_p30.tsv 'MaRaCluster' --out /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mgf/M/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addMaRaCluster.mgf
```

## 5. Generating consensus

### BEST

// It should be noted that in the ms_io.py file, when the clustering file generated by the pride is read, which is the function: def _read_clusters_spectracluster(), the original code downloaded from GitHub is read in index instead of scan, which will cause later An error is reported. The reason is that the index cannot correspond to the scan in the original spectrum (this is obviously, the scan can correspond to the scan). So here is modified. The corresponding code will also be attached to the compressed package.

```bash
python representative.py --filename_in /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mzml/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addSpectra-cluster.mzML --filename_out /mnt/h/xiyang/specpride-master/datasets/PXD004732/consensus_output/spectracluster/mgf/best/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_Spectra-clusterConsensus-best.mgf --representative_method "best_spectrum" --filename_psm /mnt/h/xiyang/specpride-master/datasets/PXD004732/msms.txt  --min_cluster_size 1
```

```bash
python representative.py --filename_in /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mzml/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addMaRaCluster.mzML --filename_out /mnt/h/xiyang/specpride-master/datasets/PXD004732/consensus_output/maracluster/mgf/best/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_MaRaClusterConsensus-best.mgf --representative_method "best_spectrum" --filename_psm /mnt/h/xiyang/specpride-master/datasets/PXD004732/msms.txt  --min_cluster_size 1
```

### MOST

consensus mgf  most_similar

```bash
python representative.py --filename_in /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mzml/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addSpectra-cluster.mzML --filename_out /mnt/h/xiyang/specpride-master/datasets/PXD004732/consensus_output/spectracluster/mgf/most/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_Spectra-clusterConsensus-most.mgf --representative_method "most_similar" --filename_psm /mnt/h/xiyang/specpride-master/datasets/PXD004732/msms.txt --min_cluster_size 1
```

```bash
python representative.py --filename_in /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mzml/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addMaRaCluster.mzML --filename_out /mnt/h/xiyang/specpride-master/datasets/PXD004732/consensus_output/maracluster/mgf/most/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_MaRaClusterConsensus-most.mgf --representative_method "most_similar" --filename_psm /mnt/h/xiyang/specpride-master/datasets/PXD004732/msms.txt --min_cluster_size 1
```

### BIN

consensus mgf  bin

```bash
python representative.py --filename_in /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mzml/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addSpectra-cluster.mzML --filename_out /mnt/h/xiyang/specpride-master/datasets/PXD004732/consensus_output/spectracluster/mgf/bin/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_Spectra-clusterConsensus-bin.mgf --representative_method "bin" --filename_psm /mnt/h/xiyang/specpride-master/datasets/PXD004732/msms.txt --min_cluster_size 1
```

```bash
python representative.py --filename_in /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mzml/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addMaRaCluster.mzML --filename_out /mnt/h/xiyang/specpride-master/datasets/PXD004732/consensus_output/maracluster/mgf/bin/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_MaRaClusterConsensus-bin.mgf --representative_method "bin" --filename_psm /mnt/h/xiyang/specpride-master/datasets/PXD004732/msms.txt --min_cluster_size 1
```



### AVERAGE

// It should be noted that the average_spectrum_clustering.py file downloaded on GitHub did not merge all the spectra in each cluster into one consensus (I don't know if there is a problem with the script here). I later modified it. Correspondingly The file is attached to the compressed package

```bash
python average_spectrum_clustering.py   /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mgf/M/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addMaRaCluster.mgf  /mnt/h/xiyang/specpride-master/datasets/PXD004732/consensus_output/maracluster/mgf/average/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_MaRaClusterConsensus-average.mgf
```

```bash
python average_spectrum_clustering.py   /mnt/h/xiyang/specpride-master/datasets/PXD004732/addCluster_output/mgf/S/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_addSpectra-cluster.mgf  /mnt/h/xiyang/specpride-master/datasets/PXD004732/consensus_output/spectracluster/mgf/average/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2_Spectra-clusterConsensus-average.mgf
```



## 6. Identification consensus spectra:

cd /mnt/h/xiyang/specpride-master

activating environment:	conda activate nf-core-proteomicstmt-1.0dev

Note: The generated mgf needs to be converted to mzml using MSconvert.

It should be noted here that some mgf files can pass through the program after being converted by MSConvert, but some cannot (very weird).

### best

```bash
#MARAcluster: best_spectrum
nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder datasets/PXD004732/consensus_output/maracluster1/mzmlC-best/ --fasta datasets/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/consensus_output/maracluster1/mzmlC-best/ -resume PXD004732 --phospho false --phospho_config datasets/PXD004732/configs/phospho_config.ini
```

```bash
#Spectra-cluster: best_spectrum
nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder datasets/PXD004732/consensus_output/spectracluster1/mzmlC-best/ --fasta datasets/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/consensus_output/spectracluster1/mzmlC-best/ -resume PXD004732 --phospho false --phospho_config datasets/PXD004732/configs/phospho_config.ini
```

### most

```bash
#MARAcluster: most_similar
nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder datasets/PXD004732/consensus_output/maracluster1/mzmlC-most/ --fasta datasets/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/consensus_output/maracluster1/mzmlC-most/ -resume PXD004732M --phospho false --phospho_config datasets/PXD004732/configs/phospho_config.ini
```

```bash
#Spectra-cluster: most_similar
nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder datasets/PXD004732/consensus_output/spectracluster1/mzmlC-most/ --fasta datasets/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/consensus_output/spectracluster1/mzmlC-most/ -resume PXD004732S --phospho false --phospho_config datasets/PXD004732/configs/phospho_config.ini
```

### bin

```bash
#MARAcluster: bin
nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder datasets/PXD004732/consensus_output/maracluster1/mzmlC-bin/ --fasta datasets/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/consensus_output/maracluster1/mzmlC-bin/ -resume PXD004732M --phospho false --phospho_config datasets/PXD004732/configs/phospho_config.ini
```

```bash
#Spectra-cluster: bin
nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder datasets/PXD004732/consensus_output/spectracluster1/mzmlC-bin/ --fasta datasets/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/consensus_output/spectracluster1/mzmlC-bin/ -resume PXD004732S --phospho false --phospho_config datasets/PXD004732/configs/phospho_config.ini
```

### average

```bash
#MARAcluster: average
nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder datasets/PXD004732/consensus_output/maracluster1/mzmlC-average/ --fasta datasets/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/consensus_output/maracluster1/mzmlC-average/ -resume PXD004732M --phospho false --phospho_config datasets/PXD004732/configs/phospho_config.ini
```

```bash
#Spectra-cluster: average
nextflow run nf-workflows/identification/main.nf -c nf-workflows/identification/nextflow.config -profile local,trace  --mzml_folder datasets/PXD004732/consensus_output/spectracluster1/mzmlC-average/ --fasta datasets/PXD004732/final_concatenated_target_decoy.fasta --id_config datasets/PXD004732/configs/msgf.ini --index_config datasets/PXD004732/configs/pi.ini --fdr_config datasets/PXD004732/configs/fdr.ini --idfilter_config datasets/PXD004732/configs/idf.ini --result_folder datasets/PXD004732/consensus_output/spectracluster1/mzmlC-average/ -resume PXD004732S --phospho false --phospho_config datasets/PXD004732/configs/phospho_config.ini
```







































conda activate my-rdkit-env





