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




