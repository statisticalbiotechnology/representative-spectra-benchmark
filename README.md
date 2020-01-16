## EuBIC Workshop on methods to merge spectra

We aim at storing code and documentation in this repository. Send an email to lukas.kall@scilifelab.se if you want to get direct write permissions.

The document is populated with pages for

- [Literature and Methods references](literature)
- [Dataset](https://github.com/ypriverol/specpride/blob/master/datasets.md)
- [Presentation](https://docs.google.com/presentation/d/1f9gMnzccAfw_EnLuwh-cbEAngYUHMVzfp19Fa_9URrc/edit?usp=sharing)

The html version of this page can be reached [here](https://statisticalbiotechnology.github.io/specpride/)

## Core contributors and Collaborators (Please add your name here):

 - [Timo Sachsenberg](sachsenb@informatik.uni-tuebingen.de), Univ. of Tübingen, Germany
 - [Eric Deutsch](edeutsch@systemsbiology.org), ISB, USA
 - [Yasset Perez-Riverol](yperez@ebi.ac.uk), EBI, UK
 - [Lukas Käll](lukas.kall@scilifelab.se), KTH, Sweden

 

### CMD execution
```
python src\convert_mgf_cluster.py mgf_add_cluster -s data\01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mgf -c data\MaRaCluster.clusters_p30.tsv -o out.mgf -a PXD004732
python src\best_spectrum.py best_spectrum -s out.mgf -o best_spectrum.mgf -m msms.txt
```