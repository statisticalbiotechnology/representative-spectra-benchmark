## Suggested methods to merge data
Different strategies have been suggested for forming representative spectra. [Frank et al. (JPR 2008)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2533155/) list five strategies, where one selects the representative spectrum to be:

- The "best spectrum”: the spectrum that maximizes a certain score, e.g., percent of explained intensity or percent of explained b/y ions.
- The “consensus spectrum”: a virtual spectrum constructed by averaging all spectra in the cluster. (Tabb et al. JASMS 2005)
- The “most similar spectrum”: the spectrum that has the highest average similarity to the other cluster members (Tabb et al. Anal Chem 2003).
- The “de novo spectrum”: the spectrum that has the highest score when submitted to de novo sequencing.
- The random spectrum: a spectrum chosen from the cluster at random.

## Tools that generate spectra libraries and the corresponding spectrum merge strategy

- Skyline [BiblioSpec](https://skyline.ms/announcements/home/support/thread.view?rowId=30508) - Best Spectrum in Cluster .
- spectra-cluster (PRIDE) - Consensus spectrum.
- [SpectraST](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2637392/) - Consensus cluster.
- [NIST](https://math.nist.gov/mcsd/Seminars/2016/2016-05-02-Sheetlin.html) - Consensus spectra

## Benchmarking datasets

Proteome tool's publication [Zolg *et al.*](https://www.nature.com/articles/nmeth.4153)
We have extracted a subset of this set to a [google drive](https://drive.google.com/drive/u/1/folders/1VO9VXTsfacZB7yna_3yw77a7AegRu34G).

## References
- [Frank et al. (JPR 2008)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2533155/)
- [Wang et al. (bioRxiv_2018)](https://www.biorxiv.org/content/10.1101/308627v2.full.pdf)
- [Saeed et al. (Trans Comput Biol Bioinform 2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6143137/)
- [Griss et al. (Nat Methods 2013), see supplementary material](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3667236/)
- [Perez‐Riverol et al.(Proteomics 2018)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099476/)
- [Denoising Peptide Tandem Mass Spectra for Spectral Libraries: A Bayesian Approach](https://pubs.acs.org/doi/10.1021/pr400080b)
