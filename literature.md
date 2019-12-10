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
