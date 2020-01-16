#!/usr/bin/env python3

# Standard library
import argparse
import gzip
import json
import os
import re
import sys
import timeit
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count


# Third party
import numpy as np
from pyteomics import mzml, auxiliary


class RepresentativeSpectrumCreator:
    """Create a representative spectrum."""

    def __init__(self, verbose=None):
        # Set verbosity
        if verbose is None:
            verbose = 0
        self.verbose = verbose

    def read_cluster_list(self, file):
        """ Read MaRacluster file """

        clusters = []
        cluster = []
        icluster = 0

        with open(file) as infile:
            for line in infile:
                line = line.rstrip()
                columns = line.split()
                if len(columns) == 0:
                    clusters.append(cluster)
                    icluster += 1
                    cluster = []
                    continue
                cluster.append(columns[1])

        return clusters

    def read_spectra(self, mzml_file, scan_list):
        """ Read spectra from mzML file for given list of scan numbers. """

        # Set up information
        t0 = timeit.default_timer()
        stats = {
            "n_spectra": 0,
            "n_ms1_spectra": 0,
            "n_ms2_spectra": 0,
            "n_HCD_spectra": 0,
            "n_IT_spectra": 0,
            "n_ETD_spectra": 0,
            "high_accuracy_precursors": "unknown",
            "fragmentation_type": "unknown",
        }

        # Show information
        n_scans = len(scan_list)
        if self.verbose >= 1:
            eprint(f"INFO: Reading {n_scans} scans from mzML file {mzml_file}")
            progress_intro = False

        # Put spectra in a list
        spectra = []

        # If the mzML is gzipped, then open with zlib, else a plain open
        match = re.search("\.gz$", mzml_file)
        if match:
            infile = gzip.open(mzml_file)
        else:
            infile = open(mzml_file, "rb")

        # Read spectra from the file
        with mzml.read(infile) as reader:

            for scan in scan_list:
                spectrum = reader.get_by_id(
                    f"controllerType=0 controllerNumber=1 scan={scan}"
                )

                # Testing. Print the data structure of the first spectrum
                # if stats['n_spectra'] == 0:
                #    auxiliary.print_tree(spectrum)

                # Set a default spectrum type
                spectrum_type = "default"
                filter_string = None

                # If the ms level is 2, then examine it for information
                if spectrum["ms level"] == 2 and "m/z array" in spectrum:
                    precursor_mz = spectrum["precursorList"]["precursor"][0][
                        "selectedIonList"
                    ]["selectedIon"][0]["selected ion m/z"]
                    precursor_charge = spectrum["precursorList"]["precursor"][0][
                        "selectedIonList"
                    ]["selectedIon"][0]["charge state"]
                    print(
                        f"INFO: Reading {scan}. Precursor m/z = {precursor_mz}. n peaks={len(spectrum['m/z array'])}"
                    )
                    peaklist = {
                        "m/z array": spectrum["m/z array"],
                        "intensity array": spectrum["intensity array"],
                        "precursor mz": precursor_mz,
                        "precursor charge": precursor_charge,
                    }
                    spectra.append(peaklist)
                else:
                    print(f"ERROR: scan {scan} is not ms_level=2! Skipping")

                # Update counters and print progress
                stats["n_spectra"] += 1

        infile.close()
        if self.verbose >= 1:
            eprint("")

        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {stats['n_spectra']} spectra from {mzml_file}")
        print(f"INFO: Elapsed time: {t1-t0}")
        print(f"INFO: Processed {stats['n_spectra']/(t1-t0)} spectra per second")
        return spectra

    def read_spectra_clustered_mgf(self, clustered_mgf_file):
        """
        Read clustered MGF file and return clusters object
        
        clusters: dict of cluster_id -> peaklists
        peaklists: list of peaklist dicts
        peaklist: dict with 'm/z array', 'intensity array', 'cluster_id', 'spectrum_usi'
        precursor mz and precursor mass
        """
        all_spectra = []
        with open(clustered_mgf_file, "rt") as mgf:
            # i = 0
            for line in mgf:
                if line[:6] == "TITLE=":
                    # Initiate new spectrum
                    # i += 1
                    peaklist = {
                        "m/z array": [],
                        "intensity array": [],
                    }
                    title = line[6:].strip()
                    peaklist["cluster_id"] = title.split(";")[0]
                    peaklist["spectrum_usi"] = title.split(";")[1]
                if line[:8] == "PEPMASS=":
                    peaklist["precursor mz"] = float(line[8:].strip())
                if line[:7] == "CHARGE=":
                    peaklist["precursor charge"] = int(line[7:].strip().strip("+"))
                if line[0].isdigit():
                    peak = line.strip().split(" ")
                    peaklist["m/z array"].append(float(peak[0]))
                    peaklist["intensity array"].append(float(peak[1]))
                if line.strip() == "END IONS":
                    # Finish up this spectrum
                    all_spectra.append(peaklist)
                    # if i > 100:
                    #    break

        # Group all spectra by cluster_id
        clusters = {}
        for peaklist in all_spectra:
            if peaklist["cluster_id"] not in clusters.keys():
                clusters[peaklist["cluster_id"]] = [peaklist]
            else:
                clusters[peaklist["cluster_id"]].append(peaklist)

        return clusters

    def bin_spectra(
        self, peaklists, minimum=100, maximum=2000, binsize=0.02, peak_quorum=0.25,
        edge_case_treshold=0.5
    ):
        """
        Combine peaklists into one representative spectrum using binning.
        
        The representative spectrum will contain the mean intensities and mean m/z
        values of all peaks in a given bin. Bins with less peaks than a given peak
        quorum (percentage of total number of spectra that will be combined) will be
        ignored. Peaks in the representative spectrum that are closer together then
        the edge case threshold will get averaged into one peak.

        Positional arguments:
        peaklists -- list of dictionaries with peak lists that are to be combined

        Keyword arguments:
        minimum -- minimum m/z value in Dalton to include in the representative spectrum
        (default: 100)
        maximum -- maximum m/z value in Dalton to include in the representative spectrum
        (default: 2000)
        binsize -- size of each bin in Dalton (default: 0.02)
        peak_quorum -- minimum percentage of peaks to the total number of spectra in the
        cluster for a peak to be included in the representative spectrum (default: 0.25)
        edge_case_threshold -- fraction of the binsize, adjacent bins with average m/z
        values closer than the threshold will be merged into one bin. Set to 0 to
        disable this function (default: 0.5)
        """

        array_size = int((maximum - minimum) / binsize) + 1
        merged_spectrum = {"minimum": minimum, "maximum": maximum, "binsize": binsize}
        merged_spectrum["intensities"] = np.zeros(array_size, dtype=np.float32)
        merged_spectrum["mzs"] = np.zeros(array_size, dtype=np.float32)
        merged_spectrum["n_peaks"] = np.zeros(array_size, dtype=np.int32)
        merged_spectrum["precursor_mzs"] = []
        merged_spectrum["precursor_charges"] = []

        # Determine how many peaks need to be present to keep a final peak
        peak_quorum_int = int(len(peaklists) * peak_quorum) + 1

        for peaklist in peaklists:
            # Convert the peak lists to np arrays
            intensity_array = np.asarray(peaklist["intensity array"])
            mz_array = np.asarray(peaklist["m/z array"])

            # Limit the np arrays to the region we're interested in
            intensity_array = intensity_array[
                (mz_array >= merged_spectrum["minimum"])
                & (mz_array < merged_spectrum["maximum"])
            ]
            mz_array = mz_array[
                (mz_array >= merged_spectrum["minimum"])
                & (mz_array < merged_spectrum["maximum"])
            ]

            # Compute their bin locations and store n_peaks and intensities
            bin_array = (
                (mz_array - merged_spectrum["minimum"]) / merged_spectrum["binsize"]
            ).astype(int)

            merged_spectrum["n_peaks"][bin_array] += 1
            merged_spectrum["intensities"][bin_array] += intensity_array
            merged_spectrum["mzs"][bin_array] += mz_array

            merged_spectrum["precursor_mzs"].append(peaklist["precursor mz"])
            merged_spectrum["precursor_charges"].append(peaklist["precursor charge"])

        # Check that all precursor charges are the same
        charges = merged_spectrum["precursor_charges"]
        assert all(
            x == charges[0] for x in charges
        ), "Not all precursor charges in cluster are equal"

        # Try to handle the case where a single peak is split on a bin boundary
        # Create a temporary array of mzs that are correct means
        merged_mzs = merged_spectrum["mzs"]
        merged_mzs[merged_mzs == 0] = np.nan
        merged_mzs = np.divide(merged_mzs, merged_spectrum["n_peaks"])

        # Subtract the mzs from their previous mz
        delta_index = np.arange(0, len(merged_mzs), 1) - 1
        delta_index[0] = 0
        mz_deltas = merged_mzs - merged_mzs[delta_index]
        mz_deltas[np.isnan(mz_deltas)] = 0

        # Find cases where the deltas are smaller than half the bin size
        small_mz_deltas = mz_deltas[(mz_deltas > 0) & (mz_deltas < binsize * edge_case_treshold)]
        # if 0:
        if len(small_mz_deltas) > 0:
            # Get a list of indexes of the split bin cases
            small_mz_deltas_mask = mz_deltas
            small_mz_deltas_mask[mz_deltas == 0] = -1
            small_mz_deltas_mask[mz_deltas >= binsize / 2] = -1
            small_mz_deltas_index = delta_index[small_mz_deltas_mask > -1]

            # Consolidate all the split mzs, n_peaks, intensities into one bin
            merged_spectrum["mzs"][small_mz_deltas_index] += merged_spectrum["mzs"][
                small_mz_deltas_index + 1
            ]
            merged_spectrum["n_peaks"][small_mz_deltas_index] += merged_spectrum[
                "n_peaks"
            ][small_mz_deltas_index + 1]
            merged_spectrum["intensities"][small_mz_deltas_index] += merged_spectrum[
                "intensities"
            ][small_mz_deltas_index + 1]
            merged_spectrum["mzs"][small_mz_deltas_index + 1] = 0
            merged_spectrum["n_peaks"][small_mz_deltas_index + 1] = 0
            merged_spectrum["intensities"][small_mz_deltas_index + 1] = 0

        # Take the mean of all peaks per bin
        merged_spectrum["intensities"][
            merged_spectrum["n_peaks"] < peak_quorum_int
        ] = np.nan
        merged_spectrum["intensities"] = np.divide(
            merged_spectrum["intensities"], merged_spectrum["n_peaks"]
        )

        # Only return non-zero intensity bins
        nan_mask = ~np.isnan(merged_spectrum["intensities"])
        merged_spectrum["intensities"] = merged_spectrum["intensities"][nan_mask]

        # Compute the mean of mz values in the bin
        merged_spectrum["mzs"][merged_spectrum["mzs"] == 0] = np.nan
        merged_spectrum["mzs"] = np.divide(
            merged_spectrum["mzs"], merged_spectrum["n_peaks"]
        )
        merged_spectrum["mzs"] = merged_spectrum["mzs"][nan_mask]
        merged_spectrum["n_peaks"] = merged_spectrum["n_peaks"][nan_mask]

        merged_spectrum["precursor_mz"] = np.mean(merged_spectrum["precursor_mzs"])
        merged_spectrum["precursor_charge"] = charges[0]

        del merged_spectrum["n_peaks"]
        del merged_spectrum["precursor_charges"]
        del merged_spectrum["precursor_mzs"]

        return merged_spectrum


    def bin_spectra_wrapper(self,cluster):
        cluster_id = cluster['cluster_id']
        peaklists = cluster['peaklists']

        rsc_spectrum = self.bin_spectra(
            peaklists, minimum=100, maximum=2000, binsize=0.02,
            peak_quorum=0.25, edge_case_treshold=0
        )
        rsc_spectrum["cluster_id"] = cluster_id
        return(rsc_spectrum)


    def write_spectrum(self, spectra, mgf_file):
        """Write spectrum to MGF file."""
        for i, spectrum in enumerate(spectra):
            mgf_tmp = f"""BEGIN IONS
TITLE={spectrum["cluster_id"]}
PEPMASS={spectrum['precursor_mz']}
CHARGE={spectrum['precursor_charge']}+
"""
            for mz, intensity in zip(spectrum["mzs"], spectrum["intensities"]):
                if not np.isnan(intensity):
                    mgf_tmp += f"{mz} {intensity}\n"
            mgf_tmp += "END IONS\n\n"
            mgf_file.write(mgf_tmp)


def main():
    """Command line interface."""

    argparser = argparse.ArgumentParser(
        description="Creates an index for an MSP spectral library file"
    )
    argparser.add_argument(
        "--verbose",
        action="count",
        help="If set, print more information about ongoing processing",
    )
    argparser.add_argument("--version", action="version", version="%(prog)s 0.5")
    # argparser.add_argument('--mara_file', action='store', help='Name of the mara clusters file')
    # argparser.add_argument('--mzml_file', action='store', help='Name of the mzml file')
    # argparser.add_argument('--cluster', action='store', help='Cluster number to combine')
    argparser.add_argument(
        "--mgf_file", action="store", help="Name of the clustered MGF file"
    )
    argparser.add_argument(
        "--out",
        action="store",
        default="merged_spectra.mgf",
        help="Name of the output mgf file",
    )
    params = argparser.parse_args()

    # Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    # Print and example if not everything is provided
    # if not params.mara_file or not params.mzml_file or not params.cluster:
    #    print("Example: representative_spectrum_creator.py --mzml_file ../data/01650b_BA5-TUM_first_pool_75_01_01-3xHCD-1h-R2.mzML --mara_file=../data/MaRaCluster.clusters_p30.tsv --cluster=1")
    #    print("Or use --help for additional usage information")
    #    sys.exit(10)

    if not params.mgf_file:
        print(
            "Example: binning.py --mgf_file=../data/clustered_maracluster.mgf"
        )
        print("Or use --help for additional usage information")
        sys.exit(10)

    # Create an Representative Spectrum Creator object
    rsc = RepresentativeSpectrumCreator(verbose=verbose)

    # Read the cluster file from mzML and mara file
    # clusters = rsc.read_cluster_list(params.mara_file)
    # peaklists = rsc.read_spectra(params.mzml_file,clusters[int(params.cluster)])

    # Read the cluster file from clustered MGF
    print("Reading spectra...")
    t0 = timeit.default_timer()
    clusters = rsc.read_spectra_clustered_mgf(params.mgf_file)
    t1 = timeit.default_timer()
    print(f"INFO: Read {len(clusters)} clusters in {t1-t0} seconds ({len(clusters)/(t1-t0)} clusters per second)")

    print(f"Create representative spectra for {len(clusters)} clusters via enhanced binning method...")
    rsc_spectra = []

    n_threads = cpu_count()

    if n_threads == 1:
        for cluster_id, peaklists in clusters.items():
            rsc_spectrum = rsc.bin_spectra(
                peaklists, minimum=100, maximum=2000, binsize=0.02
            )
            rsc_spectrum["cluster_id"] = cluster_id
            rsc_spectra.append(rsc_spectrum)
    else:
        #### Reformat for parallelism
        cluster_list = []
        for cluster_id, peaklists in clusters.items():
            cluster_list.append({'cluster_id': cluster_id, 'peaklists': peaklists})

        #### Process clusters in parallel. Note will be out of order
        pool = ThreadPool(n_threads)
        rsc_spectra = pool.imap_unordered(rsc.bin_spectra_wrapper, cluster_list)
        pool.close()
        pool.join()

    t2 = timeit.default_timer()
    print(f"INFO: Processed {len(clusters)} clusters in {t2-t1} seconds ({len(clusters)/(t2-t1)} clusters per second)")

    # Write representative spectra to new MGF file
    print("Writing result MGF")
    with open(params.out, "wt") as mgf_file:
        rsc.write_spectrum(rsc_spectra, mgf_file)
    t3 = timeit.default_timer()
    print(f"INFO: Wrote {len(clusters)} representative spectra in {t3-t2} seconds ({len(clusters)/(t3-t2)} spectra per second)")


if __name__ == "__main__":
    main()
