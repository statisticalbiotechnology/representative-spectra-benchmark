import numpy as np
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf
import sys
import argparse


def main(cluster_file, consensus_file):
    with mgf.read(consensus_file) as reader:
        for spectrum_dict in reader:
            peptide_seq = spectrum_dict['params']['title']
            precursor_mz = spectrum_dict['params']['pepmass'][0]
            precursor_charge = spectrum_dict['params']['charge'][0]
            cons_mz = spectrum_dict['m/z array']
            cons_intensity = spectrum_dict['intensity array']
            retention_time = float(spectrum_dict['params']['rtinseconds'])
            break
    cons_spec = sus.MsmsSpectrum(
        peptide_seq, precursor_mz=precursor_mz, precursor_charge=precursor_charge, mz=cons_mz, intensity=cons_intensity,
        retention_time=retention_time, peptide=peptide_seq)
    with mgf.read(cluster_file) as reader:
        for spectrum_dict in reader:
            precursor_mz = spectrum_dict['params']['pepmass'][0]
            precursor_charge = spectrum_dict['params']['charge'][0]
            mz = spectrum_dict['m/z array']
            intensity = spectrum_dict['intensity array']
            retention_time = float(spectrum_dict['params']['rtinseconds'])

        spectrum = sus.MsmsSpectrum(peptide_seq, precursor_mz=precursor_mz, precursor_charge=precursor_charge, mz=mz,
                                    intensity=intensity,
                                    retention_time=retention_time, peptide=peptide_seq)
        # Process the MS/MS spectrum.
        fragment_tol_mass = 10
        fragment_tol_mode = 'ppm'
        #    fragment_tol_mass = .5
        #    fragment_tol_mode = 'Da'
        spectrum = (spectrum.set_mz_range(min_mz=100, max_mz=1400)
                    .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
                    .filter_intensity(min_intensity=0.05, max_num_peaks=50)
                    .scale_intensity('root')
                    .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                                ion_types='aby'))
    # Generate theoretical spec
    # Plot the MS/MS spectrum.
    fig, ax = plt.subplots(figsize=(12, 6))
    #    sup.spectrum(spectrum, ax=ax)
    sup.mirror(spectrum, tspec, ax=ax)
    plt.show()
    plt.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(prog='plot_cluster_vs_consensus', usage='%(prog)s [options]',
                                 description="Program plotting mirror plots of cluster members an their representative spectrum")

    # Add the arguments to the parser
    ap.add_argument("cluster_file", nargs=1, help='The mgf file defining the cluster members')
    ap.add_argument("consensus_file", nargs=1, help='The mgf file defining the representative spectrum')
    args = vars(ap.parse_args())

    main(args["cluster_file"][0], args["consensus_file"][0])
