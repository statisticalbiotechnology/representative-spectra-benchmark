import numpy as np
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf
import sys
import pymzml


def plot_spectrum(identifier, precursor_mz, precursor_charge, mz, intensity, retention_time, peptide):
    # identifier = spectrum_dict['params']['title']
    # precursor_mz = spectrum_dict['params']['pepmass'][0]
    # precursor_charge = spectrum_dict['params']['charge'][0]
    # mz = spectrum_dict['m/z array']
    # intensity = spectrum_dict['intensity array']
    # retention_time = float(spectrum_dict['params']['rtinseconds'])
    # peptide = 'WNQLQAFWGTGK'

    # Create the MS/MS spectrum.
    spectrum = sus.MsmsSpectrum(
        identifier, precursor_mz=precursor_mz, precursor_charge=precursor_charge, mz=mz, intensity=intensity,
        retention_time=retention_time, peptide=peptide)

    # Process the MS/MS spectrum.
    #    fragment_tol_mass = 10
    #    fragment_tol_mode = 'ppm'
    fragment_tol_mass = .5
    fragment_tol_mode = 'Da'
    spectrum = (spectrum.set_mz_range(min_mz=100, max_mz=1400)
                .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
                .filter_intensity(min_intensity=0.05, max_num_peaks=50)
                .scale_intensity('root')
                .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                            ion_types='aby'))
    # Generate theoretical spec
    ts = sus._get_theoretical_peptide_fragments(peptide)
    tmz = [frag.calc_mz for frag in ts]
    ti = [1.0 for _ in tmz]
    tspec = sus.MsmsSpectrum(
        identifier, precursor_mz=precursor_mz, precursor_charge=precursor_charge, mz=tmz, intensity=ti,
        retention_time=retention_time, peptide=peptide)
    # Plot the MS/MS spectrum.
    fig, ax = plt.subplots(figsize=(12, 6))
    #    sup.spectrum(spectrum, ax=ax)
    sup.mirror(spectrum, tspec, ax=ax)
    plt.show()
    plt.close()


def main(mzml_file, cluster_file, msms_file, scan):
    scans, tmp_scans, peptide, spectra = set(), set(), "", []
    with open(cluster_file) as cluster_def:
        for line in cluster_def:
            if line.isspace():
                if scan in tmp_scans:
                    scans.update(tmp_scans)
                tmp_scans = set()
            else:
                words = line.split('\t')
                tmp_scans.add(int(words[1]))
    with open(msms_file) as peptides:
        next(peptides)  # skip header
        for line in peptides:
            words = line.split('\t')
            rscan = int(words[1])
            rpept = words[7][1:-1]
            if scan == rscan:
                peptide = rpept
                break
    print("Plotting cluster of spectra with the following scans ", scans, " for sequence ", peptide, file=sys.stderr)
    run = pymzml.run.Reader(mzml_file)
    for n, spec in enumerate(run):
        if spec.ID in scans:
            print(
                "Spectrum {0}, MS level {ms_level} @ RT {scan_time:1.2f}".format(
                    spec.ID, ms_level=spec.ms_level, scan_time=spec.scan_time_in_minutes()
                )
            )
            spectra.append(spec)
            # try:
            #    spec._selected_precursors
            # except AttributeError:
            #    spec._selected_precursors = None
            precursors = spec.selected_precursors
            print(precursors[0]["mz"], precursors[0]["charge"])
            plot_spectrum(spec.ID, precursors[0]["mz"], precursors[0]["charge"], np.array(spec.mz), np.array(spec.i),
                          spec.scan_time_in_minutes(), peptide)
    print("Parsed {0} spectra from file {1}".format(n, mzml_file))


#    average_merge(spectra)

# def average_merge(spectra):
#

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(main.__doc__)
        exit()
    mzml_file = sys.argv[1]
    cluster_file = sys.argv[2]
    peptide_file = sys.argv[3]
    scan = sys.argv[4]
    main(mzml_file, cluster_file, peptide_file, int(scan))
