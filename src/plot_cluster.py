import numpy as np
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf
import sys
import pymzml


def plot_spectrum(identifier,precursor_mz,precursor_charge, mz,intensity,retention_time,peptide):
#identifier = spectrum_dict['params']['title']
#precursor_mz = spectrum_dict['params']['pepmass'][0]
#precursor_charge = spectrum_dict['params']['charge'][0]
#mz = spectrum_dict['m/z array']
#intensity = spectrum_dict['intensity array']
#retention_time = float(spectrum_dict['params']['rtinseconds'])
#peptide = 'WNQLQAFWGTGK'
    print(type(mz))
    print(type(mz[0]))
    print(np.asarray(mz, np.float32))
    print(np.asarray(intensity, np.float32))

# Create the MS/MS spectrum.
    spectrum = sus.MsmsSpectrum(
        identifier, precursor_mz, precursor_charge, mz, intensity,
        retention_time=retention_time, peptide=peptide)

    # Process the MS/MS spectrum.
    fragment_tol_mass = 10
    fragment_tol_mode = 'ppm'
    spectrum = (spectrum.set_mz_range(min_mz=100, max_mz=1400)
            .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
            .filter_intensity(min_intensity=0.05, max_num_peaks=50)
            .scale_intensity('root')
            .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                        ion_types='aby'))
    # Plot the MS/MS spectrum.
    #(sup.spectrum(spectrum).properties(width=640, height=400)
    #                   .save('spectrum_iplot.json'))


def main(mzml_file,cluster_file,scan):
    scans, tmp_scans = set(), set()
    with open(cluster_file) as cluster_def:
        for line in cluster_def:
            if line.isspace():
                if scan in tmp_scans:
                    scans.update(tmp_scans)
                tmp_scans = set()
            else:
                words = line.split('\t')
                tmp_scans.add(int(words[1]))
    print("Plotting cluster of spectra with the following scans ", scans, file=sys.stderr)
    run = pymzml.run.Reader(mzml_file)
    for n, spec in enumerate(run):
        if spec.ID in scans:
            print(
                "Spectrum {0}, MS level {ms_level} @ RT {scan_time:1.2f}".format(
                    spec.ID, ms_level=spec.ms_level, scan_time=spec.scan_time_in_minutes()
                    )
            )
            #try:
            #    spec._selected_precursors
            #except AttributeError:
            #    spec._selected_precursors = None
            precursors = spec.selected_precursors
            plot_spectrum(spec.ID,precursors[0]["mz"],precursors[0]["charge"], spec.mz,spec.i,spec.scan_time_in_minutes(),"SLEGSSDTALPLRR")
    print("Parsed {0} spectra from file {1}".format(n, mzml_file))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(main.__doc__)
        exit()
    mzml_file = sys.argv[1]
    cluster_file = sys.argv[2]
    scan = sys.argv[3]
    main(mzml_file,cluster_file,int(scan))
