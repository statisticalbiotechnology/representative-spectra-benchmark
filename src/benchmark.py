import sys
import spectrum_utils.spectrum as sus
from pyteomics import mgf, parser

def fraction_of_by(peptide_seq, precursor_mz, precursor_charge,mz, intensity):
    if not parser.fast_valid(peptide_seq):
        print("Invalid peptide sequence encountered", file=sys.stderr)
        return 0.0
    spec = sus.MsmsSpectrum(
        peptide_seq, precursor_mz=precursor_mz, precursor_charge=precursor_charge,mz=mz, intensity=intensity,
        peptide=peptide_seq)
    fragment_tol_mass = 50
    fragment_tol_mode = 'ppm'
    spectrum = (spectrum.set_mz_range(min_mz=100, max_mz=1400)
            .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
            .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                        ion_types='by'))
    current, by_current = 0.,0.
    for ix in range(len(spectrum.intensity)):
        current += spectrum.intensity[ix]
        if spectrum.annotation[ix] != None:
            by_current += spectrum.intensity[ix]
    if current > 0.:
        return by_current/current
    else:
        return 0.0

if __name__ == "__main__":

    for spectrum_dict in mgf.read("data/clusters_maracluster.mgf"):
#    for spectrum_dict in mgf.read("data/.mgf"):
        peptide_seq = spectrum_dict['params']['title'].split(':')[-1][:-2]
        precursor_mz = spectrum_dict['params']['pepmass'][0]
        precursor_charge = spectrum_dict['params']['charge'][0]
        mz = spectrum_dict['m/z array']
        intensity = spectrum_dict['intensity array']
        break

    print (fraction_of_by(peptide_seq, precursor_mz, precursor_charge,mz, intensity))
