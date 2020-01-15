import numpy as np
import sys
import spectrum_utils.spectrum as sus
from pyteomics import mgf, parser
from scipy.stats import binned_statistic


mz_unit = 1.000508 # Space in Th between fragments
mz_space = mz_unit*.005 # Resolution approximately 0.005 Da

def bin_proc(spectrum, mz_space, max_mz):
    bins =  np.arange(-mz_space/2., max_mz, mz_space)
    # bins = np.linspace(-mz_space/2., max_mz, num=xxx)
    dig_spec, _, __ = binned_statistic(spectrum.mz, spectrum.intensity,
        statistic='sum', bins=bins)
    return dig_spec


def cos_dist(representative_spectrum, cluster_member):
    max_mz = max(representative_spectrum.mz[-1],cluster_member.mz[-1])
    discrete_a = bin_proc(representative_spectrum, mz_space, max_mz)
    discrete_b = bin_proc(cluster_member, mz_space, max_mz)
    a = np.dot(discrete_a, discrete_a)
    b = np.dot(discrete_b, discrete_b)
    ab = np.dot(discrete_a, discrete_b)
    if a==0. or b==0.:
        return 0.
    else:
        return ab/np.sqrt(a*b)


def average_cos_dist(representative_spectrum, cluster_members):
    sum_dist = 0.0
    for member in cluster_members:
        sum_dist += cos_dist(representative_spectrum,member)
    if len(cluster_members)>0:
        return sum_dist/float(len(cluster_members))
    else:
        return 0.0

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
