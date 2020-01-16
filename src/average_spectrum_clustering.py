import argparse
import numpy as np
from itertools import groupby
from pyteomics import mgf, mass

H = mass.nist_mass['H+'][0][0]

# The following constants represent arbitrarily set default values for several parameters
# Each of these parameters will affect benchmarking and should be carefully chosen.
# See also the list of optional command-line arguments, as they will affect output.

# Here is an incomplete list of things that can be implemented differently:
# - how to set the precursor mass for representative spectrum
# - whether to filter the representative spectrum for low-intensity peaks and how
#   (e.g. "dynamic range") filter currently implemented
# - whether to use peak intensities as weights for averaging within MS/MS-level clusters
# - filter peak clusters based on the number of scans they are representing
# - spectrum preprocessing (deisotoping / denoising / "demixing")
#   before clustering AND in subsequent benchmarks

DIFF_THRESH = 0.01
DYN_RANGE = 1000
MIN_FRACTION = 0.5


def average_spectrum(spectra, title='', pepmass='', rtinseconds='', charge='', **kwargs):
    '''
    Produces an average spectrum for a cluster.

    Parameters
    ----------
    spectra : Iterable
        An iterable of spectra, each spectrum is expected to be in pyteomics format.
    title : str, optional
        Title of the output spectrum
    mz_accuracy : float, keyword only, optional
        m/z accuracy used for MS/MS clustering. Default is `DIFF_THRESH`
    dyn_range : float, keyword only, optional
        Dynamic range to apply to output (peaks less than max_intensity / dyn_range
        are discarded)
    msms_avg : {'naive', 'weighted'}
        Method for calculation of MS/MS peak m/z values in representative spectrum.
        Naive: simple average of peak m/z within MS/MS-level cluster.
        Weighted: weighted average with MS/MS peak intensities as weights.
    min_fraction, float, keyword only, optional
        Minimum fraction of cluster spectra need to contain the peak.

    Returns
    -------
    out : average spectrum in pyteomics format
    '''
    mz_accuracy = kwargs.get('mz_accuracy', DIFF_THRESH)
    dyn_range = kwargs.get('dyn_range', DYN_RANGE)
    min_fraction = kwargs.get('min_fraction', MIN_FRACTION)
    msms_avg = kwargs.get('msms_avg')

    mz_arrays, int_arrays = [], []
    n = 0 # number of spectra
    for s in spectra:
        mz_arrays.append(s['m/z array'])
        int_arrays.append(s['intensity array'])
        n += 1
    if n > 1:
        mz_array_all = np.concatenate(mz_arrays)
        intensity_array_all = np.concatenate(int_arrays)

        idx = np.argsort(mz_array_all)
        mz_array_all = mz_array_all[idx]
        intensity_array_all = intensity_array_all[idx]
        diff_array = np.diff(mz_array_all)

        new_mz_array = []
        new_intensity_array = []

        ind_list = list(np.where(diff_array >= mz_accuracy)[0]+1)

        i_prev = ind_list[0]

        mz_array_sum = np.cumsum(mz_array_all)
        intensity_array_sum = np.cumsum(intensity_array_all)
        if msms_avg == 'weighted':
            mult_mz_intensity_array_sum = np.cumsum(mz_array_all*intensity_array_all)

        min_l = min_fraction * n
        if i_prev >= min_l:
            I_sum = intensity_array_sum[i_prev-1]
            if msms_avg == 'naive':
                new_mz_array.append(mz_array_sum[i_prev-1] / i_prev)
            elif msms_avg == 'weighted':
                mzI_sum = mult_mz_intensity_array_sum[i_prev-1]
                new_mz_array.append(mzI_sum / I_sum)
            new_intensity_array.append(I_sum / n)

        for i in ind_list[1:-1]:
            if i - i_prev >= min_l:
                I_sum = intensity_array_sum[i-1] - intensity_array_sum[i_prev-1]
                if msms_avg == 'naive':
                    mz_sum = mz_array_sum[i-1] - mz_array_sum[i_prev-1]
                    new_mz_array.append(mz_sum / (i-i_prev))
                elif msms_avg == 'weighted':
                    mzI_sum = mult_mz_intensity_array_sum[i-1] - mult_mz_intensity_array_sum[i_prev-1]
                    new_mz_array.append(mzI_sum / I_sum)
                new_intensity_array.append(I_sum / n)
            i_prev = i

        if (len(mz_array_sum) - i_prev) >= min_l:
            I_sum = intensity_array_sum[-1] - intensity_array_sum[i_prev-1]
            if msms_avg == 'naive':
                mz_sum = mz_array_sum[-1] - mz_array_sum[i_prev-1]
                new_mz_array.append(mz_sum / (len(mz_array_sum) - i_prev))
            elif msms_avg == 'weighted':
                mzI_sum = mult_mz_intensity_array_sum[-1] - mult_mz_intensity_array_sum[i_prev-1]
                new_mz_array.append(mzI_sum / I_sum)
            new_intensity_array.append(I_sum / n)
    else:
        new_mz_array = mz_arrays[0]
        new_intensity_array = int_arrays[0]

    new_mz_array = np.array(new_mz_array)
    new_intensity_array = np.array(new_intensity_array)

    min_i = new_intensity_array.max() / dyn_range
    idx = new_intensity_array >= min_i
    new_intensity_array = new_intensity_array[idx]
    new_mz_array = new_mz_array[idx]

    return {'params': {'title': title, 'pepmass': pepmass,
            'rtinseconds': rtinseconds, 'charge': charge},
        'm/z array': new_mz_array,
        'intensity array': new_intensity_array}


def _lower_median_mass_index(masses):
    i = np.argsort(masses)
    k = (len(masses) - 1) // 2
    idx = i[k]
    return idx, masses[idx]

def lower_median_mass(spectra):
    masses, charges = _neutral_masses(spectra)
    i, m = _lower_median_mass_index(masses)
    z = charges[i]
    return (m + z*H) / z, z

def lower_median_mass_rt(spectra):
    masses, charges = _neutral_masses(spectra)
    rts = [s['params']['rtinseconds'] for s in spectra]
    i, m = _lower_median_mass_index(masses)
    return rts[i]

def get_cluster_id(title):
    return title.split(';', 1)[0]

def naive_average_mass_and_charge(spectra):
    mzs = [s['params']['pepmass'][0] for s in spectra]
    charges = {tuple(s['params']['charge']) for s in spectra}
    if len(charges) > 1:
        raise ValueError('There are different charge states in the cluster. Cannot average precursor m/z.')
    return sum(mzs) / len(mzs), charges.pop()[0]

def _neutral_masses(spectra):
    mzs = [s['params']['pepmass'][0] for s in spectra]
    charges = [s['params']['charge'][0] for s in spectra if len(s['params']['charge']) == 1]
    masses = [(m*c-c*H) for m, c in zip(mzs, charges)]
    return masses, charges

def neutral_average_mass_and_charge(spectra):
    masses, charges = _neutral_masses(spectra)
    z = int(round(sum(charges)/len(charges)))
    avg_mass = sum(masses) / len(masses)
    return (avg_mass + z*H) / z, z

def median_rt(spectra):
    rts = [s['params']['rtinseconds'] for s in spectra]
    return np.median(rts)


def process_maracluster_mgf(fname, get_cluster=get_cluster_id,
    get_pepmass=naive_average_mass_and_charge,
    get_rt=median_rt,
    **kwargs):
    outputs = []
    with mgf.IndexedMGF(fname, read_charges=False) as fin:
        titles = fin.index.keys()
        for cluster_id, cluster_group in groupby(titles, get_cluster):
            ids = list(cluster_group)
            spectra = fin[ids]
            mz, c = get_pepmass(spectra)
            rt = get_rt(spectra)
            outputs.append(average_spectrum(
                spectra, cluster_id, pepmass=mz, charge=c, rtinseconds=rt, **kwargs))
    return outputs


def main():
    pars = argparse.ArgumentParser()
    pars.add_argument('input', help='MGF file with clustered spectra.')
    pars.add_argument('output', nargs='?', help='Output file (default is stdout).')
    pars.add_argument('--mode', choices=['single', 'encoded_clusters'], default='encoded_clusters',
        help='Operation mode. Single: input MGF is interpreted as a single cluster.'
        'encoded_clusters: cluster IDs are parsed out of spectrum titles.')
    pars.add_argument('--dyn-range', type=float, default=DYN_RANGE,
        help='Dynamic range to apply to output spectra')
    pars.add_argument('--min-fraction', type=float, default=MIN_FRACTION,
        help='Minimum fraction of cluster spectra where MS/MS peak is present.')
    pars.add_argument('--mz-accuracy', type=float, default=DIFF_THRESH,
        help='Minimum distance between MS/MS peak clusters.')
    pars.add_argument('--append', action='store_true',
        help='Append to output file instead of replacing it.')
    pars.add_argument('--rt', choices=['median', 'mass_lower_median'], default='median')
    pars.add_argument('--pepmass', choices=['naive_average', 'neutral_average',
        'lower_median'], default='lower_median')
    pars.add_argument('--msms_avg', choices=['naive', 'weighted'], default='weighted')

    args = pars.parse_args()
    if args.pepmass == 'lower_median':
        args.rt = 'mass_lower_median'
    get_rt = {'median': median_rt, 'mass_lower_median': lower_median_mass_rt}[args.rt]
    get_pepmass = {'naive_average': naive_average_mass_and_charge,
        'neutral_average': neutral_average_mass_and_charge,
        'lower_median': lower_median_mass}[args.pepmass]

    kwargs = {'mz_accuracy': args.mz_accuracy, 'dyn_range': args.dyn_range,
        'min_fraction': args.min_fraction, 'msms_avg': args.msms_avg}
    mode = 'wa'[args.append]
    if args.mode == 'single':
        spectra = list(mgf.read(args.input))
        mz, c = get_pepmass(spectra)
        rt = get_rt(spectra)
        mgf.write([average_spectrum(spectra,
            title=args.output, pepmass=mz, charge=c, rtinseconds=rt, **kwargs)],
            args.output, file_mode=mode)
    elif args.mode == 'encoded_clusters':
        mgf.write(process_maracluster_mgf(args.input,
            get_pepmass=get_pepmass, get_rt=get_rt, **kwargs), args.output, file_mode=mode)
    else:
        raise NotImplementedError('This mode is not implemented yet.')


if __name__ == '__main__':
    main()