import argparse
import numpy as np
from itertools import groupby
from pyteomics import mgf, mass

DIFF_THRESH = 0.01
DYN_RANGE = 1000
H = mass.nist_mass['H+'][0][0]


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

    Returns
    -------
    out : average spectrum in pyteomics format
    '''
    mz_accuracy = kwargs.get('mz_accuracy', DIFF_THRESH)
    dyn_range = kwargs.get('dyn_range', DYN_RANGE)
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

        new_mz_array.append(mz_array_sum[i_prev-1]/(i_prev))
        new_intensity_array.append(intensity_array_sum[i_prev-1]/n)

        for i in ind_list[1:-1]:
            new_mz_array.append((mz_array_sum[i-1]-mz_array_sum[i_prev-1])/(i-i_prev))
            new_intensity_array.append((intensity_array_sum[i-1]-intensity_array_sum[i_prev-1])/n)
            i_prev = i

        new_mz_array.append((mz_array_sum[-1] - mz_array_sum[i_prev-1])/(len(mz_array_sum) - i_prev))
        new_intensity_array.append((intensity_array_sum[-1] - intensity_array_sum[i_prev-1])/n)
    else:
        new_mz_array = mz_arrays[0]
        new_intensity_array = int_arrays[0]

    new_mz_array = np.array(new_mz_array)
    new_intensity_array = np.array(new_intensity_array)

    min_i = new_intensity_array.max()/dyn_range
    idx = new_intensity_array >= min_i
    new_intensity_array = new_intensity_array[idx]
    new_mz_array = new_mz_array[idx]

    return {'params': {'title': title, 'pepmass': pepmass,
            'rtinseconds': rtinseconds, 'charge': charge},
        'm/z array': new_mz_array,
        'intensity array': new_intensity_array}


def get_cluster_id(title):
    return title.split(';', 1)[0]

def naive_average_mass_and_charge(spectra):
    mzs = [s['params']['pepmass'][0] for s in spectra]
    charges = {tuple(s['params']['charge']) for s in spectra}
    if len(charges) > 1:
        raise ValueError('There are different charge states in the cluster. Cannot average precursor m/z.')
    return sum(mzs) / len(mzs), charges.pop()[0]

def neutral_average_mass_and_charge(spectra):
    mzs = [s['params']['pepmass'][0] for s in spectra]
    charges = [s['params']['charge'][0] for s in spectra if len(s['params']['charge']) == 1]
    masses = [(m*c-c*H) for m, c in zip(mzs, charges)]
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
    pars.add_argument('input', help='MGF file with clustered spectra')
    pars.add_argument('output', nargs='?', help='Output file (default is stdout)')
    mode = pars.add_mutually_exclusive_group(required=True)
    mode.add_argument('--single', action='store_true',
        help='If specified, input is interpreted as containing a single cluster.')
    mode.add_argument('--encodedclusters', action='store_true',
        help='Process an MGF with cluster IDs encoded in titles.')
    pars.add_argument('--dyn-range', type=float, default=DYN_RANGE,
        help='Dynamic range to apply to output spectra')
    pars.add_argument('--mz-accuracy', type=float, default=DIFF_THRESH,
        help='Minimum distance between MS/MS peak clusters.')
    pars.add_argument('--append', action='store_true',
        help='Append to output file instead of replacing it.')
    pars.add_argument('--rt', choices=['median'], default='median')
    pars.add_argument('--pepmass', choices=['naive_average', 'neutral_average'], default='naive_average')

    args = pars.parse_args()

    get_rt = {'median': median_rt}[args.rt]
    get_pepmass = {'naive_average': naive_average_mass_and_charge,
        'neutral_average': neutral_average_mass_and_charge}[args.pepmass]

    kwargs = {'mz_accuracy': args.mz_accuracy, 'dyn_range': args.dyn_range}
    mode = 'wa'[args.append]
    if args.single:
        spectra = list(mgf.read(args.input))
        mz, c = get_pepmass(spectra)
        rt = get_rt(spectra)
        mgf.write([average_spectrum(spectra,
            title=args.output, pepmass=mz, charge=c, rtinseconds=rt, **kwargs)],
            args.output, file_mode=mode)
    elif args.encodedclusters:
        mgf.write(process_maracluster_mgf(args.input,
            get_pepmass=get_pepmass, get_rt=get_rt, **kwargs), args.output, file_mode=mode)
    else:
        raise NotImplementedError('This mode is not implemented yet.')


if __name__ == '__main__':
    main()