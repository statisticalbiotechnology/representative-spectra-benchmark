import numpy as np


def bin_mean(peaklists, minimum=100, maximum=2000, binsize=0.02):

    array_size = int( (maximum - minimum ) / binsize ) + 1
    merged_spectrum = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
    merged_spectrum['intensities'] = np.zeros(array_size, dtype=np.float32)
    merged_spectrum['mzs'] = np.zeros(array_size, dtype=np.float32)
    merged_spectrum['n_peaks'] = np.zeros(array_size, dtype=np.int32)
    merged_spectrum['precursor_mzs'] = []
    merged_spectrum['precursor_charges'] = []

    for peaklist in peaklists:
        #### Convert the peak lists to np arrays
        intensity_array = np.asarray(peaklist['intensity array'])
        mz_array = np.asarray(peaklist['m/z array'])

        #### Limit the np arrays to the region we're interested in
        intensity_array = intensity_array[ ( mz_array >= merged_spectrum['minimum'] ) & ( mz_array < merged_spectrum['maximum'] ) ]
        mz_array = mz_array[ ( mz_array >= merged_spectrum['minimum'] ) & ( mz_array < merged_spectrum['maximum'] ) ]

        #### Compute their bin locations and store n_peaks and intensities
        bin_array = (( mz_array - merged_spectrum['minimum'] ) / merged_spectrum['binsize']).astype(int)

        merged_spectrum['n_peaks'][bin_array] += 1
        merged_spectrum['intensities'][bin_array] += intensity_array
        merged_spectrum['precursor_mzs'].append(peaklist['precursor mz'])
        merged_spectrum['precursor_charges'].append(peaklist['precursor charge'])

    # Check that all precursor charges are the same
    charges = merged_spectrum['precursor_charges']
    assert all(x == charges[0] for x in charges), "Not all precursor charges in cluster are equal"

    # Take the mean of all peaks per bin
    merged_spectrum['intensities'][merged_spectrum['intensities'] == 0] = np.nan
    merged_spectrum['intensities'] = np.divide(merged_spectrum['intensities'], merged_spectrum['n_peaks'])

    # Only return non-zero intensity bins
    nan_mask = ~np.isnan(merged_spectrum['intensities'])
    merged_spectrum['intensities'] = merged_spectrum['intensities'][nan_mask]
    merged_spectrum['mzs'] = np.arange(
        minimum + (binsize / 2), maximum + binsize, binsize, dtype=np.int32
    )[nan_mask]
    merged_spectrum['precursor_mz'] = np.mean(merged_spectrum['precursor_mzs'])
    merged_spectrum['precursor_charge'] = charges[0]

    del merged_spectrum['n_peaks']
    del merged_spectrum['precursor_charges']
    del merged_spectrum['precursor_mzs']

    return merged_spectrum


def write_spectrum(spectra, mgf_file):
    for i, spectrum in enumerate(spectra):
        mgf_tmp = f"""TITLE={i}
PEPMASS={spectrum['precursor_mz']}
CHARGE={spectrum['precursor_charge']}+
BEGIN IONS
"""
        for mz, intensity in zip(spectrum['mzs'], spectrum['intensities']):
            if not np.isnan(intensity):
                mgf_tmp += f"{mz} {intensity}\n"
        mgf_tmp += 'END IONS\n\n'
        mgf_file.write(mgf_tmp)


peaklists = [
    {
        'intensity array': [100, 200, 300],
        'm/z array': [100, 200, 300],
        'precursor mz': 555,
        'precursor charge': 3
    },
    {
        'intensity array': [100, 200, 300],
        'm/z array': [100, 200, 300],
        'precursor mz': 555,
        'precursor charge': 3
    }
]

merged_spectrum = bin_mean(peaklists)
merged_spectra = [merged_spectrum, merged_spectrum]

with open('merged.mgf', 'wt') as mgf_file:
    write_spectrum(merged_spectra, mgf_file)
