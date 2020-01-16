import numpy as np
import sys
import spectrum_utils.spectrum as sus
from pyteomics import mgf, parser
from scipy.stats import binned_statistic
import metrics as mx
import ms_io

def test_method_metric(member_file,representative_file,metric=mx.cos_dist):
    # with mgf.mgf(member_file) as memb_reader, mgf.mgf(representative_file) as rep_reader:
    #     m_spec = memb_reader.next()
    #     m_specs = []
    #     # for r_spec in rep_reader:
    pass

if __name__ == "__main__":
    CLUSTER_FILE = "data/clusters_maracluster.mgf"
    # REPRASENTATIVE_FILE = "data/representativespectra/most_similar_spectra__binned_dot_product.mgf"
    REPRASENTATIVE_FILE = 'data/best_sp'
    # clusters = ms_io.read_cluster_spectra(CLUSTER_FILE)
    representatives = ms_io.read_cluster_spectra(REPRASENTATIVE_FILE, usi_present=False)
    # test_method_metric(
    #     mx.cos_dist
    #     )
