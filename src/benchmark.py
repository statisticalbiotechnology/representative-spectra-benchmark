import numpy as np
import sys
import spectrum_utils.spectrum as sus
from pyteomics import mgf, parser
from scipy.stats import binned_statistic
import metrics as mx

def test_method_metric(member_file,representative_file,metric=mx.cos_dist):
    with mgf.mgf(member_file) as memb_reader, mgf.mgf(representative_file) as rep_reader:
        m_spec = memb_reader.next()
        prev_clid = "none"
        while (r_spec = rep_reader.next()):
            m_specs = []
            curr_clid = prev_clid
            while prev_clid==curr_clid:
                m_specs += m_spec
                m_spec = memb_reader.next()
                curr_clid = m_spec["title"]
                print(curr_clid)
            prev_clid = curr_clid
            print("xxx")

if __name__ == "__main__":
    test_method_metric(
        "data/clusters_maracluster.mgf",
        "data/representativespectra/most_similar_spectra__binned_dot_product.mgf",
        mx.cos_dist
        )
