from pyopenms import *
import sys, getopt
import numpy as np
import pandas as pd

# input: MGF file with all clustered spectra
# output: MGF file with one representative spectrum per cluster
# use example: python most_similar_representative.py -i data/clusters_maracluster.mgf -o data/most_similar_cluster_representatives.mgf

# definition of distance
# XQuestScores.xCorrelationPrescore: simple, binned dot product, normalized by number of peaks.
# third parameter of xCorrelationPrescore is binsize in Da
def distance(spec1, spec2, method='xcorr'):
    if (method=='xcorr'):
        xcorr = XQuestScores().xCorrelationPrescore(spec1, spec2, 0.1)
        return 1.0-xcorr
    
    else:
        return 0


def main(argv):
   # handle arguments
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:")
    except getopt.GetoptError:
        print('most_similar_representative.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('most_similar_representative.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg
        elif opt in ("-o"):
            outputfile = arg


    # load clustered spectra
    exp = MSExperiment()
    MascotGenericFile().load(inputfile, exp)

    cluster_names = []
    cluster_membership = []

    # extract cluster names from TITLE
    for i in range(0,exp.size()):
        cl_name = exp[i].getMetaValue("TITLE").decode().split(";")[0]
        cluster_membership.append(cl_name)
        if cl_name not in cluster_names:
            cluster_names.append(cl_name)

    cluster_membership = pd.Series(cluster_membership)

    export_spec = MSExperiment()
    range_start = 0

    for cl in cluster_names:
        
        print(cl)
        # collect the spectra in the current cluster
        cluster_reached = False
        cluster_ended = False
        cluster_spec = []
        for i in range(range_start, cluster_membership.size):
            if (cl == cluster_membership[i]):
                cluster_spec.append(i)
                if (not cluster_reached): # found start of cluster
                    cluster_reached = True
            else: # reached end of cluster
                if (cluster_reached):
                    range_start = i-1
                    break
        
        # if the cluster only contains one spectrum, just return that one
        print(len(cluster_spec))
        if (len(cluster_spec) == 1):
            export_spec.addSpectrum(exp[cluster_spec[0]])
            continue
            
        # if the cluster does not have any spectra, skip it (should not happen)
        if (len(cluster_spec) == 0):
            continue
        
        # initialize distance matrix
        dist_matrix = np.zeros((len(cluster_spec), len(cluster_spec)))

        # calculate pairwise distances (fill triangular matrix)
        for i in range(0, len(cluster_spec)):
            for j in range(i, len(cluster_spec)):
                dist_matrix[i][j] = distance(exp[cluster_spec[i]], exp[cluster_spec[j]])
        
        dist_matrix = pd.DataFrame(dist_matrix)

        # summarize distances from each spetrum to all others
        total_dist = np.zeros(dist_matrix.iloc[0,:].size)
        for i in range(0, dist_matrix.iloc[0,:].size):
            total_dist[i] = (dist_matrix.iloc[i,:].sum() + dist_matrix.iloc[:,i].sum()) / dist_matrix.iloc[0,:].size

        # find best spectrum with minimal distance to all others
        best_spec = np.where(total_dist == np.amin(total_dist))

        # add best spectrum to set of exported spectra
        # first option is used, if several spectra have the same total distance
        best_spec = best_spec[0]
        if (best_spec.size > 1):
            best_spec = best_spec[0]
        best_spec = best_spec.item()
        export_spec.addSpectrum(exp[cluster_spec[best_spec]])

    # write output file
    print(export_spec.size())
    MascotGenericFile().store(outputfile, export_spec)


if __name__ == "__main__":
   main(sys.argv[1:])