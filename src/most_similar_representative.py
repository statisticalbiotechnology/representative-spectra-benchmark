import pyopenms
import sys, getopt
import numpy as np

# definition of distance
# XQS.xCorrelationPrescore: simple, binned dot product. third parameter is binsize in Da
def distance(spec1, spec2, method='xcorr'):
    if (method=='xcorr'):
        XQS = pyopenms.XQuestScores()
        xcorr = XQS.xCorrelationPrescore(spec1, spec2, 0.1)
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
    exp = pyopenms.MSExperiment()
    MGF = pyopenms.MascotGenericFile()
    MGF.load(inputfile, exp)

    # compute distances for all pairs (full matrix, makes it easier to summarize but not very efficient)
    dist_matrix = np.zeros((exp.size(), exp.size()))
    for i in range(0,exp.size()):
        for j in range(0,exp.size()):
            dist_matrix[i][j] = distance(exp[i], exp[j])

    # total distance from one spectrum to all others
    total_dist = np.zeros(dist_matrix[0].size)
    for i in range(0, dist_matrix[0].size):
        total_dist[i] = dist_matrix[i].sum() / dist_matrix[0].size
    # find spectrum with minimal distance to all other spectra
    best_spec = np.where(total_dist == np.amin(total_dist))

    # write output file
    export_spec = pyopenms.MSExperiment()
    export_spec.addSpectrum(exp[best_spec[0].item()])
    MGF.store(outputfile, export_spec)


if __name__ == "__main__":
   main(sys.argv[1:])