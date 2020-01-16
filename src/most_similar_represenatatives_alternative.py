import pyopenms
import numpy as np
import scipy.spatial as scisp



# definition of distance
# XQS.xCorrelationPrescore: simple, binned dot product. third parameter is binsize in Da
# other: similarity computed on aligned pairs of spectra (mz-features within 0.02 Da range considered as same features
#        all other features considered as different) and various distance measures available in scipy
def distance(spec1, spec2, method='xcorr'):
    if (method=='xcorr'):
        XQS = pyopenms.XQuestScores()
        xcorr = XQS.xCorrelationPrescore(spec1, spec2, 0.1)
        return 1.0-xcorr
    elif(method=='other'):
        mz1, i1 = spec1.get_peaks()
        mz2, i2 = spec2.get_peaks()

        tempMZ1 = np.zeros(len(mz1)+len(mz2),float)
        tempMZ2 = np.zeros(len(mz1)+len(mz2),float)
        tempI1 = np.zeros(len(i1)+len(i2),float)
        tempI2 = np.zeros(len(i1)+len(i2),float)

        lastJ = 0
        idx = 0

        for i in range(0, len(mz1)):
            for j in range(lastJ, len(mz2)):
                if(abs(mz1[i]-mz2[j]) <= 0.02 ):
                    tempMZ1[idx] = mz1[i]
                    tempMZ2[idx] = mz1[i]
                    tempI1[idx] = i1[i]
                    tempI2[idx] = i2[j]
                    idx = idx+1
                    lastJ = lastJ+1
                    break
                elif(mz1[i]-mz2[j] < -0.02):
                    tempMZ1[idx] = mz1[i]
                    tempMZ2[idx] = mz1[i]
                    tempI1[idx] = i1[i]
                    tempI2[idx] = 0
                    idx = idx + 1
                    break
                elif(i==len(mz1)-1 & (mz1[i]-mz2[j] > 0.02)):
                    tempMZ1[idx] = mz2[j]
                    tempMZ2[idx] = mz2[j]
                    tempI1[idx] = 0
                    tempI2[idx] = i2[j]
                    idx = idx + 1
                    lastJ = lastJ + 1
                    continue
                elif((i<len(mz1)-1) & (mz1[i]-mz2[j] > 0.02)):
                        if((mz2[j] < mz1[i+1]-0.02) & (mz1[i]-mz2[j] > 0.02)):
                            tempMZ1[idx] = mz2[j]
                            tempMZ2[idx] = mz2[j]
                            tempI1[idx] = 0
                            tempI2[idx] = i2[j]
                            idx = idx + 1
                            lastJ = lastJ+1
                            continue
                else:
                    break
        #euclidean distance
        #return scisp.distance.euclidean(tempI1, tempI2)
        #correlation-based distance
        #return scisp.distance.correlation(tempI1, tempI2)
        #manhattan/cityblock distance
        return scisp.distance.cityblock(tempI1, tempI2)
    else:
        return 0



exp = pyopenms.MSExperiment()
mgf = pyopenms.MascotGenericFile()
mgf.load("C:\\UNI\Projects_BIG\\de.NBI\\Reisen\\20200113-17_EuBIC_Developers_Meeting_Nyborg\\Hackathon\\clusters_mracluster.mgf", exp)

export_spec = pyopenms.MSExperiment()
outputfile = "C:\\UNI\Projects_BIG\\de.NBI\\Reisen\\20200113-17_EuBIC_Developers_Meeting_Nyborg\\Hackathon\\clusters_mracluster_representatives.mgf"

lastJ = 0
for i in range(0, exp.size()):
    if(i <= lastJ):
        continue
        print("continue")

    subExp = pyopenms.MSExperiment()
    subExp.addSpectrum(exp[i])

    clusterID = exp[i].getMetaValue("TITLE").decode().split(";")[0]

    for j in range(i+1, exp.size()):

        lastJ = j-1

        if (i == exp.size()):
            continue
            print("last-element - continue")

        if(exp[j].getMetaValue("TITLE").decode().split(";")[0] == clusterID):
            subExp.addSpectrum(exp[j])
            #print("match!")
        else:
            break

    print("clusterID:", clusterID, ", size:", subExp.size())



    # compute distances for all pairs (full matrix, makes it easier to summarize but not very efficient)
    dist_matrix = np.zeros((subExp.size(), subExp.size()))
    for i in range(0, subExp.size()):
        for j in range(i+1, subExp.size()):
            dist_matrix[i][j] = distance(subExp[i], subExp[j], 'other')



    # total distance from one spectrum to all others
    total_dist = np.zeros(dist_matrix[0].size)
    for i in range(0, dist_matrix[0].size):
        total_dist[i] = dist_matrix[i].sum() / dist_matrix[0].size
    # find spectrum with minimal distance to all other spectra
    best_spec = np.where(total_dist == np.amin(total_dist))



    # write to output file
    if(len(best_spec[0]) == 1):
        export_spec.addSpectrum(subExp[best_spec[0].item()])
    else:
        export_spec.addSpectrum(subExp[best_spec[0][0].item()])



#write complete output file
mgf.store(outputfile, export_spec)
