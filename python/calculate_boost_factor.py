import numpy as np
import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as plt
import glob

## 0:bincenter 1:Deltasigma_num 2:Sumwls 3:Sumwlsms 4:SumwlsResp 5:DeltaSigma_selnbias_corrected
#6:Sumwls_by_sumwl 7:rmin/2+rmax/2 8:DeltaSigma_cross 9:Sum_srcwt_per_lens 10:TotalSum_srcwt
#11:SN_ErrDeltaSigma 12:SN_ErrDeltaSigma_cross 13:Total_pairs 14:m_selnbias 15:sumct_num 16: a_selnbias

#all_data=np.genfromtxt(sys.argv[1])
#pc_data=np.genfromtxt(sys.argv[2])
#cc_data_lowz=np.genfromtxt(sys.argv[3]) # zcl < 0.4
#cc_data_highz=np.genfromtxt(sys.argv[4]) #zcl >= 0.4


def getBoostFactorCov (randoms_list):
    """
    From Hironao:
    1. Make random realizations such that each realization has the same number of random points as our lens sample. If you have 1900 randoms, you could get 100 realizations.
    2. divide the random realizations into two groups such that each of the group has 50 realizations.
    3. From group B, you can compute the denominator in the boost factor, using all the 50 realizations.
    4. From group A, using a single realization, you can compute the numerator of the boost factor. You can then divide that by the denominator you computed in 3. Then you can get 50 realizations of boost factors.
    5. From that 50 realizations of boost factors, you can compute covariance, including correlation between radial bins.
    """
    # 1 randoms_list: 100 files from running weak lensing code
    # 2. Save 2:Sumwls from each file
    groupA, groupB = [], []
    for i, file_name in enumerate(randoms_list): 
        data = np.genfromtxt(file_name)
        if i < len(randoms_list)/2:
            groupA.append(data[1:,2])
        else:
            groupB.append(data[1:,2])
    groupA = np.array(groupA)
    groupB = np.array(groupB)
    #3
    demBoost = np.mean(groupA, axis=0)
    #4
    boost_vec = []
    for realization in groupB:
        boost = realization/demBoost
        boost_vec.append (boost)
    boost_vec = np.array(boost_vec)
    # shape is (50, rbins)
    # need to transpose to get a cov with shape (rbins, rbins)
    cov = np.cov(np.transpose(boost_vec), bias=True) # Use biased variance??
    corr = np.corrcoef (np.transpose(boost_vec))

    return cov, corr

signal=np.genfromtxt(sys.argv[1])
randoms=np.genfromtxt(sys.argv[2])

r = signal[1:, 7]
NClusters = 19
NRandoms = 1900 #1564 #len(randoms)
print ("NRandoms: " , NRandoms)

boost = (signal[1:, 2]/NClusters)/(randoms[1:, 2]/NRandoms)

# After runningthe weak lensign code 100 times, p-cuts, around random points (each file: 19 random points)
path_to_randoms = "/project/plazas/WORK/HSC/weaklens_pipeline/output_WL/100RANDOMS_2021JUL01/random_points_s19a*/Dsigma.dat"
randoms_list = glob.glob(path_to_randoms)

covBoost, corrBoost = getBoostFactorCov (randoms_list)

print ("boost: ", boost)
print ("error: ", np.sqrt(np.diag(covBoost)))

## CHI2

chi2 = np.dot(boost - 1., np.dot(np.linalg.inv(covBoost), boost-1))
print ("Chi2, all scales", chi2, chi2/len(boost))


chi2 = np.dot(boost[0:5] - 1., np.dot(np.linalg.inv(covBoost[0:5, 0:5]), boost[0:5]-1))
print ("Chi2, only small scales (first 5 points)", chi2, chi2/len(boost[0:5]))

print(boost[0:6])

from scipy.stats import chisquare
chi2_test  = chisquare(boost, f_exp=np.ones(len(boost)))
print ("Chi2 test, with p value: ", chi2_test)

chi2_test  = chisquare(boost[0:5], f_exp=np.ones(len(boost[0:5])))
print ("Chi2 test, with p value (first 5 points): ", chi2_test)


pp = PdfPages('boostFactor.pdf')

fig = plt.figure()
plt.errorbar (r, boost, yerr=np.sqrt(np.diag(covBoost)), fmt='o', color='green', label='boost factor', ls='none')
plt.hlines(1.0, np.min(r), np.max(r), colors='k', linestyles='--')
plt.suptitle('Boost factor')
plt.xscale('log')
plt.xlabel('R [Mpc/h]')
plt.ylabel('Boost')
plt.ylim([0.8, 1.2])

pp.savefig(fig)


fig = plt.figure()
plt.imshow(covBoost, origin='lower', cmap=plt.cm.BuPu_r)
plt.colorbar()
plt.suptitle("Covariance")
pp.savefig(fig)

fig = plt.figure()
plt.imshow(corrBoost, origin='lower', cmap=plt.cm.BuPu_r)
plt.colorbar()
plt.suptitle("Correlation")
pp.savefig(fig)


pp.close()


