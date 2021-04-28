import numpy as np
import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as plt

## 0:bincenter 1:Deltasigma_num 2:Sumwls 3:Sumwlsms 4:SumwlsResp 5:DeltaSigma_selnbias_corrected
#6:Sumwls_by_sumwl 7:rmin/2+rmax/2 8:DeltaSigma_cross 9:Sum_srcwt_per_lens 10:TotalSum_srcwt
#11:SN_ErrDeltaSigma 12:SN_ErrDeltaSigma_cross 13:Total_pairs 14:m_selnbias 15:sumct_num 16: a_selnbias

all_data=np.genfromtxt(sys.argv[1])
pc_data=np.genfromtxt(sys.argv[2])
cc_data_lowz=np.genfromtxt(sys.argv[3]) # zcl < 0.4
cc_data_highz=np.genfromtxt(sys.argv[4]) #zcl >= 0.4


def delta_sigma_combine_cc_cuts_files (lowzFile, highzFile):
    """ delta sigma, column #5
     (sumdsig_num[i]/sumwls_ms[i]/2./Resp-sumct_num[i]/sumwls_ms[i] -
     aselnbias*sumdsig_psf_num[i]/sumwls[i])/(1.0+mselnbias) 
    - sumdsig_num[i], is column # 1
    - sumwls_ms[i] is columns #3
    - Resp = sumwls_resp[i]/sumwls[i]
    - sumwls_resp[i] is column #4
    - sumwls[i] is columns #2
    - sumct_num[i] is column #15 
    - The whole combination (aselnbias*sumdsig_psf_num[i]/sumwls[i]) is column #16 
    - mselnbias is column #14
    """
    sumdsig_num = lowzFile[1:,1] + highzFile[1:,1]
    sumwls_ms = lowzFile[1:,3] + highzFile[1:,3]
    Resp = (lowzFile[1:,4] + highzFile[1:,4])/ (lowzFile[1:,2] + highzFile[1:,2])
    sumct_num = lowzFile[1:,15] + highzFile[1:,15]
    aselnbiasFactor = lowzFile[1:,16] + highzFile[1:,16] # aselnbias*sumdsig_psf_num[i]/sumwls[i]
    mselnbias = lowzFile[1:,14] + highzFile[1:,14]

    newDeltaSigma = (sumdsig_num/sumwls_ms/2./Resp-sumct_num/sumwls_ms - aselnbiasFactor)/(1.0+mselnbias)

    # Add errors in quadrature
    SN_ErrDeltaSigma = np.sqrt(lowzFile[1:,11]**2 + lowzFile[1:,11]**2)
    SN_ErrDeltaSigma_cross = np.sqrt(lowzFile[1:,12]**2 + lowzFile[1:,12]**2)

    return newDeltaSigma, SN_ErrDeltaSigma, SN_ErrDeltaSigma_cross


# Copy one of the cc files and then overwrite columns

cc_data = np.copy(cc_data_lowz)
newDeltaSigma, SN_ErrDeltaSigma, SN_ErrDeltaSigma_cross = delta_sigma_combine_cc_cuts_files (cc_data_lowz,
                                                          cc_data_highz)

cc_data[1:,5] = newDeltaSigma
cc_data[1:,11] = SN_ErrDeltaSigma
cc_data[1:,12] = SN_ErrDeltaSigma_cross


nRows=2
nCols=1

legendFontSize = 6.5
labelFontSize = 16
titleFontSize = 9
supTitleFontSize = 18
markerSize = 25


pp=PdfPages("delta_sigma_planck_s19a.pdf")

f, a = plt.subplots(nrows=nRows, ncols=nCols, sharex='col', sharey='row', figsize=(13, 10),
                    gridspec_kw={'height_ratios': [3, 1]})

a[0].set_xlabel(r'R ($h^{-1}$ Mpc)', fontsize=labelFontSize)
a[0].set_ylabel(r'$\Delta \Sigma_{+} \times R$ ($Mpc\ M_{\odot}/pc^2$)', fontsize=labelFontSize)
a[0].tick_params(labelsize=labelFontSize)
a[0].set_xscale('log')
a[0].set_yscale('linear')

#a[0].errorbar (all_data[1:,7] + 0.14*all_data[1:,7], all_data[1:,7]*all_data[1:,5], all_data[1:,11], fmt='o', color='green', label="No cut", ls='none')#, marker='.')
#a[0].errorbar (pc_data[1:,7] + 0.07*pc_data[1:,7], pc_data[1:,7]*pc_data[1:,5], pc_data[1:,11], fmt='o', color='orange', label="P(z) cut", ls='none')#, marker='.')
#a[0].errorbar (cc_data[1:,7], cc_data[1:,7]*cc_data[1:,5], cc_data[1:,11], fmt='o', color='blue', label="CC cut", ls='none')#, marker='.')

a[0].errorbar (all_data[1:,7]+ 0.14*all_data[1:,7] , all_data[1:,7]*all_data[1:,5], yerr = np.fabs(all_data[1:,7])*all_data[1:,11], fmt='o', color='green', label="No cut", ls='none')#, marker='.')
a[0].errorbar (pc_data[1:,7]+ 0.07*pc_data[1:,7] , pc_data[1:,7]*pc_data[1:,5], yerr= np.fabs(pc_data[1:,7])*pc_data[1:,11], fmt='o', color='orange', label="P(z) cut", ls='none')#, marker='.')
a[0].errorbar (cc_data[1:,7], cc_data[1:,7]*cc_data[1:,5], yerr=np.fabs(cc_data[1:,7])*cc_data[1:,11], fmt='o', color='blue', label="CC cut", ls='none')#, marker='.')


a[0].set_xlim(0.1, 100)
#a[0].set_ylim(1., 1000)
a[0].legend(loc='upper right', prop={'size': 12})

a[1].set_xlabel(r'R ($h^{-1}$ Mpc)', fontsize=labelFontSize)
a[1].set_ylabel(r'$\Delta \Sigma_{x} \times R$ ($Mpc\ M_{\odot}/pc^2$)', fontsize=labelFontSize)
a[1].tick_params(labelsize=labelFontSize)
a[1].set_xscale('log')
a[1].set_yscale('linear')
a[1].errorbar (all_data[1:,7]+ 0.14*all_data[1:,7], all_data[1:,7]*all_data[1:,8],
        yerr=np.fabs(all_data[1:,7])*all_data[1:,12], fmt='o', color='green', label="No cut", ls='none')#, marker='.')
a[1].errorbar (pc_data[1:,7]+ 0.07*pc_data[1:,7], pc_data[1:,7]*pc_data[1:,8],
        yerr=np.fabs(pc_data[1:,7])*pc_data[1:,12], fmt='o', color='orange', label="P(z) cut", ls='none')#, marker='.')
a[1].errorbar (cc_data[1:,7], cc_data[1:,7]*cc_data[1:,8], yerr=np.fabs(cc_data[1:,7])*cc_data[1:,12], fmt='o', color='blue', label="CC cut", ls='none')#, marker='.')
a[1].set_xlim(0.1, 100)
a[1].set_ylim(-100, 100)

f.tight_layout()
pp.savefig()


f, a = plt.subplots(nrows=nRows, ncols=nCols, sharex='col', sharey='row', figsize=(13, 10),
                    gridspec_kw={'height_ratios': [3, 1]})

a[0].set_xlabel(r'R ($h^{-1}$ Mpc)', fontsize=labelFontSize)
a[0].set_ylabel(r'$\Delta \Sigma_{+} ($ M_{\odot}/pc^2$)', fontsize=labelFontSize)
a[0].tick_params(labelsize=labelFontSize)
a[0].set_xscale('log')
a[0].set_yscale('log')

#a[0].errorbar (all_data[1:,7] + 0.14*all_data[1:,7], all_data[1:,7]*all_data[1:,5], all_data[1:,11], fmt='o', color='green', label="No cut", ls='none')#, marker='.')
#a[0].errorbar (pc_data[1:,7] + 0.07*pc_data[1:,7], pc_data[1:,7]*pc_data[1:,5], pc_data[1:,11], fmt='o', color='orange', label="P(z) cut", ls='none')#, marker='.')
#a[0].errorbar (cc_data[1:,7], cc_data[1:,7]*cc_data[1:,5], cc_data[1:,11], fmt='o', color='blue', label="CC cut", ls='none')#, marker='.')

a[0].errorbar (all_data[1:,7]+ 0.14*all_data[1:,7], all_data[1:,5], yerr = all_data[1:,11], fmt='o', color='green', label="No cut", ls='none')#, marker='.')
a[0].errorbar (pc_data[1:,7]+ 0.07*pc_data[1:,7], pc_data[1:,5], yerr= pc_data[1:,11], fmt='o', color='orange', label="P(z) cut", ls='none')#, marker='.')
a[0].errorbar (cc_data[1:,7], cc_data[1:,5], yerr=cc_data[1:,11], fmt='o', color='blue', label="CC cut", ls='none')#, marker='.')


a[0].set_xlim(0.1, 100)
#a[0].set_ylim(1., 1000)
a[0].legend(loc='upper right', prop={'size': 12})

a[1].set_xlabel(r'R ($h^{-1}$ Mpc)', fontsize=labelFontSize)
a[1].set_ylabel(r'$\Delta \Sigma_{x} $ ($M_{\odot}/pc^2$)', fontsize=labelFontSize)
a[1].tick_params(labelsize=labelFontSize)
a[1].set_xscale('log')
a[1].set_yscale('log')
a[1].errorbar (all_data[1:,7]+ 0.14*all_data[1:,7], all_data[1:,8],
        yerr=all_data[1:,12], fmt='o', color='green', label="No cut", ls='none')#, marker='.')
a[1].errorbar (pc_data[1:,7]+ 0.07*pc_data[1:,7], pc_data[1:,8],
        yerr=pc_data[1:,12], fmt='o', color='orange', label="P(z) cut", ls='none')#, marker='.')
a[1].errorbar (cc_data[1:,7], cc_data[1:,8], yerr=cc_data[1:,12], fmt='o', color='blue', label="CC cut", ls='none')#, marker='.')
a[1].set_xlim(0.1, 100)
a[1].set_ylim(-100, 100)

f.tight_layout()
pp.savefig()











pp.close()

# S/N calculation

for x in [all_data, pc_data, cc_data]:
    dsigma  = x[1:,7] 
    dsigma_err = x[1:, 11]
    SNR = np.sqrt(np.sum(dsigma**2/dsigma_err**2))
    print ("SNR: ", SNR)



