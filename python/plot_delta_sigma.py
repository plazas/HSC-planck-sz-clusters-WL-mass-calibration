import numpy as np
import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as plt

## 0:bincenter 1:Deltasigma_num 2:Sumwls 3:Sumwlsms 4:SumwlsResp 5:DeltaSigma_selnbias_corrected
#6:Sumwls_by_sumwl 7:rmin/2+rmax/2 8:DeltaSigma_cross 9:Sum_srcwt_per_lens 10:TotalSum_srcwt
#11:SN_ErrDeltaSigma 12:SN_ErrDeltaSigma_cross 13:Total_pairs 14:m_selnbias 15:sumct_num 16: a_selnbias

#all_data=np.genfromtxt(sys.argv[1])
#pc_data=np.genfromtxt(sys.argv[2])
#cc_data_lowz=np.genfromtxt(sys.argv[3]) # zcl < 0.4
#cc_data_highz=np.genfromtxt(sys.argv[4]) #zcl >= 0.4


planck_centers_p_cuts=np.genfromtxt(sys.argv[1])
camira_centers_p_cuts=np.genfromtxt(sys.argv[2])
planck_centers_cc_cuts_lowz=np.genfromtxt(sys.argv[3])
planck_centers_cc_cuts_highz=np.genfromtxt(sys.argv[4])
camira_centers_cc_cuts_lowz=np.genfromtxt(sys.argv[5])
camira_centers_cc_cuts_highz=np.genfromtxt(sys.argv[6])


def delta_sigma_combine_cc_cuts_files (lowzFile, highzFile):
    """
    Output of Surhud's WL code. The CC cuts by Medezinski
    have different parameters for z_cl< 0.4 and z_cl >= 0.4. 
   
    To cover all possibilities, the WL code shoudl be run twice
    for each type of clusters, and then the desired output
    combined with this function.

    We are interested in delta_sigma and delta_sigma_x (the B mode).

    delta sigma, column #5
     (sumdsig_num[i]/sumwls_ms[i]/2./Resp-sumct_num[i]/sumwls_ms[i] -
     aselnbias*sumdsig_psf_num[i]/sumwls[i])/(1.0+mselnbias) 
    - sumdsig_num[i], is column # 1
    - sumwls_ms[i] is column #3
    - Resp = sumwls_resp[i]/sumwls[i]
      - sumwls_resp[i] is column #4
      - sumwls[i] is columns #2
    - sumct_num[i] is column #15 
    - aselnbiasFactor:
      - The whole combination (aselnbias*sumdsig_psf_num[i]/sumwls[i]) is column #16 
    - mselnbias is column #14

    delta sigma cross, column 8 
       sumdcross_num[i]/sumwls_ms[i]/2./Resp-sumcx_num[i]/sumwls_ms[i]-
       aselnbias*sumdcross_psf_num[i]/sumwls[i]

     Needed to add in the output columns 17, 18, 19
     - sumdcross_num[i] is column #18   (added in this branch) 
     - sumwls_ms[i] is column #3
     - Resp = sumwls_resp[i]/sumwls[i]
       - sumwls_resp[i] is column #4
       - sumwls[i] is columns #2
     - sumcx_num[i] is column #19  (added in this branch)
     - sumwls_ms[i] is column #3
     - aselnbiasFactorCross: 
        - The whole combination (aselnbias*sumdcross_psf_num[i]/sumwls[i]) is column #17 (added inthis branch)
     

    """
    sumdsig_num = lowzFile[1:,1] + highzFile[1:,1]
    sumwls_ms = lowzFile[1:,3] + highzFile[1:,3]
    Resp = (lowzFile[1:,4] + highzFile[1:,4])/ (lowzFile[1:,2] + highzFile[1:,2])
    sumct_num = lowzFile[1:,15] + highzFile[1:,15]
    aselnbiasFactor = lowzFile[1:,16] + highzFile[1:,16] # aselnbias*sumdsig_psf_num[i]/sumwls[i]
    mselnbias = lowzFile[1:,14] + highzFile[1:,14]

    newDeltaSigma = (sumdsig_num/sumwls_ms/2./Resp-sumct_num/sumwls_ms - aselnbiasFactor)/(1.0+mselnbias)

    #cross component
    sumdcross_num = lowzFile[1:,18] + highzFile[1:,18]
    sumcx_num = lowzFile[1:,19] + highzFile[1:,19]
    aselnbiasFactorCross =  lowzFile[1:,17] + highzFile[1:,17]

    newDeltaSigmaCross = sumdcross_num/sumwls_ms/2./Resp-sumcx_num/sumwls_ms - aselnbiasFactorCross


    # Add errors in quadrature
    SN_ErrDeltaSigma = np.sqrt(lowzFile[1:,11]**2 + lowzFile[1:,11]**2)
    SN_ErrDeltaSigma_cross = np.sqrt(lowzFile[1:,12]**2 + lowzFile[1:,12]**2)

    return newDeltaSigma, newDeltaSigmaCross, SN_ErrDeltaSigma, SN_ErrDeltaSigma_cross


# Copy one of the cc files and then overwrite columns
camira_centers_cc_cuts = np.copy(camira_centers_cc_cuts_lowz)

newDeltaSigma, newDeltaSigmaCross, SN_ErrDeltaSigma, SN_ErrDeltaSigma_cross = delta_sigma_combine_cc_cuts_files (camira_centers_cc_cuts_lowz, camira_centers_cc_cuts_highz)
camira_centers_cc_cuts[1:,5] = newDeltaSigma
camira_centers_cc_cuts[1:,8] = newDeltaSigmaCross
camira_centers_cc_cuts[1:,11] = SN_ErrDeltaSigma
camira_centers_cc_cuts[1:,12] = SN_ErrDeltaSigma_cross

planck_centers_cc_cuts = np.copy(planck_centers_cc_cuts_lowz)

newDeltaSigma, newDeltaSigmaCross, SN_ErrDeltaSigma, SN_ErrDeltaSigma_cross = delta_sigma_combine_cc_cuts_files (planck_centers_cc_cuts_lowz, planck_centers_cc_cuts_highz)
planck_centers_cc_cuts[1:,5] = newDeltaSigma
planck_centers_cc_cuts[1:,8] = newDeltaSigmaCross
planck_centers_cc_cuts[1:,11] = SN_ErrDeltaSigma
planck_centers_cc_cuts[1:,12] = SN_ErrDeltaSigma_cross


#    a[0].errorbar (pc_data[1:,7]+ 0.07*pc_data[1:,7] , pc_data[1:,7]*pc_data[1:,5], yerr=
#    np.fabs(pc_data[1:,7])*pc_data[1:,11], fmt='o', color='orange', label="P(z) cut", ls='none')#,
#    marker='.')

#   a[1].errorbar (pc_data[1:,7]+ 0.07*pc_data[1:,7], pc_data[1:,7]*pc_data[1:,8],
#   yerr=np.fabs(pc_data[1:,7])*pc_data[1:,12], fmt='o', color='orange', label="P(z) cut", ls='none')#,
#   marker='.')


pp=PdfPages("delta_sigma_planck_s19a.pdf")

def plotDeltaSigma (file1, file2, pdfPages, legend1='', legend2='', multByR=False, title=''):
    """
    y1 : DeltaSigma
    y1Cross : DeltaSigmaCross
    """
    nRows=2
    nCols=1

    legendFontSize = 6.5
    labelFontSize = 16
    titleFontSize = 9
    supTitleFontSize = 18
    markerSize = 25
    
    r1, y1, y1err = file1[1:,7], file1[1:,5], file1[1:,11] 
    y1Cross, y1CrossErr = file1[1:,8], file1[1:,12]
    
    r2, y2, y2err = file2[1:,7], file2[1:,5], file2[1:,11]
    y2Cross, y2CrossErr = file2[1:,8], file2[1:,12]


    xLabel=r'R ($h^{-1}$ Mpc)'
    yLabel=r'$\Delta \Sigma_{+}$ ($M_{\odot}/pc^2$)'
    yLabelCross=r'$\Delta \Sigma_{x}$ ($M_{\odot}/pc^2$)'

    if multByR:
        y1*=r1
        y1err*=np.fabs(r1)

        y1Cross*=r1
        y1CrossErr*=np.fabs(r1)

        y2*=r2
        y2err*=np.fabs(r2)

        y2Cross*=r2
        y2CrossErr*=np.fabs(r2)

        xLabel=r'R ($h^{-1}$ Mpc)'
        yLabel=r'$\Delta \Sigma_{+} \times R$ ($Mpc\ M_{\odot}/pc^2$)'
        yLabelCross=r'$\Delta \Sigma_{x} \times R$ ($Mpc\ M_{\odot}/pc^2$)'


    fig, ax = plt.subplots(nrows=nRows, ncols=nCols, sharex='col', sharey='row', figsize=(13, 10),
                           gridspec_kw={'height_ratios': [3, 1]})

    legend1=""
    legend2=""
    # Delta Sigma
    ax[0].set_xlabel(xLabel, fontsize=labelFontSize)
    ax[0].set_ylabel(yLabel, fontsize=labelFontSize)
    ax[0].tick_params(labelsize=labelFontSize)
    ax[0].set_xscale('log')
    ax[0].set_yscale('linear')
    ax[0].errorbar (r1+ 0.14*r1, y1, yerr=y1err, fmt='o', color='green', label=legend1, ls='none')#, marker='.')
    ax[0].errorbar (r2+ 0.07*r2, y2, yerr=y2err, fmt='o', color='orange', label=legend2, ls='none')#, marker='.')
    #ax[0].set_xlim(0.2, 5)
    #ax[0].set_ylim(0, 400)
    ax[0].legend(loc='upper right', prop={'size': 12})

    # Cross component
    ax[1].set_xlabel(xLabel, fontsize=labelFontSize)
    ax[1].set_ylabel(yLabelCross, fontsize=labelFontSize)
    ax[1].tick_params(labelsize=labelFontSize)
    ax[1].set_xscale('log')
    ax[1].set_yscale('linear')
    ax[1].errorbar (r1+ 0.14*r1, y1Cross, yerr=y1CrossErr, fmt='o', color='green', label="", ls='none')#, marker='.')
    ax[1].errorbar (r2+ 0.07*r2, y2Cross, yerr=y2CrossErr, fmt='o', color='orange', label="", ls='none')#, marker='.')
    #ax[1].set_xlim(0.2, 5)
    #ax[1].set_ylim(-100, 200)

    fig.tight_layout()
    #fig.suptitle(title)
    pp.savefig(fig)
    plt.clf()
    multByR = False


# Planck vs Camira centers for P-cuts
filePairs = {1: [planck_centers_cc_cuts, camira_centers_cc_cuts, 'Planck centers - CC cut', 'Camira centers - CC cut'], 
             2: [camira_centers_p_cuts,  planck_centers_p_cuts, 'Camira centers - P cut', 'Planck centers - P cut'], 
             3: [planck_centers_cc_cuts, planck_centers_p_cuts, 'Planck centers - CC cut', 'Planck centers - P cut'], 
             4: [camira_centers_p_cuts, camira_centers_cc_cuts, 'Camira centers - P cut', 'Camira centers - CC cut'],
            }

filePairs = {4: [camira_centers_p_cuts, camira_centers_cc_cuts, 'Camira centers - P cut', 'Camira centers - CC cut']}
filePairs = {3: [planck_centers_cc_cuts, planck_centers_p_cuts, 'Planck centers - CC cut', 'Planck centers - P cut']}
filePairs = {2: [camira_centers_p_cuts,  planck_centers_p_cuts, 'Camira centers - P cut', 'Planck centers - P cut']}
filePairs = {1: [planck_centers_cc_cuts, camira_centers_cc_cuts, 'Planck centers - CC cut', 'Camira centers - CC cut']}

filePairs = {2: [camira_centers_p_cuts,  planck_centers_p_cuts, 'Camira centers - P cut, Planck redshift',
'Camira centers - CC cut, Planck redshift']}


filePairs = {2: [camira_centers_p_cuts,  planck_centers_p_cuts, 'Medezinski+18 centers, redshift - P cut',
    'Medezinski+18 centers, redshift - CC cut']}

# Planck vs Camira centers for CC-cuts
# CC-cuts vs P-cuts for Planck centers
# CC-cuts vs P-cuts for Camira centers

for key in filePairs:
    for multByR in [False, True]:
        print (key, multByR, filePairs[key][2], filePairs[key][3])
        plotDeltaSigma (filePairs[key][0], filePairs[key][1], pp,
                        legend1=filePairs[key][2], legend2=filePairs[key][3], multByR=multByR,
                        title=f"{filePairs[key][2]} and {filePairs[key][3]}. Medezinski+18")
pp.close()

# S/N calculation
for x in [planck_centers_cc_cuts, camira_centers_cc_cuts, planck_centers_p_cuts, camira_centers_p_cuts]:
    dsigma  = x[1:,5] 
    dsigma_err = x[1:, 11]
    SNR = np.sqrt(np.sum(dsigma**2/dsigma_err**2))
    print ("SNR: ", SNR)

