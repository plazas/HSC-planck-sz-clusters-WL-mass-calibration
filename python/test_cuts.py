import numpy as np
import astropy.io.fits as pf
import ColorColor as cc
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages  # AP

# Color-Color cuts from Medezinski+18
def SelectSources(table, r, g1, g2, b, CCparams, rlim, blim, ext='mag_forced_cmodel'):
    """
    SelectSources(table, r,g1,g2,b, CCparams, rlim, blim)
    select red and blue background galaxies in color-color space from a sample of galaxies, 
    given parameters that describe the red sequence in CCparams, 
    and a set of cuts defined in rlim and blim, respectively.
  
    y-axis is b-g2 color, x-axis is g1-r colors, where r,g1,g2,b define the 
    filter column name first letter (e.g., 'z','i','r','g') and the rest of the column name goes in ext, 
    e.g. ext='mag_forced_cmodel'.
    e.g.,
    aCM, bCM, aCC, bCC = [-0.025, 1.6, 2.276, -0.152]
    rlim = [0.5, 21, 28, -0.5, 2.]
    blim = [-3, -0.8, 22, 28, 0.5, 0.5]

    r, g1, g2, b---> z, r, i, g

    Parameters
    ----------
    table : `astropy.io.fits.fitsrec.FITS_rec`
        Input catalog
    r : `str`
        Color prefix: y-axis is b-g2 color, x-axis is g1-r colors.
  
    g1 : `str`
        Color prefix: y-axis is b-g2 color, x-axis is g1-r colors.
    g2 : `str`
        Color prefix: y-axis is b-g2 color, x-axis is g1-r colors.
    b : `str`
        Color prefix: y-axis is b-g2 color, x-axis is g1-r colors.
    CCParams : `list`
        List with four components: aCM, bCM, aCC, bCC 
            aCM : Slope of line that follows red-sequence in g1-r as a function
                  of "r" magnitude (Eq. A3 in Medezinski+18)
            bCm : Ordinate intercept of line line that follows red-sequence in 
                  g1-r as a function of "r" magnitude (Eq. A3 in Medezinski+18).
            aCC : Slope of line that follows the red-sequence in "b-g2" and 
                  "g1-r" color space (Eq. A1 in Medezinski+18).
            bCC : Ordinate intercept of the line that follows the red-sequence in "b-g2" and 
                  "g1-r" color space (Eq. A1 in Medezinski+18).
    rlim : `list`
        List with cuts in CC space for "red" galaxies, see Appendix A1 in
        Medezinski+18
    
    blim : `list`
        List with cuts in CC space for "blue" galaxies, see Appendix A1 in
        Medezinski+18
    Return
    ------
       red : `astropy.io.fits.fitsrec.FITS_rec`
          Subset of input catalog with "red" background galaxies.
       blue : `astropy.io.fits.fitsrec.FITS_rec`
          Subset of input catalog wiht "blue" background galaxies.
       blueOrRedflag : `list`[`bool`]
           List with flags that identify galaxies that are "red" or "blue".
          
    Notes
    -----
    The original code is from https://github.com/elinor-lev/SourceSelection,
    implementing cuts in the appendix of Medezinski+18 (1706.00427).
    For filters like in HSC, "g1-r" is "r-z", and "b-g2" is "g-i"
    (RZ and BR in the code below, respectively).
    Medezinski+18 define two sets of "rlim" and "blim" cuts for clusters
    at zl >=0.4 and zl < 0.4.
    """
    aCM, bCM, aCC, bCC = CCparams
    
    # Correct for extinction. It is the columns 'a_[grizy]' in the catalog.
    table[r+ext]-=table['a_'+r.split('_')[1]]
    table[g1+ext]-=table['a_'+g1.split('_')[1]]
    table[g2+ext]-=table['a_'+g2.split('_')[1]]
    table[b+ext]-=table['a_'+b.split('_')[1]]
    
    Z = table[r+ext] #- table['a_'+]
    RZ = table[g1+ext] - table[r+ext]    # r-z
    BR = table[b+ext] - table[g2+ext]    # g-i
    
    f1 = aCM*Z + bCM # color-mag R-I seq
    seqdif1 = RZ - f1 # Eq. 34
    CCf = aCC * RZ + bCC          #color-color cluster sequence
    CCf2 = -1./aCC*RZ - bCC/aCC**2 # line perpendicular to color-color sequence
    CCdif = BR-CCf                 # B-R -CC Sequence
    CCdif2 = (BR-CCf2)/(1.+1./aCC**2)# B-R -perCC Sequence, blue\Deltacolor#2
    
    r_rzlim    = ( RZ > rlim[0] ) # RZ lower limit, separate red from blue / Eq. 28
    r_maglim   =  ( Z > rlim[1] ) & ( Z < rlim[2] ) # magnitude limit, redderst band  / Eq. 29
    r_CCseqlim = ( CCdif < rlim[3] ) &  ( CCdif2 < rlim[4] ) # Eq. 26 and Eq. 27
    red = table[ r_maglim & r_rzlim & r_CCseqlim ]
    # (seqdif1 < rlim[1] # doesn't exist...
    
    b_CMseqlim = (seqdif1 < blim[1]) & (seqdif1 > blim[0]) & (BR<4)
    b_maglim   =  ( Z > blim[2] ) & ( Z < blim[3] ) # magnitude limit, redderst band
    b_rzlim    = ( RZ < blim[4]) # RZ upper limit, separate red from blue
    b_CCseqlim = ( CCdif2 < blim[5] )
    blue =  table[  b_rzlim & (b_CCseqlim |  b_CMseqlim )  & b_maglim ]
    
    redFlag = r_maglim & r_rzlim & r_CCseqlim
    blueFlag = b_rzlim & (b_CCseqlim |  b_CMseqlim )  & b_maglim
    
    blueOrRedFlag = redFlag | blueFlag
    
    return red, blue, blueOrRedFlag


"""
ColDefs(
    name = 'object_id'; format = 'K'
    name = 'parent_id'; format = 'K'
    name = 'ira'; format = 'D'
    name = 'idec'; format = 'D'
    name = 'iflux_kron'; format = 'D'
    name = 'iflux_kron_err'; format = 'D'
    name = 'iflux_kron_flags'; format = 'L'
    name = 'imag_kron'; format = 'E'
    name = 'imag_kron_err'; format = 'E'
    name = 'imag_cmodel'; format = 'E'
    name = 'imag_cmodel_err'; format = 'E'
    name = 'iflux_cmodel_flags'; format = 'L'
    name = 'iflux_cmodel'; format = 'D'
    name = 'iflux_cmodel_err'; format = 'D'
    name = 'merge_measurement_i'; format = 'L'
    name = 'a_g'; format = 'E'
    name = 'a_r'; format = 'E'
    name = 'a_i'; format = 'E'
    name = 'a_z'; format = 'E'
    name = 'a_y'; format = 'E'
    name = 'gmag_forced_kron'; format = 'E'
    name = 'gmag_forced_kron_err'; format = 'E'
    name = 'gflux_forced_kron_flags'; format = 'L'
    name = 'rmag_forced_kron'; format = 'E'
    name = 'rmag_forced_kron_err'; format = 'E'
    name = 'rflux_forced_kron_flags'; format = 'L'
    name = 'imag_forced_kron'; format = 'E'
    name = 'imag_forced_kron_err'; format = 'E'
    name = 'iflux_forced_kron_flags'; format = 'L'
    name = 'zmag_forced_kron'; format = 'E'
    name = 'zmag_forced_kron_err'; format = 'E'
    name = 'zflux_forced_kron_flags'; format = 'L'
    name = 'ymag_forced_kron'; format = 'E'
    name = 'ymag_forced_kron_err'; format = 'E'
    name = 'yflux_forced_kron_flags'; format = 'L'
    name = 'gmag_forced_cmodel'; format = 'E'
    name = 'gmag_forced_cmodel_err'; format = 'E'
    name = 'gflux_forced_cmodel'; format = 'D'
    name = 'gflux_forced_cmodel_err'; format = 'D'
    name = 'gflux_forced_cmodel_flags'; format = 'L'
    name = 'rmag_forced_cmodel'; format = 'E'
    name = 'rmag_forced_cmodel_err'; format = 'E'
    name = 'rflux_forced_cmodel'; format = 'D'
    name = 'rflux_forced_cmodel_err'; format = 'D'
    name = 'rflux_forced_cmodel_flags'; format = 'L'
    name = 'imag_forced_cmodel'; format = 'E'
    name = 'imag_forced_cmodel_err'; format = 'E'
    name = 'iflux_forced_cmodel'; format = 'D'
    name = 'iflux_forced_cmodel_err'; format = 'D'
    name = 'iflux_forced_cmodel_flags'; format = 'L'
    name = 'zmag_forced_cmodel'; format = 'E'
    name = 'zmag_forced_cmodel_err'; format = 'E'
    name = 'zflux_forced_cmodel'; format = 'D'
    name = 'zflux_forced_cmodel_err'; format = 'D'
    name = 'zflux_forced_cmodel_flags'; format = 'L'
    name = 'ymag_forced_cmodel'; format = 'E'
    name = 'ymag_forced_cmodel_err'; format = 'E'
    name = 'yflux_forced_cmodel'; format = 'D'
    name = 'yflux_forced_cmodel_err'; format = 'D'
    name = 'yflux_forced_cmodel_flags'; format = 'L'
    name = 'ishape_hsm_regauss_e1'; format = 'E'
    name = 'ishape_hsm_regauss_e2'; format = 'E'
    name = 'ishape_hsm_regauss_sigma'; format = 'E'
    name = 'ishape_hsm_regauss_resolution'; format = 'E'
    name = 'ishape_hsm_regauss_derived_weight'; format = 'D'
    name = 'ishape_hsm_regauss_derived_sigma_e'; format = 'E'
    name = 'ishape_hsm_regauss_derived_rms_e'; format = 'D'
    name = 'ishape_hsm_regauss_derived_bias_m'; format = 'D'
    name = 'ishape_hsm_regauss_derived_bias_c1'; format = 'D'
    name = 'ishape_hsm_regauss_derived_bias_c2'; format = 'D'
    name = 'ishape_sdss_ixx'; format = 'E'
    name = 'ishape_sdss_iyy'; format = 'E'
    name = 'ishape_sdss_ixy'; format = 'E'
    name = 'ishape_sdss_psf_ixx'; format = 'E'
    name = 'ishape_sdss_psf_iyy'; format = 'E'
    name = 'ishape_sdss_psf_ixy'; format = 'E'
    name = 'weak_lensing_flag'; format = 'L'
    name = 'tract'; format = 'J'
    name = 'patch'; format = 'J'
    name = 'merge_peak_g'; format = 'L'
    name = 'merge_peak_r'; format = 'L'
    name = 'merge_peak_i'; format = 'L'
    name = 'merge_peak_z'; format = 'L'
    name = 'merge_peak_y'; format = 'L'
    name = 'gcountinputs'; format = 'I'
    name = 'rcountinputs'; format = 'I'
    name = 'icountinputs'; format = 'I'
    name = 'zcountinputs'; format = 'I'
    name = 'ycountinputs'; format = 'I'
    name = 'iflags_pixel_bright_object_center'; format = 'L'
    name = 'iflags_pixel_bright_object_any'; format = 'L'
    name = 'iblendedness_flags'; format = 'L'
    name = 'iblendedness_abs_flux'; format = 'E'
)
"""

"""
Indeed the color cuts will be applied in source_select function.

All the magnitudes required for these selections should be in the weak
lensing catalog files. In particular you need the iband cmodel
magnitude, and then the forced cmodel magnitude for all the colors.
Please take care to correct for the extinctions.

DataStore/S16A_v2.0/AEGIS_tracts/16821.fits

has the required columns.
"""

#filename="16821.fits"
#filename = "hsc-unblinded-Aug2017/WIDE12H_calibrated.fits"
#filename = "/project/plazas/WORK/HSC/weaklens_pipeline/DataStore/Calibrated_S16A_v2.0/XMM_calibrated.fits"

# S19A; different column name format: 'forced_z_cmodel_mag'
filename="/project/plazas/WORK/HSC/weaklens_pipeline/DataStore/S19A_v2.0/GAMA09H_tracts/10047_no_m.fits"
data = pf.open(filename)[1].data

#print (len(data))
CCParams = [-0.0248, 1.604, 2.276, -0.152] # Slope and intercept of lines in Eq. 23 and 21
#rlim = [0.5, 21, 28, -0.5, 2.]
rlim = [0.5, 21, 28, -0.7, 4.] # For zl < 0.4, Eq. 26-27 of appendix of Medezinski+18
#blim = [-3, -0.8, 22, 28, 0.5, 0.5]
blim = [-3, -0.8, 22, 28, 0.5, 0.5] # For zl < 0.4, Eq. 35 of appendix of Medezinski+18


#CCparams [-0.039063287566281635, 1.4441117764471074, -1.595881595881596, -0.16631566631566574]

red, blue, redOrBlueFlag = SelectSources (data, 'forced_z', 'forced_r', 'forced_i', 'forced_g', CCParams, rlim, blim,
        ext="_cmodel_mag")


newSources = data[redOrBlueFlag]

print ("Data, new sources: ", len(data), len(newSources))


#print ("RED: ", red)
#import ipdb; ipdb.set_trace()
#columns = [ 'z', 'r', 'g', 'i']
#cc.run_CCselection (data, columns, )

#print ("RA: ", np.min(data['ira']), np.max(data['ira']))
#print ("DEC: ", np.min(data['idec']), np.max(data['idec']))


# An ACT cluster in XMM: 37.9292, -4.8222 (approx)
# Exact: see Miyatake+19 (ACT-CL J0231.7-0452 2:31:43.63 −4:52:56.16)

#CCparams = [-0.06939717631541527, 2.1186868686868676, -1.5184544892337102, -0.5815448080058467]
#columns = [ 'z', 'r', 'g', 'i', 'ira', 'idec']
#center = [37.93179166666667, -4.882266666666667]
#scale = 3600 # deg to arcmin???

#red, blue, back, CCParams = cc.run_CCselection (data, columns, center, scale=scale, CClim = np.array([[-1.5, 2.5],[-3.5,
#    5.0]]), CCparams = CCparams)

#print (len(red), len(blue), len(back), len(CCParams), CCParams)

CClim = np.array([[-1.5, 2.5],[-3.5,5.0]])

singleFigure=False
pp = PdfPages("./test_cc_cuts.pdf")

#%% CC plot
plt.figure()
plotargs = { 'facecolors':'0.25','edgecolor':'0.25'}
obj1 = cc.CCplot(newSources,'forced_z', 'forced_r', 'forced_i', 'forced_g', pp, singleFigure=singleFigure,
        axis=CClim,plotpoints=True, CCparams=CCParams, 
                ext="_cmodel_mag", **plotargs)
plotargs = { 'facecolors':'r','edgecolor':'r'}
obj2 = cc.CCplot(red,'forced_z', 'forced_r', 'forced_i', 'forced_g', pp, singleFigure=singleFigure, axis=CClim,plotpoints=True,
                 ext="_cmodel_mag", CCparams=CCParams, **plotargs)
plotargs = { 'facecolors':'b','edgecolor':'b'}
obj3 = cc.CCplot(blue,'forced_z', 'forced_r', 'forced_i', 'forced_g', pp, singleFigure=singleFigure, axis=CClim,plotpoints=True,
                 ext="_cmodel_mag", CCparams=CCParams, **plotargs)
if singleFigure:
    pp.savefig()

pp.close()
