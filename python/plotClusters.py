import numpy as np
import matplotlib.pyplot as plt 
import healpy as hp 
import astropy.io.fits as pf

dirRoot = "/project/plazas/WORK/HSC/sz_clusters_mass_cal/"

def selectGalaxyCatalog(d, star_mask = "S18A", d_sm = None): 
    print("Number of galaxies before any selection: %d" % len(d))  
    # apply FDFC mask 
    m = hp.read_map(dirRoot + "s19a_fdfc_hp_contarea_izy-gt-5_trimmed.fits", nest = True, dtype = bool) 
    indices_map = np.where(m)[0] 
    nside = hp.get_nside(m) 
    phi = np.radians(d["RA"]) 
    theta = np.radians(90. - d["DEC"]) 
    indices_obj = hp.ang2pix(nside, theta, phi, nest = True) 
    sel = np.in1d(indices_obj, indices_map) 
    print("number of galaxies after all selections %d" % sel.sum()) 
    return sel  

# Planck 20215 SZ Union cluster catalog
hdulist = pf.open(dirRoot+ "HFI_PCCS_SZ-union_R2.08.fits")

data = hdulist[1].data 
sel = selectGalaxyCatalog(data)
new_data = data[sel]
redshiftSelClusters = new_data['REDSHIFT']
print (redshiftSelClusters)
masSelClusters = data[sel]['MSZ']
print (np.max(data['MSZ']), np.min(data['MSZ']), masSelClusters)

# print some properties
print (" number NAME RA[deg] DEC[deg] REDSHIFT MSZ (10^14 M_sun)")
planckNAME, planckRA, planckDEC, planckREDSHIFT, planckMSZ = [], [], [], [], []
numberPlackClusters = 0
for i, row in enumerate(new_data):
    print (i+1, row['NAME'], row['RA'], row['DEC'], row['REDSHIFT'], row['MSZ'])
    planckNAME.append(row['NAME'])
    planckRA.append(row['RA'])
    planckDEC.append(row['DEC'])
    planckREDSHIFT.append(row['REDSHIFT'])
    planckMSZ.append(row['MSZ'])
    numberPlackClusters+=1
"""
# Plot
sc = plt.scatter(new_data['RA'], new_data['DEC'], s=10*masSelClusters, c=redshiftSelClusters, vmin=0.05,
                 vmax=0.51)
plt.scatter(data[~sel]['RA'], data[~sel]['DEC'], s=data[~sel]['MSZ'], c='k')
#plt.axes().set_aspect(1./2)
plt.xlim ([370,-10]) 
plt.ylim ([-10,50])
plt.xlabel('RA [deg]')
plt.ylabel('DEC [deg]')
plt.grid(True)
plt.colorbar(sc)
plt.show()
"""


### --- Read camira S20 cluster catalog 
camiraData = np.genfromtxt (dirRoot + "camira_s20a_wide_v1.dat")
ra_cam, dec_cam, zcl_cam, mlog_cam  = camiraData[:,0], camiraData[:,1], camiraData[:,2], camiraData[:,4]

# with cut=0.1, all have only 1 candidate except ID=0 (and ID=16, that has zplanck=-1)

"""
339.06524 1.329462 0.1856 
225.2846 42.407843 0.2736
216.814982 44.10775 0.5043
213.431212 43.644742 0.1
351.935008 0.868511 0.2592
354.415554 0.271375 0.2643
355.920075 0.312438 0.284 
0.912199 1.825302 0.3446 
33.66538 -4.53105 0.1463 
37.921574 -4.882613 0.1858  
131.365685 3.460808 0.3222
134.474796 3.176485 0.188  
140.531878 3.766339 0.252 
139.038475 -0.404508 0.304 
145.102438 2.477641 0.167 
143.805679 0.811196 0.3532
146.506134 0.470659 1.166
180.105635 3.347094 0.1322
185.70092 2.10474 0.2245


All except 0, 7, 16 have < 0.2 in all 3 (ra, dec, z)

0 —> cut: 0.15 for all 3
16 —> cut: 0.05 for RA, DEC
7—> Two candidates, but not too close in redshift: 
Plack:  0.9450593131603874 2.0536454341381765 0.0924

0.912199 1.825302 0.3446
0.914498 2.235773 0.4224
"""
# File made by hand with data from above
centersCamira = np.genfromtxt (dirRoot + "centers_from_camira.dat")
candidateRA, candidateDEC, candidateZ = centersCamira[:,0], centersCamira[:,1], centersCamira[:,2] 
                                         
cut = 0.1
counter=1
for (r, d, z, m) in zip(planckRA, planckDEC, planckREDSHIFT, planckMSZ):
    diffRA, diffDEC, diffZ = np.fabs(ra_cam - r), np.fabs(dec_cam - d), np.fabs(zcl_cam - z)
    
    mask = (diffRA <= cut) & (diffDEC <= cut) #& ((diffZ <= cut))

    print (" Number: ", counter)
    print ("Plack: ", r, d, z, m)
    print ("Candidates Camira: ")
    for (r2, d2, z2, m2) in zip(ra_cam[mask], dec_cam[mask], zcl_cam[mask], mlog_cam[mask]):
        print (r2, d2, z2, m2, 10**m2)
    counter+=1

##### Print PLANCK and CAMIRA Centers at the same time
print (" number        NAME       RA[deg]   DEC[deg]   REDSHIFT   MSZ(10^14 M_sun)  RA_CAMIRA[deg]      DEC_CAMIRA[deg]     REDSHIT_CAMIRA ")

for (name, ra_p, dec_p, z_p, m_p, ra_camira, dec_camira, z_camira) in zip (planckNAME, planckRA, planckDEC,
        planckREDSHIFT, planckMSZ, candidateRA, candidateDEC, candidateZ):
    #print (name, ra_p, dec_p, z_p, m_p, ra_camira, dec_camira, z_camira)
    print (ra_p, dec_p, z_p, ra_camira, dec_camira, z_camira)

##############  Try with DR8 redmapper SDDS
"""
ColDefs(
    name = 'ID'; format = 'J'
    name = 'NAME'; format = '20A'
    name = 'RA'; format = 'D'
    name = 'DEC'; format = 'D'
    name = 'Z_LAMBDA'; format = 'E'
    name = 'Z_LAMBDA_ERR'; format = 'E'
    name = 'LAMBDA'; format = 'E'
    name = 'LAMBDA_ERR'; format = 'E'
    name = 'S'; format = 'E'
    name = 'Z_SPEC'; format = 'E'
    name = 'OBJID'; format = 'K'
    name = 'IMAG'; format = 'E'
    name = 'IMAG_ERR'; format = 'E'
    name = 'MODEL_MAG_U'; format = 'E'
    name = 'MODEL_MAGERR_U'; format = 'E'
    name = 'MODEL_MAG_G'; format = 'E'
    name = 'MODEL_MAGERR_G'; format = 'E'
    name = 'MODEL_MAG_R'; format = 'E'
    name = 'MODEL_MAGERR_R'; format = 'E'
    name = 'MODEL_MAG_I'; format = 'E'
    name = 'MODEL_MAGERR_I'; format = 'E'
    name = 'MODEL_MAG_Z'; format = 'E'
    name = 'MODEL_MAGERR_Z'; format = 'E'
    name = 'ILUM'; format = 'E'
    name = 'P_CEN'; format = '5E'
    name = 'RA_CEN'; format = '5E'
    name = 'DEC_CEN'; format = '5E'
    name = 'ID_CEN'; format = '5K'
    name = 'PZBINS'; format = '21E'
    name = 'PZ'; format = '21E'
)
"""

# From http://risa.stanford.edu/redmapper/
hdulist = pf.open(dirRoot + "redmapper_dr8_public_v6.3_catalog.fits.gz")
dataSDSS = hdulist[1].data
ra_sdss, dec_sdss, zcl_sdss = dataSDSS['RA'], dataSDSS['DEC'], dataSDSS['Z_LAMBDA']

#import ipdb; ipdb.set_trace()

print ("ra_sdss", ra_sdss, ra_sdss.shape)
counter=1
for (r, d, z, m) in zip(planckRA, planckDEC, planckREDSHIFT, planckMSZ):
    #print (r, d, z, m)
    diffRA, diffDEC, diffZ = np.fabs(ra_sdss - r), np.fabs(dec_sdss - d), np.fabs(zcl_sdss - z)
    mask = (diffRA <= cut) & (diffDEC <= cut) #& ((diffZ <= cut))

    print ("Number: ", counter)
    print ("Plack: ", r, d, z, m)
    print ("Candidates from SDSS: ")
    for (r2, d2, z2) in zip(ra_sdss[mask], dec_sdss[mask], zcl_sdss[mask]): # , mlog_cam[mask]):
        print (r2, d2, z2,) #m2, 10**m2)
    counter+=1
