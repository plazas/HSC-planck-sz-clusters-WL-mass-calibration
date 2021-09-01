import glob
import astropy.io.fits as pf
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

## Create files iwht random points to calculate weak lensing signal around them


files=glob.glob("/project/plazas/WORK/HSC/weaklens_pipeline/DataStore/xiangchong.li/S19ACatalogs/v2/catalog_others/*rand*")

downscale=80000 #100000
allRa=[]
allDec=[]
for fieldFile in files:
    data = pf.open(fieldFile)[1].data[::downscale]
    print (f"{fieldFile}")
    print (len(data))
    mask1 = (data["g_mask_s18a_bright_objectcenter"] & data["r_mask_s18a_bright_objectcenter"] &
             data["i_mask_s18a_bright_objectcenter"] & data["z_mask_s18a_bright_objectcenter"] &
             data["y_mask_s18a_bright_objectcenter"])
    mask2 = (data["g_pixelflags_bright_objectcenter"] & data["i_pixelflags_bright_objectcenter"] &
             data["r_pixelflags_bright_objectcenter"] & data["z_pixelflags_bright_objectcenter"] &
             data["y_pixelflags_bright_objectcenter"])
    mask3 = (data["g_mask_brightstar_any"] & data["r_mask_brightstar_any"] &
             data["i_mask_brightstar_any"] & data["z_mask_brightstar_any"] &
             data["y_mask_brightstar_any"])
    finalMask = ~(mask1 & mask2 & mask3)
    #print (finalMask)
    data = data[finalMask]
    print ("After mask: ", len(data))
    #print (data['ra'])
    #stop
    allRa.append(data['ra'])
    allDec.append(data['dec'])

# read redshits from 19 clusters
fileClusters="/project/plazas/WORK/HSC/weaklens_pipeline/DataStore/centers_from_camira.dat"
redshiftClusters = np.genfromtxt(fileClusters)[:,2]
print ("redshiftClusters", redshiftClusters)
# Plot
pp=PdfPages("randoms_PSZ2_hsc_s19a.pdf")

#fig = plt.figure(figsize=(20,10))

raVec, decVec, zVec = [],[],[]
for (ra, dec, fileName) in zip(allRa, allDec, files):
    fig = plt.figure(figsize=(20,10))
    fieldName = fileName.split('/')[-1].split('.')[0]
    #redshiftRandoms = np.random.random(len(ra))
    redshiftRandoms = np.random.choice(redshiftClusters, size=len(ra))
    sc = plt.scatter(ra, dec, s=4, c=redshiftRandoms, vmin=0.001,
                     vmax=np.max(redshiftClusters), cmap='viridis')
    for (r,d,z) in zip(ra, dec, redshiftRandoms):
        raVec.append(r)
        decVec.append(d)
        zVec.append(z)
    #plt.scatter(data[~sel]['RA'], data[~sel]['DEC'], s=data[~sel]['MSZ'], c='k')
    #plt.axes().set_aspect(1./2)
    #plt.xlim ([370,-10]) 
    #plt.ylim ([-10,50])
    plt.xlabel('RA [deg]')
    plt.ylabel('DEC [deg]')
    plt.grid(True)
    plt.suptitle(f'{fieldName}')
    plt.colorbar(sc, label='Redshift')
    pp.savefig(fig)

print ("Final number of randoms: ", len(raVec))
raVec=np.array(raVec[:1900])
decVec=np.array(decVec[:1900])
zVec=np.array(zVec[:1900])
assert (len(raVec) == len(decVec))
assert (len(decVec) == len(zVec))
weight = np.ones(len(raVec))

print (raVec, decVec, zVec)
print ("Types: ", type(raVec), type(raVec[1]))
print ("Types: ", type(decVec), type(decVec[1]))
print ("Types: ", type(zVec), type(zVec[1]))
print ("Types: ", type(weight), type(weight[1]))

outputData = np.array([raVec, decVec, zVec, weight]).astype(np.float64)
print ("Types: ", type(outputData))
print (outputData)
print (outputData.shape)
outputDataTranspose = np.transpose(outputData)
# Save singel file
fileNameRandoms = "/project/plazas/WORK/HSC/weaklens_pipeline/DataStore/random_points_s19a_test.dat"
np.savetxt(fileNameRandoms, np.transpose(outputData))

# Split the 1900 random points into 100 files of 19 points each. Save the 100 files in disk. 

part = np.split(outputDataTranspose, 100)
root="/project/plazas/WORK/HSC/weaklens_pipeline/DataStore/randoms_S19A/"
for i, randomRealization in enumerate(part):
    fileName=root+f"random_points_s19a_{i}.dat"
    np.savetxt(fileName, randomRealization)

#pp.savefig()
pp.close()
#plt.show()
