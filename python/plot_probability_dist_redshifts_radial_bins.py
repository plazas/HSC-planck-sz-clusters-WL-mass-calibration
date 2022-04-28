import numpy as np
import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as plt

# Can you stack the P(z)s of sources in radial bins around the clusters
# for samples selected by P(z) cuts and those by color cuts? The boost
# factor may not capture the problems if obscuration is an issue.

# You will have to setup an array of size rbins*Pzbins, then update this as soon as you find a lens source
# pair. ---> surhud. he implemented this in his code.


# rr, zz, pz_rbins[write_element]

data = np.genfromtxt(sys.argv[1], dtype=np.float64)

dsigma = np.genfromtxt(sys.argv[2], dtype=np.float64)

npairs = dsigma[:,13]

# For every radial bin, plot P(z)
# there’s 20 bins (the first column repeats itself every 20th entry, a 100 times). I’m stacking the P(z) in
# every radial bin along those 100

rbins = 10**(np.mean(data[:,0].reshape((100,20)), axis=0))
pofz = data[:,2].reshape((100,20))

print ("rbins: ", rbins)

pofz_radial_bins = {}
radial_bins=[]

#for rr in rbins:
#    key=f"{rr}"
#    print ("rr", rr)
#    pofz_radial_bins[key]={'z':[], 'p_of_z':[]}

for row in data:
    rr, _, _ = row
    pofz_radial_bins[10**rr]={'z':[], 'p_of_z':[]}

print ("pofz_radial_bins: ", pofz_radial_bins)

mean_pofz = []
for row in data:
    rr, zz, pofz = row
    mean_pofz
    rr=10**rr
    pofz_radial_bins[rr]['z'].append(zz)
    pofz_radial_bins[rr]['p_of_z'].append(pofz)

mean_pofz = []
for rr in pofz_radial_bins:
    pofz_array = pofz_radial_bins[rr]['p_of_z']
    mean_pofz.append(pofz_array)

mean_pofz = np.array(mean_pofz)
print ("before mean: ", mean_pofz)
mean_pofz = np.mean(mean_pofz, axis=0, dtype=np.float64)
print ("HOLA mean: ", mean_pofz)

#for x in mean_pofz:
#    print (f"{x}")

#assert(len(rbins) == len(pofz))
pp = PdfPages ("pofz_radial_bins.pdf")

fig = plt.figure(figsize=(15,15))
for i, rr in enumerate(pofz_radial_bins):
    #fig = plt.figure(figsize=(10,10))
    #fig.add_subplot(1,1,1)#4,5,i+1)
    z, pofz = pofz_radial_bins[rr]['z'], pofz_radial_bins[rr]['p_of_z']
    if i == 0:
        print (rr, pofz , mean_pofz, (pofz - mean_pofz)/mean_pofz)
    plt.plot(z, pofz, '-', label=f"rbin: {rr:.3f} \n Npairs: {npairs[i]}")
    plt.plot(z, mean_pofz, 'k-', label=f"Mean")
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.xlabel('z')
    plt.ylabel('p(z)')
    plt.xlim([1e-5,1])
    plt.xlim([0,1])
    plt.suptitle('P-cut', fontsize=15)

plt.tight_layout()
pp.savefig(fig)

fig = plt.figure(figsize=(15,15))
for i, rr in enumerate(pofz_radial_bins):
    #fig = plt.figure(figsize=(10,10))
    #fig.add_subplot(1,1,1)#4,5,i+1)
    z, pofz = pofz_radial_bins[rr]['z'], pofz_radial_bins[rr]['p_of_z']
    if i == 0:
        print (rr, pofz , mean_pofz, (pofz - mean_pofz)/mean_pofz)
    # only z<1
    z = np.array(z)
    mask = z < 1.0
    delta = (pofz - mean_pofz)/mean_pofz
    delta = delta[mask]
    z = z[mask]
    plt.plot(z, delta, '-', label=f"rbin: {rr:.3f} \n Max: {np.max(delta)}")
    #plt.plot(z, mean_pofz, 'k-', label=f"Mean")
    plt.yscale('linear')
    plt.legend(loc='upper right')
    plt.xlabel('z')
    plt.ylabel('p(z) - mean_pofz / mean_pofz')
    plt.ylim([-1, 1])
    plt.xlim([0,1])
    plt.suptitle('P-cut', fontsize=20)

plt.tight_layout()
pp.savefig(fig)




#fig = plt.figure(figsize=(20,10))
#for r, p in zip(rbins, pofz):
#    plt.plot(r, p, '-.')
#    plt.xlabel('radial bin')
#    plt.ylabel('P(z)')   
#pp.savefig(fig)
pp.close()
