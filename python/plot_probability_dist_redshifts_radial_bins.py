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

for row in data:
    rr, zz, pofz = row
    rr=10**rr
    pofz_radial_bins[rr]['z'].append(zz)
    pofz_radial_bins[rr]['p_of_z'].append(pofz)


#assert(len(rbins) == len(pofz))
pp = PdfPages ("pofz_radial_bins.pdf")

fig = plt.figure(figsize=(15,15))
for i, rr in enumerate(pofz_radial_bins):
    #fig = plt.figure(figsize=(10,10))
    fig.add_subplot(4,5,i+1)
    z, pofz = pofz_radial_bins[rr]['z'], pofz_radial_bins[rr]['p_of_z']
    plt.plot(z, pofz, '-', label=f"rbin: {rr:.3f} \n Npairs: {npairs[i]}")
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.xlabel('z')
    plt.ylabel('p(z)')
    plt.xlim([0,4])
    plt.suptitle('P-cut', fontsize=15)

plt.tight_layout()
pp.savefig(fig)


#fig = plt.figure(figsize=(20,10))
#for r, p in zip(rbins, pofz):
#    plt.plot(r, p, '-.')
#    plt.xlabel('radial bin')
#    plt.ylabel('P(z)')   
#pp.savefig(fig)
pp.close()
