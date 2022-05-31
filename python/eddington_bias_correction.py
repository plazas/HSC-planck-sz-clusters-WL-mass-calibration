import matplotlib.pyplot as plt
import numpy as np

from colossus.cosmology import cosmology
cosmology.setCosmology('planck15')

from colossus.lss import mass_function

from scipy.interpolate import UnivariateSpline

from matplotlib.backends.backend_pdf import PdfPages



def get_alpha(M, z):
    M_array = 10**np.arange(11.0, 15.5, 0.1)
    mfunc = mass_function.massFunction(M_array, z, mdef = '200m', model = 'tinker08', q_out = 'dndlnM')
    # Interpolate
    m_func_interp = UnivariateSpline(M_array, mfunc, s=0, k=3)
    # Define new function
    def dn_dM (x):
        return m_func_interp(x)*(1./x)

    dn_dM_interp = UnivariateSpline(M_array, dn_dM(M_array), s=0, k=3)
    dn_dM_derivative = dn_dM_interp.derivative(n=1)

    def dn2_d2lnM (x):
        return dn_dM_derivative(x)*x**2

    def alpha(x):
        return dn2_d2lnM(x)/m_func_interp(x)

    return alpha(M)


# Read Planck Masses and redshifts
mass_planck=np.array([1.9290442, 5.8693595, 5.225707, 2.0749116, 4.9161544, 7.3293004, 5.1834126, 2.2387743,
    4.3063836, 3.961237, 5.271713, 4.5982094,  5.781798, 4.867756, 3.43227, 5.4773803,  5.538075, 4.098035 ])*10**14

redshift = np.array([0.0594, 0.2917, 0.5020798, 0.089, 0.302061 , 0.2779,  0.27, 0.0924, 0.1393, 0.1843 ,
    0.3269, 0.203, 0.2701, 0.33240458, 0.1508, 0.35682267, 0.1339, 0.22898613])


def get_corrected_mass(Mobs, zobs, sigma_ln_M  = 0.3):
    temp = np.log(Mobs) + 0.5*get_alpha(Mobs, zobs)*sigma_ln_M**2
    return np.exp(temp)


for m, z in zip(mass_planck, redshift):
    mcorr = get_corrected_mass (m, z)
    factor = 1e14
    print (f"{m/factor}, {mcorr/factor}")


pp = PdfPages("mass_function.pdf")

# https://bdiemer.bitbucket.io/colossus/_static/tutorial_lss_mass_function.html
M = 10**np.arange(11.0, 15.5, 0.1)
print (M)

z = 0.05
mfunc = mass_function.massFunction(M, z, mdef = '200m', model = 'tinker08', q_out = 'dndlnM') #, q_out = 'dndlnM')
print (mfunc)


plt.figure()
plt.xlabel('M200m')
#plt.ylabel('dn/dln(M)')
plt.ylabel('M^2/rho_0 dn/d(M)')
plt.loglog()
plt.xlim(1E11, 4E15)
plt.ylim(1E-7, 1E-1)
plt.plot(M, mfunc, '-', label = f'z = {z}')
plt.legend()
pp.savefig()


# Interpolate
m_func_interp = UnivariateSpline(M, mfunc, s=0, k=3)

# Define new function
def dn_dM (M):
    return m_func_interp(M)*(1./M)

dn_dM_interp = UnivariateSpline(M, dn_dM(M), s=0, k=3)
dn_dM_derivative = dn_dM_interp.derivative(n=1)

def dn2_d2lnM (M):
    return -dn_dM_derivative(M)*M**2

def alpha(M):
    return dn2_d2lnM(M)/m_func_interp(M)

print ("alpha: ", alpha(M))

plt.figure()
plt.xlabel('M200m')
plt.ylabel('alpha(M)')
#plt.loglog()
plt.semilogx()
plt.xlim(1E11, 4E15)
plt.ylim(0, -5)
plt.gca().invert_yaxis()
plt.plot(M[:-2], -alpha(M)[:-2], '-', label = f'z = {z}')
plt.legend()
pp.savefig()

"""

