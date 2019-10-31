import ncdf2dict as nc
from pyEquilibrium.equilibrium import equilibrium
import f90nml
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.interpolate import interp1d
import sys

class Parameter:
    def __init__(self, value):
        self.value = value

    def set(self, value):
        self.value = value

    def __call__(self):
        return self.value


def fit(function, parameters, y, x=None):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        # return np.sum(abs(y - function(x)))
        return y - function(x)

    if x is None:
        x = np.arange(y.shape[0])
    p = [param() for param in parameters]
    # return optimize.minimize(f, p, method='Nelder-Mead')
    return optimize.least_squares(f, p)


def miller_rz(theta, kappa, delta, rcen, rmin):

    R = rcen + rmin * np.cos(theta + np.arcsin(delta)*np.sin(theta))
    Z = kappa * rmin * np.sin(theta)

    return R, Z


def miller_bp(theta):

    term1 = dpsidr_fit() / (kap_par() * R_par())

    term2 = np.sqrt(np.sin(theta + x_par() * np.sin(theta))**2 *
                    (1 + x_par() * np.cos(theta))**2 + (kap_par() * np.cos(theta))**2)

    term3 = np.cos(x_par() * np.sin(theta)) + shift_fit() * np.cos(theta)

    term4 = (s_kappa_fit() - s_delta_fit()*np.cos(theta) + (1 + s_kappa_fit()) * x_par() *
             np.cos(theta)) * np.sin(theta) * np.sin(theta + x_par() * np.sin(theta))
    Bp = term1 * term2 / (term3 + term4)

    return Bp


args = sys.argv[1:]

# Readin GEQDSK file name
gfile = args[0]

runname = gfile[:-7]

# Load in equilbrium
eq = equilibrium()
eq.load_geqdsk(gfile)


# Select surface in rho_psi coordinates
rho_psi = 0.5

psibnd = eq.psi_bnd

psi = rho_psi * psibnd

q = eq.qpsi(psi)

f = eq.fpol(psi)

bdy = eq.get_fluxsurface(psiN=1.0)
Rbdy = bdy.R
Zbdy = bdy.Z
amin = (max(Rbdy) - min(Rbdy)) / 2

print('Minor radius = {:.3f} m'.format(amin))

flxsur = eq.get_fluxsurface(psiN=rho_psi)

R = flxsur.R
Z = flxsur.Z
R_par = Parameter(R)

bpol = [eq.Bp(R[i], Z[i]) for i in range(len(R))]

btor = eq.fpol(psi) / R

bppts = np.squeeze(np.array(bpol))

bmag = np.sqrt(btor**2 + bppts**2)

Rmax = max(R)
Rmin = min(R)
Zmax = max(Z)
Zmin = min(Z)
Zmid = (Zmax+Zmin)/2

if Zmid > 1.e-4:
    print('#### Zmid != 0 ###')
    print('Zmid = {:.3f}'.format(Zmid))

Zind = np.argmax(Z)
Rupper = R[Zind]

Rmaj = (Rmax + Rmin) / 2
rmin = (Rmax - Rmin) / 2
rhoc = rmin/amin

r = R - Rmaj
z = Z

kappa = max([Zmax, -Zmin])/rmin

kap_par = Parameter(kappa)

delta = (Rmaj - Rupper)/rmin
x = np.arcsin(delta)
x_par = Parameter(x)

psi_high = psi * 1.000001
psi_low = psi / 1.000001

Rhigh, Zhigh = eq.get_fluxsurface(psiN=psi_high/psibnd)
Rlow, Zlow = eq.get_fluxsurface(psiN=psi_low/psibnd)


rhigh = (max(Rhigh) - min(Rhigh))/2
rlow = (max(Rlow) - min(Rlow))/2

dr = rhigh - rlow

q_high = eq.qpsi(psi_high)
q_low = eq.qpsi(psi_low)
dq = q_high - q_low

p = eq.pres(psi)
p_high = eq.pres(psi_high)
p_low = eq.pres(psi_low)
dp = p_high - p_low


shat = rmin/q * dq/dr
shat = shat.flatten()[0]

dpdr = dp/dr

theta_geo = np.arctan(z/r)
for i in range(len(theta_geo)):
    if r[i] < 0 and z[i] >= 0:
        theta_geo[i] += np.pi
    elif r[i] < 0 and z[i] < 0:
        theta_geo[i] = - np.pi + theta_geo[i]


theta = np.arcsin(z/(kappa*rmin))

nan = np.argwhere(np.isnan(theta))

for i in range(len(theta)):
    if R[i] < Rupper:
        if z[i] >= 0:
            theta[i] = np.pi - theta[i]
        elif z[i] < 0:
            theta[i] = -np.pi - theta[i]

bpol_mid = [eq.Bp(Rval, 0) for Rval in R]
btor_mid = [eq.Bt(Rval, 0)/Rval for Rval in R]

bpol_mid = np.squeeze(bpol_mid)
btor_mid = np.squeeze(btor_mid)

bmag_mid = np.sqrt(btor_mid**2 + bpol_mid**2)

Bcen = eq.Bt(Rmaj, Zmid)[0][0]
print('B0 = {:.3f}'.format(Bcen))
Rcen = Rmaj
print('R0 = {:.3f}'.format(Rcen))

# Assume half pressure comes from electrons
p_e = p/2

beta = p_e/Bcen**2 * 8 * np.pi * 1.0e-7

beta_prime = amin/p * dpdr * beta * 2

print("a/Lp             = {:.4f}".format(amin/p * dpdr))
print('beta             = {:.4f}'.format(beta))
print('beta_prime_input = {:.4f}'.format(beta_prime))


cdffile = runname+'.cdf'
data = nc.Dataset(cdffile)

psiN = data['rho_psi'][:]**2

con = np.argmin(abs(psiN - rho_psi))

# This doesn't match well with GACODE
#amin = data['amin'][0].data

s_kappa = data['TGLF_S_KAPPA'][con]
s_delta = data['TGLF_S_DELTA'][con]
dpsidr = data['dPsidrho'][con] / amin
shift = data['TGLF_DRMAJDX'][con]

print(' ')
print('Psi normalised = {:.3f}'.format(rho_psi))
print('Major radius   = {:.3f}'.format(Rmaj))
print('Minor radius   = {:.3f}'.format(rmin))
print('rhoc           = {:.3f}'.format(rhoc))

print('Elongation     = {:.3f}'.format(kappa))
print('Triangularity  = {:.3f}'.format(delta))



# Generate species dependant parameters
rho_a = np.flip(data['TGLF_RMIN'][:].data, axis=0)

# Temp gradients
tprims = np.flip(data['TGLF_RLTS_e'][:].data, axis=0)
tprim = np.interp(rhoc, rho_a, tprims)

# Density gradients
fprims = np.flip(data['TGLF_RLNS_e'][:].data, axis=0)
fprim = np.interp(rhoc, rho_a, fprims)

# Collisionality
vnuis = np.flip(data['vnui'][:].data, axis=0)
vnui= np.interp(rhoc, rho_a, vnuis)

vnues = np.flip(data['vnue'][:].data, axis=0)
vnue = np.interp(rhoc, rho_a, vnues)


### Print Output ###
print(' ')
print('Initial guess for gradients')
print('S_KAPPA        = {:.3f}'.format(s_kappa))
print('S_DELTA        = {:.3f}'.format(s_delta))

print('Shift          = {:.3f}'.format(shift))
print('dpsidr         = {:.3f}'.format(dpsidr))


Rmil, Zmil = miller_rz(theta, kappa, delta, Rmaj, rmin)

plt.plot(R, Z, linewidth=3, label='SCENE GEQDSK')
plt.plot(Rmil, Zmil, linewidth=3, label='Fit')
plt.xlabel('R(m)')
plt.ylabel('Z(m)')
plt.legend()
ax = plt.gca()
ax.set_aspect('equal')
plt.show()


s_kappa_fit = Parameter(s_kappa)

s_delta_fit = Parameter(s_delta)

shift_fit = Parameter(shift)

dpsidr_fit = Parameter(dpsidr)


kap_par = Parameter(kappa)

x_par = Parameter(x)


params = [s_kappa_fit, s_delta_fit, shift_fit,
          dpsidr_fit]

fits = fit(miller_bp, params, bppts, theta)

fit_params = fits.x
#fit_params = fits[0]

print(' ')
print('S_KAPPA fitted  = {:.3f}'.format(fit_params[0]))
print('S_DELTA fitted  = {:.3f}'.format(fit_params[1]))

print('Shift fitted    = {:.3f}'.format(fit_params[2]))
print('dpsidr fitted   = {:.3f}'.format(fit_params[3]))

akappri = fit_params[0] * kappa/(rhoc)

tripri = fit_params[1] * np.sqrt(1-delta**2) / (rhoc)

shift = fit_params[2]

print(' ')
print('GS2 Parameters')
print('!      psi_norm  = {:.4f}'.format(rho_psi))
print('rhoc             = {:.4f}'.format(rhoc))
print('rmaj             = {:.4f}'.format(Rmaj/amin))
print('r_geo            = {:.4f}'.format(Rmaj/amin))
print('qinp             = {:.4f}'.format(q))
print('shat             = {:.4f}'.format(shat))
print('shift            = {:.4f}'.format(shift))
print('akappa           = {:.4f}'.format(kappa))
print('akappri          = {:.4f}'.format(akappri))
print('tri              = {:.4f}'.format(delta))
print('tripri           = {:.4f}'.format(tripri))
print('s_hat_input      = {:.4f}'.format(shat))
print('beta             = {:.4f}'.format(beta))
print('beta_prime_input = {:.4f}'.format(beta_prime))
print('tprim            = {:.4f}'.format(tprim))
print('fprim            = {:.4f}'.format(fprim))

plt.plot(theta/np.pi, miller_bp(theta), 'gd', label='Python Miller fit')
plt.plot(theta/np.pi, bppts, 'ro', label='SCENE Actual')
plt.xlabel(r'$\theta / \pi$')
plt.ylabel(r'$B_p$')
plt.legend()
plt.savefig('scene_miller.png')
plt.show()


# shift_fit.set(-0.8)
# bmag_high = np.sqrt(btor**2 + miller_bp(theta))

# shift_fit.set(-0.3)
# bmag_low = np.sqrt(btor**2 + miller_bp(theta))


# plt.plot(theta/np.pi, bmag, 'x', label='shift = 0.51')
# plt.plot(theta/np.pi, bmag_high, '^', label='shift = 0.8')
# plt.plot(theta/np.pi, bmag_low, '^', label='shift = 0.3')
# plt.xlabel(r"$\theta (\pi) $")
# plt.ylabel(r"$|B|$")
# plt.legend()
# plt.title(r"Magnitude of B for $\rho_\psi = {}$".format(rho_psi))
# plt.show()


# GS2 namelist template
nml = f90nml.read("/home/userfs/b/bbp501/SCENEv2/scripts/gs2_template.in")


# Geometry terms
nml['theta_grid_parameters']['rhoc'] = rhoc

nml['theta_grid_parameters']['akappa'] = kappa
nml['theta_grid_parameters']['akappri'] = akappri

nml['theta_grid_parameters']['tri'] = delta
nml['theta_grid_parameters']['tripri'] = tripri

nml['theta_grid_parameters']['qinp'] = q[()]
nml['theta_grid_parameters']['shat'] = shat

nml['theta_grid_parameters']['r_geo'] = Rmaj/amin
nml['theta_grid_parameters']['rmaj'] = Rmaj/amin

nml['theta_grid_parameters']['shift'] = shift

nml['theta_grid_eik_knobs']['s_hat_input'] = shat
nml['theta_grid_eik_knobs']['beta_prime_input'] = beta_prime

nml['parameters']['beta'] = beta


# Species terms

# Ions first
nml['species_parameters_1']['tprim'] = tprim
nml['species_parameters_1']['fprim'] = fprim
nml['species_parameters_1']['vnewk'] = vnui


# Electrons
nml['species_parameters_2']['tprim'] = tprim
nml['species_parameters_2']['fprim'] = fprim
nml['species_parameters_2']['vnewk'] = vnue


nmlfile = runname+'_gs2.in'
nml.write(nmlfile, force=True)

