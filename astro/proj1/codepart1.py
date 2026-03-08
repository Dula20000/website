import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


G    = 6.674e-11
h    = 6.626e-34
m_e  = 9.109e-31
m_H  = 1.673e-27
M_sun = 1.989e30


#   P = K * rho^(5/3)       using rho = 2m_H*n_e for Non Relativistic Degenerate
#^^^^^check dulran.com/astro/proj1/image1.heic for derivation
#   P = (h^2/5m_e) * (3/8pi)^(2/3) * n_e^(5/3)
#   rho = 2*m_H*n_e
#   n_e = rho/(2*m_H)
#   K = (h^2/5m_e) * (3/8pi)^(2/3) * (1/(2*m_H))^(5/3)
K_si = (h**2 / (5*m_e)) * (3/(8*np.pi))**(2/3) * (1/(2*m_H))**(5/3)



#Code units
rho0  = 1e9              # We chose a typical density for a whitedwardf
R0    = np.sqrt(K_si / (G * rho0**(1/3))) #check dulran.com/astro/proj1/image2.heic for derivation
M0    = rho0 * R0**3
P0    = K_si * rho0**(5/3)
print(f"  rho0 = {rho0:.3e} kg/m^3")
print(f"  R0   = {R0/1e6:.3f} Mm  ({R0:.3e} m)")
print(f"  M0   = {M0/M_sun:.4f} M_sun  ({M0:.3e} kg)")
print(f"  P0   = {P0:.3e} Pa")
print()


def P_from_rho(rho_t):
    return rho_t**(5/3)

def rho_from_P(P_t):
    if P_t <= 0:
        return 0.0
    return P_t**(3/5)
##^^^^check dulran.com/astro/proj1/image3.heic for derivation


def integrate_WD(rho_c_tilde, dr=1e-4, return_profile=False): #to integrate stellar struc from center to surface
    #returns the arrays of r, m, and P from center to surfaces when return profile is on
    #normally returns just the full mass and and surface radius
    P_c = P_from_rho(rho_c_tilde)
    r = dr #to avoid zero division error
    m = (4/3) * np.pi * r**3 * rho_c_tilde
    P = P_c

    if return_profile:
        rs, ms, Ps = [r], [m], [P]

    max_steps = 200000
    for i in range(max_steps):# iterate 200k times until convergance
        rho = rho_from_P(P)

        dPdr = -rho * m / r**2          # hydrostatic equilibrium in code units
        dmdr = 4 * np.pi * r**2 * rho   # mass conservation in code units

        P_new = P + dPdr * dr #iterate to get p new
        m_new = m + dmdr * dr # iterate to get m new
        r_new = r + dr # r new by iteration

        if P_new <= 0:
            frac  = P / (P - P_new)     #if we have gone past P=0 we use this to get the fraction of the step to get to P=0
            #check dulran.com/astro/proj1/image4.heic
            r_surf = r + frac * dr
            m_surf = m + dmdr * frac * dr
            if return_profile:
                return r_surf, m_surf, np.array(rs), np.array(ms), np.array(Ps)
            return r_surf, m_surf

        P, m, r = P_new, m_new, r_new #updates the P m r to the iteration

        if return_profile:
            rs.append(r); ms.append(m); Ps.append(P)

    # If we never hit P=0
    raise RuntimeError(f"Integration did not end for rho_c={rho_c_tilde:.3e}")


#Sweeping through masses from 0.2 to 1.4 Mnught and solving to get the needed density
#solve using brenq solver
def find_rho_c(M_target_tilde, dr=1e-4,
               rho_lo=1e-3, rho_hi=1e4):#goes between rho_lo and rho_hi to get the rho for which you get the target M
    def residual(rho_c):
        _, m_surf = integrate_WD(rho_c, dr=dr)#uses the integrate function from earlier to get only the Mass and gets the difference from the taget mass
        return m_surf - M_target_tilde


    rho_c_sol = brentq(residual, rho_lo, rho_hi, xtol=1e-6, rtol=1e-6)#uses brentq to minimize the m_surf - M target (residual) such that we havr the rho needed for the Mass taget.
    return rho_c_sol


###getting the densities from masses
masses_sun = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                        0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4])

masses_tilde  = masses_sun * M_sun / M0# convert to code units
dr            = 5e-5
print("central densities needed for taget mass and Radius","\n\n")
radii_tilde = []
rho_c_list  = []

for M_t in masses_tilde:
    rho_c = find_rho_c(M_t, dr=dr)
    R_t, _ = integrate_WD(rho_c, dr=dr)##gets the radius from the density using the integrate function
    radii_tilde.append(R_t)
    rho_c_list.append(rho_c)
    print(f"  M = {M_t*M0/M_sun:.2f} M_sun  |  "
          f"rho_c = {rho_c*rho0:.2e} kg/m^3  |  "
          f"R = {R_t*R0/1e6:.2f} Mega m  ({R_t*R0/6.371e6:.3f} R_earth)")

radii_tilde = np.array(radii_tilde)
radii_Mm    = radii_tilde * R0 / 1e6
radii_Rearth = radii_tilde * R0 / 6.371e6


print("\n\nConvergence Check as Steps size decreases for M = 0.6")
M_test = 0.6 * M_sun / M0
for dr_test in [1e-4, 5e-5, 2.5e-5, 2e-5]:#1/200000 = 2e-5 , 200k was the step limit
    rho_c = find_rho_c(M_test, dr=dr_test, rho_lo=0.1, rho_hi=1e4)
    R_t, _ = integrate_WD(rho_c, dr=dr_test)
    print(f"  dr={dr_test:.1e}  =>  R = {R_t*R0/1e6:.4f} Mm")


#Plotting

# Plot 1: Mass and Radius relation vs R ~= M^-1/3
fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(masses_sun, radii_Rearth, 'o-', color='steelblue',
        linewidth=2, markersize=6, label='Numerical (NR degenerate EOS)')

M_ref  = 0.2 #choose M = 0.2 as the anchor point to plot the R = M^-1/3 as it doesn't have a reference condition/ initial conditon. The relation is only a shape you need an initial cond
R_ref  = radii_Rearth[masses_sun == M_ref][0]
M_anal = np.linspace(0.18, 1.42, 200)
R_anal = R_ref * (M_anal / M_ref)**(-1/3)

ax.plot(M_anal, R_anal, '--', color='tomato',
        linewidth=1.8, label=r'Analytical scaling $R \propto M^{-1/3}$')

ax.set_xlabel(r'Mass  [$M_\odot$]', fontsize=13)
ax.set_ylabel(r'Radius  [$R_\oplus$]', fontsize=13)
ax.set_title('White Dwarf Mass–Radius Relation', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
plt.close()
