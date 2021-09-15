import numpy as np
import astropy.units as u
from astropy import constants as const
from scipy.special import lambertw

# ---- General classes for coronal models ----

class density_model():
    @u.quantity_input
    def __init__(self, r: u.km, model, nfold=1):
        self.model = model
        self.nfold = nfold
        self.r = r
        self.Ne = nfold * self.density(r, model)
    
    def density(self, r, model):
        r"""
        Return the electron density profile from different empirical coronal models.
        """

        if model == 'Newkirk':
            N0 = 4.2 * 10**4
            Ne = N0 * 10**(4.32*u.R_sun/r)
        elif model == 'Baumbach':  # Baumbach-Allen Density model (From Aswaden page 82)
            Ne = (2.99*(r/u.R_sun)**(-16) + 1.55*(r/u.R_sun)**(-6) + 0.036*(r/u.R_sun)**(-1.5))*10**8
        elif model == 'Saito':  # Saito Density model
            Ne = (3.09*(r/u.R_sun)**(-16) + 1.58*(r/u.R_sun)**(-6) + 0.0251*(r/u.R_sun)**(-2.5))*10**8
        elif model == 'Leblanc':  # Leblanch Density model
            Ne = (0.8*10**8*(r/u.R_sun)**(-6) + 0.41*10**7*(r/u.R_sun)**(-4) + 0.33*10**6*(r/u.R_sun)**(-2))
        else:
            raise ValueError("Non supported density model.")

        # All the models should give a density in cm^-3 before returning.
        Ne = ((Ne * u.cm**-3).decompose()).to(u.cm**-3)

        return Ne

@u.quantity_input
class magnetic_model():
    @u.quantity_input
    def __init__(self, r: u.km, model, topo,
                 B0: u.gauss, usw: (u.km / u.second) = None, alpha=1, ):
        self.model = model
        self.topo = topo
        self.B0 = B0
        self.r = r
        self.usw = usw
        self.B = self.magnetic_field(r, model, B0, alpha)
        self.br_angle = self.configuration(topo)

    def configuration(self, topo):
        if topo == 'Radial':
            br_angle = np.zeros_like(self.r.data) * u.rad
        if topo == 'Parker spiral(*)':
            sun_omega = 14.713 * (u.deg / u.day)
            usw = self.usw
            if usw is not None:
                br_angle = (sun_omega / usw) * self.r
            else:
                br_angle = None

        return br_angle

    def magnetic_field(self, r, model, B0, alpha):
        r"""
        Return the total magnetic field profile from coronal magnetic
        field models. Only quiet sun models are considered with
        recomended photosheric magnetic field of B0 = 2.2 gauss.
        The 1/r^2 model is recomended.
        """
        if model == '1/r^1':
            B = B0 * (1*u.R_sun/r)**1
        elif model == '1/r^2':
            B = B0 * (1*u.R_sun/r)**2
        elif model == '1/r^3':
            B = B0 * (1*u.R_sun/r)**3
        elif model == '1/r^a':
            B = B0 * (1*u.R_sun/r)**alpha
        else:
            raise ValueError("Non supported magnetic field model.")

        B = (B.decompose()).to(u.gauss)

        return B

@u.quantity_input
class solarwind_speed_model():
    @u.quantity_input
    def __init__(self, r: u.km, model, T: u.Kelvin):
        self.model = model
        self.T = T
        self.r = r
        self.vsw = self.solarwind_speed(r, model, T)

    def solarwind_speed(self, r, model, T):    
        r"""
        Return the solar wind speed from parker solar wind
        solution using the lambertw function. Only this model is
        considered for the moment.
        """

        Rm = const.k_B / const.m_p
        mu = 0.6
        G = const.G
        M_sun = const.M_sun

        if model == 'Parker':
            # From here: https://arxiv.org/pdf/astro-ph/0406176.pdf
            v_c = (Rm*T/mu)**(1/2)
            r_c = G*M_sun/(2*v_c**2)
            D = (r/r_c)**-4 * np.exp(4*(1-r_c/r)-1)
            u_sw = np.zeros(r.shape) * u.m/u.s
            W_0 = lambertw(-(D.decompose()).value, 0)
            u_sw[r < r_c] = np.sqrt(-v_c**2 * W_0[r < r_c])
            W_m1 = lambertw(-(D.decompose()).value, -1)
            u_sw[r >= r_c] = np.sqrt(-v_c**2 * W_m1[r >= r_c])
        else:
            raise ValueError("Non supported solar wind model.")

        u_sw = u_sw.to(u.km/u.s)

        return u_sw           
