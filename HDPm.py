import warnings
import numpy as np
import astropy.units as u
from astropy import constants as const
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import LineCollection
from scipy.special import lambertw
import streamlit as st

__all__ = ['disturbance_model']


class disturbance_model:
    r"""
    Examples
    --------

    Example calculation and ploting using default values

    >>> import astropy.units as u
    >>> from coroshock.HDPmt import HDPm
    >>> smodel = HDPm.disturbance_model(0*(u.km/u.second**2),800*(u.km/u.second))
    >>> x,y = smodel.propagate(smodel.products.time_propagation)
    >>> smodel.plot_propagation()
    >>> smodel.connection_time_to_point(45*u.deg,2*u.R_sun)
    >>> smodel.first_connection_times(45*u.deg)
    >>> smodel.calculate_theta_BN(45*u.deg,2*u.R_sun)
    >>> smodel.density_models(2*u.R_sun,model='Saito')
    >>> smodel.plot_phivstime()
    """
    @u.quantity_input
    def __init__(self, 
                 a0: (u.km/u.second**2) = 0*(u.km/u.second**2),
                 V0: (u.km/u.second) = 800*(u.km/u.second),
                 alpha = 0.5,
                 epsilon = 1,
                 density_model = 'Saito',
                 NFold = 1,
                 magnetic_model = '1/r^2',
                 B0: (u.gauss) = 2.2 * u.gauss,
                 sw_model = 'Parker',
                 T0: (u.Kelvin) = 1.4 * 10**6 * u.Kelvin):
        self.a0 = a0  # Disturbance expansion acceleration [km/s^2]
        self.V0 = V0  # Disturbance expansion speed [km/s^2]
        self.alpha = alpha
        self.epsilon = epsilon
        self.coronal_models = disturbance_model.coronal_models()
        self.coronal_models.density_model.type = density_model
        self.coronal_models.density_model.NFold = NFold
        self.coronal_models.magnetic_model.type = magnetic_model
        self.coronal_models.magnetic_model.B0 = B0
        self.coronal_models.sw_model.type = sw_model
        self.coronal_models.sw_model.T0 = T0
        # self.products = disturbance_model.products()

    class connection_points:
        phi = np.linspace(0, np.pi/2, 4*91)*u.rad
        r = np.logspace(np.log10(1*u.R_sun.in_units(u.km)),
                        np.log10(215*u.R_sun.in_units(u.km)), 250)*u.km
        phi_mesh, r_mesh = np.meshgrid(phi, r)

    class coronal_models:
        class density_model:
            type = []
            NFold = []

        class magnetic_model:
            type = []
            B0 = []

        class sw_model:
            type = []
            T0 = []

    class _make_disturbance:
        time_propagation = np.linspace(0, 120, 25)*(60*u.second)
        class coordinates:
            x = []
            y = []
        coords = coordinates()

    class _make_parameters:
        time = []
        time_fcs = []
        time_fc = []
        MA = []
        Theta_BN = []
        Vshn = []

    class constants:
        R_sun = const.R_sun
        k_B = const.k_B
        m_p = const.m_p
        G = const.G
        M_sun = const.M_sun

    @property
    def parameters(self):
        parameters = disturbance_model._make_parameters()
        parameters.time_fcs, _ = self.first_connection_times(surface = True)
        parameters.time_fc, _  = self.first_connection_times()
        parameters.MA, parameters.Theta_BN, parameters.Vshn, parameters.time = self.shock_parameters(do_mesh=True);
        return parameters

    @parameters.setter
    def parameters(self, parameters):
        self._parameters = parameters

    @property
    def disturbance(self):
        disturbance = disturbance_model._make_disturbance()
        disturbance.coords.x, disturbance.coords.y = self.propagate(disturbance.time_propagation)
        return disturbance
        
    @disturbance.setter
    def disturbance(self, disturbance):
        self._disturbance = disturbance

    @u.quantity_input
    def dsh(self, time: u.second):        
        R_sun = self.constants.R_sun
        Rsh = 0.5*self.a0*time**2 + self.V0*time
        dsh = R_sun + self.alpha * Rsh
        return dsh

    @u.quantity_input
    def propagate(self, time: u.second):
        # Calculate the x,y[km] coordinates of the disturbance at
        # a given time[sec]. Input time can be a vecor also, in this
        # case the x,y coordinates will be matrices, time per row.

        R_sun = self.constants.R_sun

        omega = np.linspace(0, 2*np.pi, 361)*u.rad
        x = np.zeros((time.size, omega.size))*u.km
        y = np.zeros((time.size, omega.size))*u.km
        
        time = np.atleast_1d(time)
        
        for i in range(time.size):
            Rsh_x = 0.5 * self.a0 * time[i]**2 + self.V0 * time[i]
            Rsh_y = self.epsilon * Rsh_x
            Rrc = self.alpha * Rsh_x
            x[i, :] = np.cos(omega) * Rsh_x + (R_sun+Rrc)
            y[i, :] = np.sin(omega) * Rsh_y

        return x, y

    @u.quantity_input
    def connection_time_to_point(self, phi: u.rad, r: u.km):
        r"""
        Calculate the connection time of the disturbance to a point
        with polar coordinates (phi[rad],r[km])
        """

        R_sun = self.constants.R_sun

        # Since the point-(r,phi) are known, I use in this calculation
        # the ellipse equation in polar cordinates and I solve for Rsh.
        # The initial eq. is e^2 (r cos(phi)-dsh)^2 + (r sin(phi))^2 = e^2 Rsh^2
        a_ = self.epsilon**2 * (self.alpha**2-1)
        b_ = 2*self.epsilon**2 * self.alpha * (R_sun - r*np.cos(phi))
        c_ = r**2 * (self.epsilon**2 * np.cos(phi)**2 + np.sin(phi)**2) + self.epsilon**2 * (R_sun**2 - 2 * R_sun * r * np.cos(phi))
        rc_pp = (-b_ - np.sqrt(b_**2 - 4*a_*c_))/(2*a_)

        if self.a0 == 0:
            time = rc_pp/self.V0
        else:
            time = (-self.V0 + np.sqrt(self.V0**2 + 2 * self.a0 * rc_pp))/self.a0

        time[time.imag.value != 0] = np.nan

        return time.decompose()

    @u.quantity_input
    def first_connection_times(self, phi: u.rad = None, surface=False):
        r"""
        Calculate the time and height that the disturbance connects
        to a field line for the first time. Phi[rad] is the separation angle of the
        field line with respect to the origin of the disturbance
        """
        
        if phi is None:
            phi = self.connection_points.phi

        R_sun = self.constants.R_sun

        if surface:
            r = R_sun
        else:
            # The following calculation it is based on finding the
            # tangent line to an elipse. Start from the ellipse
            # equation in cartesian coordinates, e^2 (x-dsh)^2 + (y)^2
            # = e^2 Rsh^2 and the equation of line y = mx.
            # Solve for x, and calculate Rsh by seting the discriminant
            # equal to zero. Find x from -b/2a and then r from r =
            # sqrt(x^2+y^2) seting y=mx and x=-b/2a.
            # (a^2 - k/tan(phi)^2) Rsh^2 + 2 a R_sun Rsh + R_sun^2 = 0
            k = np.tan(phi)**2 + self.epsilon**2
            l = (self.alpha**2 - k/np.tan(phi)**2)
            Rsh = R_sun * (- np.sqrt(self.alpha**2 - l) - self.alpha) / l
            dsh = R_sun + self.alpha * Rsh
            x = 2 * self.epsilon**2 * dsh / (2*k)
            r = np.abs(x) * np.sqrt(1+np.tan(phi)**2)

            r[r.imag.value != 0] = R_sun  # (imag(r)~=0|isnan(r))
            r[r < R_sun] = R_sun

        time = self.connection_time_to_point(phi, r)

        return time.decompose(), r.decompose()

    @u.quantity_input
    def calculate_theta_BN(self, phi: u.rad, r: u.km, time: u.second = None):
        r"""
        Calculate the shock geometry (Theta_BN angle) at a point
        with polar coordinates (phi[rad],r[km])
        """

        if time is None:
            time = self.connection_time_to_point(phi, r)

        R_sun = self.constants.R_sun
        Rsh = 0.5*self.a0*time**2 + self.V0*time
        dsh = R_sun + self.alpha * Rsh
        omega = np.arctan2(1, 1/(r*np.sin(phi)/(r*np.cos(phi)-dsh)))
        omega = np.arctan2(1, 1/(self.epsilon**-2 * np.tan(omega)))
        theta_BN = omega - phi
        theta_BN[theta_BN>(np.pi/2)*u.rad] = -(theta_BN[theta_BN>(np.pi/2)*u.rad] - (np.pi)*u.rad)
        
        return theta_BN.decompose(), omega.decompose(), time.decompose()

    @u.quantity_input
    def calculate_vshn(self, phi: u.rad = None, r: u.km = None, 
                       time: u.second = None, theta_BN: u.rad = None):

        if phi is None:
            phi = self.connection_points.phi_mesh
        if r is None:
            r = self.connection_points.r_mesh
        if time is None:
            time = self.connection_time_to_point(phi, r)  
        if theta_BN is None:
            theta_BN = self.calculate_theta_BN(phi, r, time = time)

        omega = np.abs(phi) + theta_BN;
        V_ = self.a0 * time + self.V0
        Vshn = (V_ * np.sqrt(np.cos(phi)**2 +
                (self.epsilon)**2 * np.sin(phi)**2) +
                V_ * self.alpha * np.cos(omega))        

        Vshn[Vshn < 0] = np.nan

        return Vshn.decompose()

    @u.quantity_input
    def shock_parameters(self, phi: u.rad = None, r: u.km = None, do_mesh = False):
        r"""
        Calculate the shock parameters (MA, Theta_BN, time) at a point
        with polar coordinates (phi[rad],r[m])
        """

        if phi is None:
            phi = self.connection_points.phi
        if r is None:
            r = self.connection_points.r

        R_sun = self.constants.R_sun
        m_i = self.constants.m_p

        # TODO: Do we need to control here when to mesh grid?
        if do_mesh is True:
            phi_, r_ = np.meshgrid(phi, r)
        else:
            phi_, r_ = phi, r
        
        time = self.connection_time_to_point(phi_, r_)
        theta_BN, _, time = self.calculate_theta_BN(phi_, r_, time=time)
        
        Ne = self.density_models(r_)
        B = self.magnetic_models(r_)
        VA = B / np.sqrt(const.mu0*m_i*(1.92*Ne))
        Vshn = self.calculate_vshn(phi, r, time = time, theta_BN = theta_BN);
        Vsw = self.solar_wind_models(r_)
        MA = np.abs(Vshn-Vsw*np.cos(theta_BN))/VA

        return MA.decompose(), theta_BN.decompose(), Vshn.decompose(), time.decompose()

    # ---- General Functions for Coronal Models ----
    @u.quantity_input
    def density_models(self, r: u.km, model=None, NFold=None):
        r"""
        Return the electron density profile from different empirical coronal models.
        """

        R_sun = self.constants.R_sun

        if model is None:
            model = self.coronal_models.density_model.type
        
        if NFold is None:
            NFold = self.coronal_models.density_model.NFold

        if model == 'Newkirk':
            N0 = 4.2 * 10**4
            Ne = N0 * 10**(4.32*R_sun/r)
        elif model == 'Baumbach':  # Baumbach-Allen Density model (From Aswaden page 82)
            Ne = (2.99*(r/R_sun)**(-16) + 1.55*(r/R_sun)**(-6) + 0.036*(r/R_sun)**(-1.5))*10**8
        elif model == 'Saito':  # Saito Density model
            Ne = (3.09*(r/R_sun)**(-16) + 1.58*(r/R_sun)**(-6) + 0.0251*(r/R_sun)**(-2.5))*10**8
        elif model == 'Leblanc':  # Leblanch Density model
            Ne = (0.8*10**8*(r/R_sun)**(-6) + 0.41*10**7*(r/R_sun)**(-4) + 0.33*10**6*(r/R_sun)**(-2))
        else:
            raise ValueError("Non supported density model.")

        # All the models should give a density in cm^-3 before returning.
        Ne = NFold * ((Ne * u.cm**-3).decompose()).to(u.cm**-3)

        return Ne

    @u.quantity_input
    def magnetic_models(self, r: u.km, B0: u.gauss = None, model=None, alpha=1):
        r"""
        Return the total magnetic field profile from coronal magnetic
        field models. Only quiet sun models are considered with
        recomended photosheric magnetic field of B0 = 2.2 gauss.
        The 1/r^2 model is recomended.
        """

        R_sun = self.constants.R_sun

        if B0 is None:
            B0 = self.coronal_models.magnetic_model.B0

        if model is None:
            model = self.coronal_models.magnetic_model.type

        if model == '1/r^1':
            B = B0 * (R_sun/r)**1
        elif model == '1/r^2':
            B = B0 * (R_sun/r)**2
        elif model == '1/r^3':
            B = B0 * (R_sun/r)**3
        elif model == '1/r^a':
            B = B0 * (R_sun/r)**alpha
        else:
            raise ValueError("Non supported magnetic field model.")

        B = (B.decompose()).to(u.gauss)

        return B

    @u.quantity_input
    def solar_wind_models(self, r: u.km, T: u.Kelvin = None, model=None):
        r"""
        Return the solar wind speed from parker solar wind
        solution using the lambertw function. Only this model is
        considered for the moment.
        """

        k_B = self.constants.k_B
        m_p = self.constants.m_p
        Rm = k_B / m_p
        mu = 0.6
        G = self.constants.G
        M_sun = self.constants.M_sun

        if T is None:
            T = self.coronal_models.sw_model.T0

        if model is None:
            model = self.coronal_models.sw_model.type

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

    # ---- Plot functions ----
    @u.quantity_input
    def plot_propagation(self, plt_type='disturbance', app = False, extra = None,
                         xlim=(0, 12), ylim=(-6, 6), cmmode='Time'):
        r"""
        This function plots the propagation of the disturbance
        """
        @u.quantity_input
        def overplot_(self, phi: u.rad, r: u.km, mode = 0):
            if mode == 0:
                tpp_ = self.connection_time_to_point(phi, r)
            elif mode == 1:
                tpp_, r = self.first_connection_times(phi)

            th_BNpp, omega, time = self.calculate_theta_BN(phi, r, time=tpp_)

            x, y = self.propagate(tpp_)
            plt.plot((x[0,:]).to_value(u.R_sun), (y[0,:]).to_value(u.R_sun), color=(0, 0, 0), linewidth=1.0)

            x = r.to_value(u.R_sun) * np.cos(phi.to_value(u.rad))
            y = r.to_value(u.R_sun) * np.sin(phi.to_value(u.rad))
            dsh = (self.dsh(tpp_)).to_value(u.R_sun)
            plt.plot([0, x], [0, y], color=(0, 0, 0), linestyle = '--', linewidth=0.8, zorder=105);
            plt.plot([dsh, x],[0, y], color=(0, 0, 0), linestyle = '--', linewidth=0.8, zorder=105);

            plt.scatter(0, 0, s=25, marker='o', color=(0, 0, 0), zorder=105)
            plt.scatter(x, y, s=25, marker='o', color=(0, 0, 0), zorder=105)
            plt.scatter(dsh, 0, s=25, marker='o', color=(0, 0, 0), zorder=105)
        
            un = x / np.sqrt(x**2+y**2)
            wn = y / np.sqrt(x**2+y**2)
            plt.quiver(x, y, un, wn, scale=4.5, zorder=105)
            un = (x - dsh) / np.sqrt((x - dsh)**2 + y**2)
            wn = y / np.sqrt((x - dsh)**2 + y**2)
            plt.quiver(x, y, un, wn, scale=4.5, zorder=105)
            if th_BNpp + phi < np.pi*u.rad/2:
                signe_ = 1
            elif th_BNpp + phi > np.pi*u.rad/2:
                signe_ = -1
            else:
                signe_ = 0

            un = signe_ / np.sqrt( 1**2 + np.tan(omega.to_value(u.rad))**2)
            wn = signe_ * np.tan(omega.to_value(u.rad)) / np.sqrt( 1**2 + np.tan(omega.to_value(u.rad))**2)
            plt.quiver(x, y, un, wn, scale = 4.5, color = (1, 0, 0), zorder = 105)
            Vshn = self.calculate_vshn(phi, r, time, th_BNpp)
            plt.text(x+2.7*un,y+2.7*wn,
                              '$\Theta_{Bn} = %2.1f^\\circ$ \n $V_{Shn}=%.0f$ km/s' % (
                              th_BNpp.to_value(u.deg), Vshn.to_value(u.km/u.second)),
                              color=(1, 0, 0), fontsize=11, verticalalignment='center',
                              clip_on=True)

        tdist_ =  self.disturbance
        time = tdist_.time_propagation
        x = tdist_.coords.x
        y = tdist_.coords.y

        fig = plt.figure(dpi=150)
        ax = fig.add_subplot(111, aspect='equal')
        
        if plt_type == 'disturbance':
            if cmmode=='Time':
                cmap = plt.get_cmap("jet_r", len(time))
                norm = colors.BoundaryNorm(np.arange(len(time)+1), len(time))
                sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
                cbar = plt.colorbar(sm, ticks=np.arange(0.5, len(time) ))
                cbar.ax.set_yticklabels(time.to_value(u.min))
                cbar.set_label('Time [minutes]')
                for i in range(time.size):
                    plt.plot((x[i, :]).to_value(u.R_sun),
                             (y[i, :]).to_value(u.R_sun),
                             color=cmap(i), linewidth=1.0)
            else:
                phi, _, r = cart2sph((x).to_value(u.km), (y).to_value(u.km), 0*(x).to_value(u.km))
                MA, ThBN, Vshn, _ = self.shock_parameters(phi * u.rad, r * u.km)
                if cmmode=='Theta_BN':
                    norm = plt.Normalize(0, 90)
                    parameter = (ThBN.reshape(x.shape)).to_value(u.deg)
                    label = '$\Theta_{Bn}$ [degrees]'
                if cmmode=='Mach Alfven':
                    norm = plt.Normalize(1, np.ceil(np.nanquantile(MA, 0.80)))
                    parameter = MA.reshape(x.shape)
                    label = 'Mach Alfven'
                    
                for i in range(ThBN.shape[0]):
                    points = np.array([(x[i,:]).to_value(u.R_sun), 
                                       (y[i,:]).to_value(u.R_sun)
                                       ]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)
                    lc = LineCollection(segments, cmap='jet', norm=norm)
                    lc.set_array(parameter[i,:])
                    lc.set_linewidth(1)
                    line = ax.add_collection(lc)
                cbar = plt.colorbar(line, ax=ax)
                cbar.set_label(label)
                
        elif plt_type == 'fieldlines':
            if cmmode=='Phi':
                phi = np.linspace(-90,90,19)
                top = plt.cm.get_cmap('jet_r', 9)
                bottom = plt.cm.get_cmap('jet', 9)

                newcolors = np.vstack((top(np.linspace(0, 1, 9)),
                                       (0,0,0,1),
                                       bottom(np.linspace(0, 1, 9))))
                cmap = colors.ListedColormap(newcolors)
                norm = colors.BoundaryNorm(np.arange(len(phi)+1), len(phi))
                sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
                for i, phi_ in enumerate(phi):
                    x, y, z = sph2cart(phi_, 0, 20, from_degrees=True)
                    plt.plot((0,x), (0,y),
                             color=cmap(i), linewidth=1.0)
                cbar = plt.colorbar(sm, ticks=np.arange(0.5, len(phi) ))
                cbar.ax.set_yticklabels(phi)
                cbar.set_label('Phi [degrees]')
            else:
                phi = np.linspace(-90, 90, 19) * u.deg
                r = np.linspace(  0, 20, 120) * u.R_sun
                phi, r = np.meshgrid(phi, r)
                MA, ThBN, Vshn, _ = self.shock_parameters(phi, r)
                x, y, _ = sph2cart(phi, 0*phi,r)
                if cmmode=='Theta_BN':
                    norm = plt.Normalize(0, 90)
                    parameter = (ThBN.reshape(x.shape)).to_value(u.deg)
                    label = '$\Theta_{Bn}$ [degrees]'
                if cmmode=='Mach Alfven':
                    norm = plt.Normalize(1, np.ceil(np.nanquantile(MA, 0.80)))
                    parameter = MA.reshape(x.shape)
                    label = 'Mach Alfven'
                for i in range(ThBN.shape[1]):
                    points = np.array([(x[:,i]).to_value(u.R_sun), 
                                       (y[:,i]).to_value(u.R_sun)
                                       ]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)
                    lc = LineCollection(segments, cmap='jet', norm=norm)
                    lc.set_array(parameter[:,i])
                    lc.set_linewidth(1)
                    line = ax.add_collection(lc)
                cbar = plt.colorbar(line, ax=ax)
                cbar.set_label(label)
      
        circle = plt.Circle((0, 0), facecolor=[1.00, 0.41, 0.16],
                            edgecolor='black', radius=1, zorder=100)
        plt.gca().add_patch(circle)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xlabel('x-axis [Rsun]')
        plt.ylabel('y-axis [Rsun]')
        plt.grid(b=True)
        title_text = 'Model Param.: $V_0=%1.0f\\,%s, a_0=%0.2f\\,%s, \\alpha=%1.2f, \\epsilon=%1.2f$' % (self.V0.value, self.V0.unit, self.a0.value,  self.a0.unit, self.alpha, self.epsilon)
        plt.title(title_text, fontsize=10)
        
        if extra is not None:
            overplot_(self, extra['phi'], extra['r'], extra['mode'])
        
        plt.tight_layout()
        if app is True:
            st.pyplot(fig)
        else:
            plt.show()

        return plt

    @u.quantity_input
    def plot_phivstime(self, plot_MA = True, app = False, extra = None):
        r"""
        This function plots the phi versus time
        """
        
        @u.quantity_input
        def overplot_(self, phi: u.rad, r: u.km, mode = 0):
            if mode == 0:
                tpp_ = self.connection_time_to_point(phi,r)
            elif mode == 1:
                tpp_, r = self.first_connection_times(phi)
            
            plt.plot([phi.to_value(u.deg), phi.to_value(u.deg)],
                     [1, 1e3], color=(0, 0, 0), linestyle = '--',
                               linewidth=1, zorder=105)
            plt.plot([0, 90],
                     [tpp_.to_value(u.min), tpp_.to_value(u.min)],
                     color=(0, 0, 0), linestyle = '--',
                     linewidth=1, zorder=105)
            plt.scatter(phi.to_value(u.deg),
                        tpp_.to_value(u.min),
                        s=25, marker='o', color=(0, 0, 0), zorder=105)
        
            MA, theta_BN, Vshn, time = self.shock_parameters(phi, r, do_mesh=True)
            plt.text(60,2, 'Parameters: \n $T_{con.} = %2.1f$ min. \n $\Phi = %2.1f^\\circ$ \n $\Theta_{Bn} =     %2.1f^\\circ$ \n' % (
                           tpp_.to_value(u.min), phi.to_value(u.deg), theta_BN.to_value(u.deg)),
                           color=(0, 0, 0), fontsize=11, verticalalignment='center')
        
        phi = self.connection_points.phi
        phi_mesh = self.connection_points.phi_mesh

        sparameters = self.parameters
        time_cs = sparameters.time_fcs
        time_cf = sparameters.time_fc
        time = sparameters.time
        theta_BN = sparameters.Theta_BN
        MA = sparameters.MA

        fig = plt.figure(dpi = 150)
        ax = fig.add_subplot(111)
        levels = np.linspace(0, 90, 91)
        cmap = plt.get_cmap("turbo", len(levels))
        norm = colors.BoundaryNorm(levels, len(levels))
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.contour(phi_mesh.to_value(u.deg),
                    time.to_value(u.minute),
                    theta_BN.to_value(u.deg),
                    levels=levels, cmap=cmap, linewidths=0.6)
        plt.contour(phi_mesh.to_value(u.deg),
                    time.to_value(u.minute),
                    theta_BN.to_value(u.deg),
                    levels=[45], colors='black',
                    linestyles='-.', linewidths=1.0)
        plt.plot(phi.to_value(u.deg), time_cs.to_value(u.minute),
                 color='black', linestyle='--', linewidth=1.0)
        plt.plot(phi.to_value(u.deg), time_cf.to_value(u.minute),
                 color='blue', linestyle='-', linewidth=1.0)
        if plot_MA:
            cs = plt.contour(phi_mesh.to_value(u.deg),
                             time.to_value(u.minute),
                             MA,
                             levels=[0.5, 1, 1.5, 2, 2.7, 4, 6, 8, 10],
                             colors='black', linestyles='-', linewidths=0.7)
            ax.clabel(cs, cs.levels, inline=True, fmt=r'%1.1f', fontsize=9)

        plt.yscale('log')
        plt.ylabel('Connection Time [minutes]')
        plt.xlabel('Longitudinal Separation Angle [degrees]')
        plt.xlim(0, 90)
        plt.ylim(1, 1000)
        cbar = plt.colorbar(sm,ticks=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
        # cbar.ax.set_yticklabels(time.to_value(u.min))
        cbar.set_label('$\\Theta_{BN}$ [degrees]')
        title_text = 'Model Param.: $V_0=%1.0f\\,%s, a_0=%0.2f\\,%s, \\alpha=%1.2f, \\epsilon=%1.2f$' % (self.V0.value, self.V0.unit, self.a0.value,  self.a0.unit, self.alpha, self.epsilon)
        plt.title(title_text, fontsize=10)
        
        if extra is not None:
            overplot_(self, extra['phi'], extra['r'], extra['mode'])
        
        plt.tight_layout()
        if app is True:
            st.pyplot(fig)
        else:
            plt.show()

        return plt, sparameters

    @u.quantity_input
    def plot_parameters_profile(self, parameter_mode = 'MA', xmode = 'height' , app = False, extra = None):

        @u.quantity_input
        def overplot_(self, phi: u.rad, r: u.km, mode = 0, xmode = 'height'):
            if mode == 0:
                tpp_ = self.connection_time_to_point(phi,r)
            elif mode == 1:
                tpp_, r = self.first_connection_times(phi)

            if xmode == 'height':
                plt.plot([r.to_value(u.R_sun), r.to_value(u.R_sun)],
                          plt.gca().get_ylim(), color=(0, 0, 0), linestyle = '--',
                          linewidth=1, zorder=105)
            elif xmode == 'time':
                plt.plot([tpp_.to_value(u.minute), tpp_.to_value(u.minute)],
                          plt.gca().get_ylim(), color=(0, 0, 0), linestyle = '--',
                          linewidth=1, zorder=105)

        phi = self.connection_points.phi
        r_mesh = self.connection_points.r_mesh

        sparameters = self.parameters
        time_cs = sparameters.time_fcs
        time_cf = sparameters.time_fc
        time = sparameters.time
        theta_BN = sparameters.Theta_BN
        MA = sparameters.MA

        if xmode == 'height':
            xp = r_mesh.to_value(u.R_sun)
            xlm = (1, 215)
            xlbl = 'Connection Height [Rs]'
        elif xmode == 'time':
            xp = time.to_value(u.minute)
            xlm = (1, 1e3)
            xlbl = 'Connection time [minutes]'

        if parameter_mode == 'MA':
            parameter = MA.to_value()
            yscl = 'log'
            ylbl = 'Mach Alfven'
            ylm = (0.5, 100)
        elif parameter_mode == 'ThBN':
            parameter = theta_BN.to_value(u.degree)
            yscl = 'linear'
            ylbl = '$\Theta_{BN}$ [degrees]'
            ylm = (0, 90)

        siz = (xp.shape)[1]
        fig = plt.figure(dpi=150)
        ax = fig.add_subplot(111)
        cmap = plt.get_cmap("turbo", siz)
        norm = colors.BoundaryNorm(phi.to_value(u.degree), siz)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

        for i in range(0,siz):
            plt.plot(xp[:,i], parameter[:,i], color=cmap(i), linewidth=0.5)

        plt.xscale('log')
        plt.yscale(yscl)
        plt.ylabel(ylbl)
        plt.xlabel(xlbl)
        plt.xlim(xlm)
        plt.ylim(ylm)
        cbar = plt.colorbar(sm,ticks=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
        cbar.set_label('$\\Phi$ [degrees]')
        title_text = 'Model Param.: $V_0=%1.0f\\,%s, a_0=%0.2f\\,%s, \\alpha=%1.2f, \\epsilon=%1.2f$' % (self.V0.value, self.V0.unit, self.a0.value,  self.a0.unit, self.alpha, self.epsilon)
        plt.title(title_text, fontsize=10)

        if extra is not None:
            MA, theta_BN, Vshn, time = self.shock_parameters(extra['phi'], self.connection_points.r, do_mesh=False)
            if xmode == 'height':
                xp = (self.connection_points.r).to_value(u.R_sun)
            elif xmode == 'time':
                xp = time.to_value(u.minute)
            if parameter_mode == 'MA':
                parameter = MA.to_value()
            elif parameter_mode == 'ThBN':
                parameter = theta_BN.to_value(u.degree)
            plt.plot(xp, parameter, color=(0,0,0,1), linewidth=1)
            overplot_(self, extra['phi'], extra['r'], extra['mode'], xmode)

        plt.tight_layout()

        if app is True:
            st.pyplot(fig)
        else:
            plt.show()

        return plt

    def plot_coronal_models(self, app = False):
        
        r = self.connection_points.r

        fig = plt.figure(figsize=(8,7.5),dpi=100)
        ax1 = fig.add_subplot(221)
        models_ = ('Newkirk','Baumbach','Saito','Leblanc')
        for model_ in models_:
            Ne = self.density_models(r,model=model_,NFold=1)
            plt.plot(r.to_value(u.R_sun),Ne,color=[0, 0.45, 0.74],linewidth=1)
        Ne = self.density_models(r,
                                 model=self.coronal_models.density_model.type,
                                 NFold=self.coronal_models.density_model.NFold)
        h1 = plt.plot(r.to_value(u.R_sun),Ne,color=[0.64, 0.08, 0.18],linewidth=2)
        plt.xlim(1, 215)
        plt.xlabel('Distance from solar center [Rsun]')
        plt.ylabel('Electron Density $[cm^{-1}$]')
        plt.yscale('log')
        ax1.legend(h1, [str(self.coronal_models.density_model.NFold) + 'x ' +
                        self.coronal_models.density_model.type])

        ax2 = fig.add_subplot(222)
        models_ = ('1/r^1','1/r^2','1/r^3')
        for model_ in models_:
            B = self.magnetic_models(r, model=model_)
            plt.plot(r.to_value(u.R_sun), B, color=[0, 0.45, 0.74], linewidth=1)
        B = self.magnetic_models(r,model=self.coronal_models.magnetic_model.type)
        h2 = plt.plot(r.to_value(u.R_sun), B, color=[0.64, 0.08, 0.18], linewidth=2)
        plt.xlim(1, 215)
        plt.xlabel('Distance from solar center [Rsun]')
        plt.ylabel('Total Magnetic Field [gauss]')
        plt.yscale('log')
        ax2.legend(h2, [self.coronal_models.magnetic_model.type +
                        ' ($B_0$=' + str(self.coronal_models.magnetic_model.B0) + ')'])
       
        ax3 = fig.add_subplot(223)
        models_ = ('Parker',)
        for model_ in models_:
            usw = self.solar_wind_models(r, model=model_)
            plt.plot(r.to_value(u.R_sun), usw, color=[0, 0.45, 0.74], linewidth=1)
        usw = self.solar_wind_models(r, model=self.coronal_models.sw_model.type)
        h3 = plt.plot(r.to_value(u.R_sun), usw, color=[0.64, 0.08, 0.18], linewidth=2)
        plt.xlim(1, 215)
        plt.xlabel('Distance from solar center [Rsun]')
        plt.ylabel('Solar Wind Speed [km/s]')
        ax3.legend(h3, [self.coronal_models.sw_model.type +
                       ' ($T_0$=' + str((self.coronal_models.sw_model.T0).to(u.megaKelvin)) + ')'])

        plt.tight_layout()
        if app is True:
            st.pyplot(fig)
        else:
            plt.show()
            
        return plt

def sph2cart(azimuth,elevation,r,from_degrees=False):
    #convert coords from spheric to cartesian coordinates
    if from_degrees:
    	elevation = elevation*np.pi/180.
    	azimuth = azimuth*np.pi/180.
    x = r*np.cos(elevation)*np.cos(azimuth)
    y = r*np.cos(elevation)*np.sin(azimuth)
    z = r*np.sin(elevation)
    return x,y,z

def cart2sph(x,y,z,to_degrees=False,wrap=False):
    hypotxy = np.sqrt(x**2+y**2)
    r = np.sqrt(hypotxy**2+z**2)
    elevation = np.arctan2(z,hypotxy)
    azimuth = np.arctan2(y,x)
    if wrap:
        azimuth = azimuth % (2*np.pi)
    
    if to_degrees:
        azimuth = np.rag2deg(azimuth)
        elevation = np.rag2deg(elevation)
    return azimuth, elevation, r
