"""
Heliospheric Disturbance Propagation Model & Tool (HDPmt)
An open-source software package that built to model the propagation of
heliospheric disturbances.

Copyright (C) 2021  Athanasios Kouloumvakos

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import warnings
import numpy as np
import astropy.units as u
from astropy import constants as const
from astropy.coordinates import BaseCoordinateFrame, SphericalRepresentation, SphericalDifferential
import streamlit as st

import coronal_models

__all__ = ['disturbance_model']


class connection_points(BaseCoordinateFrame):
    default_representation = SphericalRepresentation
    default_differential = SphericalDifferential
    _wrap_angle = 360 * u.deg  # for longitude in spherical coordinates

    def __init__(self, lon: (u.degree), distance: (u.R_sun)):
        super().__init__(lon, 0*u.rad, distance)

class coronal_parameters():
    def __init__(self, connection_points, options):
        
        r = connection_points.distance
        
        self.options = options

        self.density = coronal_models.density_model(r, model=options['density_model']['model'],
                                                    nfold=options['density_model']['nfold'])
        self.magnetic_field = coronal_models.magnetic_model(r, model=options['magnetic_model']['model'],
                                                            topo=options['magnetic_model']['topo'],
                                                            B0=options['magnetic_model']['B0'],
                                                            usw=options['magnetic_model']['usw'])
        self.sw_speed = coronal_models.solarwind_speed_model(r, model=options['sw_model']['model'],
                                                             T=options['sw_model']['T'])

class disturbance:
    r"""
     This is a class for the disturbance
    """
    
    @u.quantity_input
    def __init__(self,
                 a0: (u.km/u.second**2) = 0*(u.km/u.second**2),
                 V0: (u.km/u.second) = 800*(u.km/u.second),
                 alpha = 0.5,
                 epsilon = 1):

        self.a0 = a0  # Disturbance expansion acceleration [km/s^2]
        self.V0 = V0  # Disturbance expansion speed [km/s^2]
        self.alpha = alpha
        self.epsilon = epsilon
        self.connection_points = []
        self.coronal_parameters = []
        self.parameters = self.parameters()

    class parameters:
        time_con = []
        MA = []
        Theta_BN = []
        Vshn = []

    class constants:
        R_sun = const.R_sun
        k_B = const.k_B
        m_p = const.m_p
        G = const.G
        M_sun = const.M_sun

    def set_parameters(self, connection_points, coronal_parameters):
        self.connection_points = connection_points
        self.coronal_parameters = coronal_parameters
        self.parameters.time_con = self.connection_time_to_point(connection_points)
        self.parameters.MA, \
        self.parameters.Theta_BN, \
        self.parameters.Vshn, \
        self.parameters.time_con = self.shock_parameters(connection_points, coronal_parameters)

        return self

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
    def dsh(self, time: u.second):        
        R_sun = self.constants.R_sun
        Rsh = 0.5 * self.a0 * time**2 + self.V0 * time
        dsh = 1 * u.R_sun + self.alpha * Rsh
        return dsh

    @u.quantity_input
    def connection_time_to_point(self, connection_points):
        r"""
        Calculate the connection time of the disturbance to a point
        with polar coordinates (phi[rad],r[km])
        """

        phi = connection_points.lon
        r = connection_points.distance
        
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
    def calculate_theta_BN(self, connection_points, br_angle: u.rad, time: u.second = None):
        r"""
        Calculate the shock geometry (Theta_BN angle) at a point
        with polar coordinates (phi[rad], r[km])
        """

        phi = connection_points.lon
        r = connection_points.distance

        if time is None:
            time = self.connection_time_to_point(connection_points)
        
        if time is None:
            time = self.connection_time_to_point(connection_points)

        R_sun = self.constants.R_sun
        Rsh = 0.5*self.a0*time**2 + self.V0*time
        dsh = R_sun + self.alpha * Rsh
        omega = np.arctan2(1, 1/(r*np.sin(phi)/(r*np.cos(phi)-dsh)))
        omega = np.arctan2(1, 1/(self.epsilon**-2 * np.tan(omega)))
        theta_NR = omega - phi
        theta_NR[theta_NR<0*u.rad] = -(theta_NR[theta_NR<0*u.rad] + (np.pi)*u.rad)
        omega[phi>(np.pi)*u.rad] = phi[phi>(np.pi)*u.rad] - theta_NR[phi>(np.pi)*u.rad]

        theta_BN = (theta_NR + br_angle) % (np.pi*u.rad)

        return theta_BN.decompose(), omega.decompose(), time.decompose()

    @u.quantity_input
    def calculate_vshn(self, connection_points, time: u.second = None, omega: u.rad = None):

        phi = connection_points.lon
        r = connection_points.distance

        if time is None:
            time = self.connection_time_to_point(connection_points)

        if omega is None:
            #  I put here br_angle = 0 because it doesn't mater on the calculation. This will
            #  change in the future.
            _, omega, _ = self.calculate_theta_BN(connection_points, 0*u.rad, time = time) 

        V_ = self.a0 * time + self.V0
        Vshn = (V_ * np.sqrt(np.cos(phi)**2 +
                (self.epsilon)**2 * np.sin(phi)**2) +
                V_ * self.alpha * np.cos(omega))        

        Vshn[Vshn < 0] = np.nan

        return Vshn.decompose()

    @u.quantity_input
    def shock_parameters(self, connection_points, coronal_parameters,
                         time: u.second = None, theta_BN: u.rad = None,
                         omega: u.rad = None):
        r"""
        Calculate the shock parameters (MA, Theta_BN, time) at a point
        with polar coordinates (phi[rad],r[m])
        """

        phi = connection_points.lon
        r = connection_points.distance

        R_sun = self.constants.R_sun
        m_i = self.constants.m_p

        if time is None:
            time = self.connection_time_to_point(connection_points)

        if (theta_BN is None) or (omega is None):
            theta_BN, omega, _ = self.calculate_theta_BN(connection_points, 
                                                         coronal_parameters.magnetic_field.br_angle,
                                                         time = time)

        Ne = coronal_parameters.density.Ne
        B = coronal_parameters.magnetic_field.B
        Vsw = coronal_parameters.sw_speed.vsw

        VA = B / np.sqrt(const.mu0*m_i*(1.92*Ne))

        Vshn = self.calculate_vshn(connection_points, time = time, omega = omega);

        MA = np.abs(Vshn-Vsw*np.cos(theta_BN))/VA

        return MA.decompose(), theta_BN.decompose(), Vshn.decompose(), time.decompose()

    @u.quantity_input
    def first_connection_times(self, phi: u.rad, surface=False):
        r"""
        Calculate the time and height that the disturbance connects
        to a field line for the first time. Phi[rad] is the separation angle of the
        field line with respect to the origin of the disturbance
          #time_fcs, _ = self.first_connection_times(phi, surface = True)
          #time_fc, height_fc = self.first_connection_times(phi)
        """
        
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

        time = self.connection_time_to_point(connection_points(phi, r))

        return time.decompose(), r.decompose()
