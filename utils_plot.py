import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import LineCollection
import streamlit as st

import HDPm
from coronal_models import *

# ---- Plot functions ----
@u.quantity_input
def plot_propagation(smodel, plt_type='disturbance', app = False, extra = None,
                     xlim=(0, 12), ylim=(-6, 6), cmmode='Time'):
    r"""
    This function plots the propagation of the disturbance
    """
    @u.quantity_input
    def overplot_(smodel, extra):
        phi = extra['phi']
        mode = extra['mode']

        if mode == 0: # Connection to point
            r = extra['r']
            cpoints = HDPm.connection_points(phi, r)
            tpp_ = smodel.connection_time_to_point(cpoints)
        elif mode == 1: # Connection at Surface
            r = extra['r'][0]
            cpoints = HDPm.connection_points(phi, r)
            tpp_ = smodel.connection_time_to_point(cpoints)
        elif mode == 2: # First connection to fieldline
            tpp_, r = smodel.first_connection_times(phi)
            cpoints = HDPm.connection_points(phi, r)

        th_BNpp, omega, time = smodel.calculate_theta_BN(cpoints, time=tpp_)
        Vshn = smodel.calculate_vshn(cpoints, time = time, theta_BN = th_BNpp)

        x, y = smodel.propagate(tpp_)
        plt.plot((x[0,:]).to_value(u.R_sun), (y[0,:]).to_value(u.R_sun), color=(0, 0, 0), linewidth=1.0)

        x = r.to_value(u.R_sun) * np.cos(phi.to_value(u.rad))
        y = r.to_value(u.R_sun) * np.sin(phi.to_value(u.rad))
        dsh = (smodel.dsh(tpp_)).to_value(u.R_sun)
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

        plt.text(x+2.7*un,y+2.7*wn,
                 '$\Theta_{Bn} = %2.1f^\\circ$ \n $V_{Shn}=%.0f$ km/s' % (
                 th_BNpp.to_value(u.deg), Vshn.to_value(u.km/u.second)),
                 color=(1, 0, 0), fontsize=11, verticalalignment='center',
                 clip_on=True)
    
    time = np.linspace(0, 120, 25)*(60*u.second)
    x, y =  smodel.propagate(time)

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
            phi, _, r = cart2sph((x).to_value(u.km), (y).to_value(u.km), 0*(x).to_value(u.km), wrap=True)
            cpoints = HDPm.connection_points(phi*u.rad, r*u.km)
            cparams = HDPm.coronal_parameters(cpoints, smodel.coronal_parameters.options)
            MA, ThBN, Vshn, _ = smodel.shock_parameters(cpoints, cparams)
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
            cpoints = HDPm.connection_points(phi, r)
            cparams = HDPm.coronal_parameters(cpoints, smodel.coronal_parameters.options)
            MA, ThBN, Vshn, _ = smodel.shock_parameters(cpoints, cparams)

            x, y, _ = sph2cart(phi, 0*phi, r)
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
    title_text = 'Model Param.: $V_0=%1.0f\\,%s, a_0=%0.2f\\,%s, \\alpha=%1.2f, \\epsilon=%1.2f$' % (smodel.V0.value, smodel.V0.unit, smodel.a0.value,  smodel.a0.unit, smodel.alpha, smodel.epsilon)
    plt.title(title_text, fontsize=10)

    if extra is not None:
        overplot_(smodel, extra)

    plt.tight_layout()
    if app is True:
        st.pyplot(fig)
    else:
        plt.show()

    return plt

@u.quantity_input
def plot_phivstime(smodel, plot_MA = True, app = False, extra = None):
    r"""
    This function plots the phi versus time
    """

    @u.quantity_input
    def overplot_(extra):
        phi = extra['phi']
        mode = extra['mode']

        if mode == 0: # Connection to point
            r = extra['r']
            cpoints = HDPm.connection_points(phi, r)
            tpp_ = smodel.connection_time_to_point(cpoints)
        elif mode == 1: # Connection at Surface
            r = extra['r'][0]
            cpoints = HDPm.connection_points(phi, r)
            tpp_ = smodel.connection_time_to_point(cpoints)
        elif mode == 2: # First connection to fieldline
            tpp_, r = smodel.first_connection_times(phi)
            cpoints = HDPm.connection_points(phi, r)

        cparams = HDPm.coronal_parameters(cpoints, smodel.coronal_parameters.options)
        MA, ThBN, Vshn, _ = smodel.shock_parameters(cpoints, cparams)

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

        plt.text(60,2, 'Parameters: \n $T_{con.} = %2.1f$ min. \n $\Phi = %2.1f^\\circ$ \n $\Theta_{Bn} =     %2.1f^\\circ$ \n' % (
                       tpp_.to_value(u.min),
                       phi.to_value(u.deg),
                       ThBN.to_value(u.deg)),
                       color=(0, 0, 0), fontsize=11, verticalalignment='center')

    phi = smodel.connection_points.lon
    time = smodel.parameters.time_con
    theta_BN = smodel.parameters.Theta_BN
    MA = smodel.parameters.MA

    fig = plt.figure(dpi = 150)
    ax = fig.add_subplot(111)
    levels = np.linspace(0, 90, 91)
    cmap = plt.get_cmap("turbo", len(levels))
    norm = colors.BoundaryNorm(levels, len(levels))
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    plt.contour(phi.to_value(u.deg),
                time.to_value(u.minute),
                theta_BN.to_value(u.deg),
                levels=levels, cmap=cmap, linewidths=0.6)
    plt.contour(phi.to_value(u.deg),
                time.to_value(u.minute),
                theta_BN.to_value(u.deg),
                levels=[45], colors='black',
                linestyles='-.', linewidths=1.0)
    if plot_MA:
        cs = plt.contour(phi.to_value(u.deg),
                         time.to_value(u.minute),
                         MA,
                         levels=[0.5, 1, 1.5, 2, 2.7, 4, 6, 8, 10],
                         colors='black', linestyles='-', linewidths=0.7)
        ax.clabel(cs, cs.levels, inline=True, fmt=r'%1.1f', fontsize=9)

    phi_ = np.linspace(0, np.pi/2, 4*91)*u.rad
    time_fcs, _ = smodel.first_connection_times(phi_, surface = True)
    time_fc, height_fc = smodel.first_connection_times(phi_)
    plt.plot(phi_.to_value(u.deg), time_fcs.to_value(u.minute),
             color='black', linestyle='--', linewidth=1.0)
    plt.plot(phi_.to_value(u.deg), time_fc.to_value(u.minute),
             color='blue', linestyle='-', linewidth=1.0)

    plt.yscale('log')
    plt.ylabel('Connection Time [minutes]')
    plt.xlabel('Longitudinal Separation Angle [degrees]')
    plt.xlim(0, 90)
    plt.ylim(1, 1000)
    cbar = plt.colorbar(sm,ticks=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    # cbar.ax.set_yticklabels(time.to_value(u.min))
    cbar.set_label('$\\Theta_{BN}$ [degrees]')
    title_text = 'Model Param.: $V_0=%1.0f\\,%s, a_0=%0.2f\\,%s, \\alpha=%1.2f, \\epsilon=%1.2f$' % (smodel.V0.value, smodel.V0.unit, smodel.a0.value,  smodel.a0.unit, smodel.alpha, smodel.epsilon)
    plt.title(title_text, fontsize=10)

    if extra is not None:
        overplot_(extra)

    plt.tight_layout()
    if app is True:
        st.pyplot(fig)
    else:
        plt.show()
    return plt

@u.quantity_input
def plot_parameters_profile(smodel, parameter_mode = 'MA', xmode = 'height' , app = False, extra = None):
    @u.quantity_input
    def overplot_(smodel, extra, xmode = 'height'):
        phi = extra['phi']
        mode = extra['mode']

        if mode == 0: # Connection to point
            r = extra['r']
            cpoints = HDPm.connection_points(phi, r)
            tpp_ = smodel.connection_time_to_point(cpoints)
        elif mode == 1: # Connection at Surface
            r = extra['r'][0]
            cpoints = HDPm.connection_points(phi, r)
            tpp_ = smodel.connection_time_to_point(cpoints)
        elif mode == 2: # First connection to fieldline
            tpp_, r = smodel.first_connection_times(phi)
            cpoints = HDPm.connection_points(phi, r)

        if xmode == 'height':
            plt.plot([r.to_value(u.R_sun), r.to_value(u.R_sun)],
                      plt.gca().get_ylim(), color=(0, 0, 0), linestyle = '--',
                      linewidth=1, zorder=105)
        elif xmode == 'time':
            plt.plot([tpp_.to_value(u.minute), tpp_.to_value(u.minute)],
                      plt.gca().get_ylim(), color=(0, 0, 0), linestyle = '--',
                      linewidth=1, zorder=105)

    phi_ = smodel.connection_points.lon[0, :]
    r = smodel.connection_points.distance

    sparameters = smodel.parameters
    time = sparameters.time_con
    theta_BN = sparameters.Theta_BN
    MA = sparameters.MA

    if xmode == 'height':
        xp = r.to_value(u.R_sun)
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
    norm = colors.BoundaryNorm(phi_.to_value(u.degree), siz)
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
    title_text = 'Model Param.: $V_0=%1.0f\\,%s, a_0=%0.2f\\,%s, \\alpha=%1.2f, \\epsilon=%1.2f$' % (smodel.V0.value, smodel.V0.unit, smodel.a0.value,  smodel.a0.unit, smodel.alpha, smodel.epsilon)
    plt.title(title_text, fontsize=10)

    if extra is not None:
        st.write()
        cpoints = HDPm.connection_points(extra['phi'], smodel.connection_points.distance[:, 0])
        cparams = HDPm.coronal_parameters(cpoints, smodel.coronal_parameters.options)
        MA, theta_BN, Vshn, time = smodel.shock_parameters(cpoints, cparams)
        if xmode == 'height':
            xp = (smodel.connection_points.distance).to_value(u.R_sun)
        elif xmode == 'time':
            xp = time.to_value(u.minute)
        if parameter_mode == 'MA':
            parameter = MA.to_value()
        elif parameter_mode == 'ThBN':
            parameter = theta_BN.to_value(u.degree)
        plt.plot(xp, parameter, color=(0,0,0,1), linewidth=1)
        overplot_(smodel, extra, xmode)

    plt.tight_layout()

    if app is True:
        st.pyplot(fig)
    else:
        plt.show()

    return plt

def plot_coronal_models(smodel, app = False):

    r = smodel.connection_points.distance[:, 0]

    fig = plt.figure(figsize=(8,7.5),dpi=100)
    ax1 = fig.add_subplot(221)
    models_ = ('Newkirk','Baumbach','Saito','Leblanc')
    
    for model_ in models_:
        model = density_model(r, model=model_, nfold=1)
        plt.plot(r.to_value(u.R_sun), model.Ne, color=[0, 0.45, 0.74],linewidth=1)
    model = density_model(r, model=smodel.coronal_parameters.options['density_model']['model'],
                          nfold=smodel.coronal_parameters.options['density_model']['nfold'])
    h1 = plt.plot(r.to_value(u.R_sun), model.Ne, color=[0.64, 0.08, 0.18], linewidth=2)
    plt.xlim(1, 215)
    plt.xlabel('Distance from solar center [Rsun]')
    plt.ylabel('Electron Density $[cm^{-1}$]')
    plt.yscale('log')
    ax1.legend(h1, ["{:.2f}".format(smodel.coronal_parameters.options['density_model']['nfold']) + 'x ' +
                    smodel.coronal_parameters.options['density_model']['model']])

    ax2 = fig.add_subplot(222)
    models_ = ('1/r^1','1/r^2','1/r^3')
    for model_ in models_:
        model = magnetic_model(r, model=model_, B0=smodel.coronal_parameters.options['magnetic_model']['B0'])
        plt.plot(r.to_value(u.R_sun), model.B, color=[0, 0.45, 0.74], linewidth=1)
    model = magnetic_model(r, B0=smodel.coronal_parameters.options['magnetic_model']['B0'],
                       model=smodel.coronal_parameters.options['magnetic_model']['model'])
    h2 = plt.plot(r.to_value(u.R_sun), model.B, color=[0.64, 0.08, 0.18], linewidth=2)
    plt.xlim(1, 215)
    plt.xlabel('Distance from solar center [Rsun]')
    plt.ylabel('Total Magnetic Field [gauss]')
    plt.yscale('log')
    ax2.legend(h2, [smodel.coronal_parameters.options['magnetic_model']['model'] +
                    ' ($B_0$=' + "{:.2f}".format(smodel.coronal_parameters.options['magnetic_model']['B0']) + ')'])

    ax3 = fig.add_subplot(223)
    models_ = ('Parker',)
    for model_ in models_:
        model = solarwind_speed_model(r, model=model_, T=smodel.coronal_parameters.options['sw_model']['T'])
        plt.plot(r.to_value(u.R_sun), model.vsw, color=[0, 0.45, 0.74], linewidth=1)
    model = solarwind_speed_model(r, T=smodel.coronal_parameters.options['sw_model']['T'],
                                model=smodel.coronal_parameters.options['sw_model']['model'])
    h3 = plt.plot(r.to_value(u.R_sun), model.vsw, color=[0.64, 0.08, 0.18], linewidth=2)
    plt.xlim(1, 215)
    plt.xlabel('Distance from solar center [Rsun]')
    plt.ylabel('Solar Wind Speed [km/s]')
    ax3.legend(h3, [smodel.coronal_parameters.options['sw_model']['model'] +
                    ' ($T_0$=' + "{:.2f}".format((smodel.coronal_parameters.options['sw_model']['T']).to(u.megaKelvin)) + ')'])

    plt.tight_layout()
    if app is True:
        st.pyplot(fig)
    else:
        plt.show()
    
    return plt

def sph2cart(azimuth, elevation, r, from_degrees=False):
    #convert coords from spheric to cartesian coordinates
    if from_degrees:
    	elevation = elevation*np.pi/180.
    	azimuth = azimuth*np.pi/180.
    x = r*np.cos(elevation)*np.cos(azimuth)
    y = r*np.cos(elevation)*np.sin(azimuth)
    z = r*np.sin(elevation)
    return x,y,z

def cart2sph(x,y,z, to_degrees=False, wrap=False):
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
