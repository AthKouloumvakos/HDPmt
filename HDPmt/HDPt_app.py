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


import io
import json
import zipfile

import astropy.units as u
import coronal_models
import HDPm
import numpy as np
import streamlit as st
import utils_app
import utils_plot
from astropy.utils.misc import JsonCustomEncoder
from config import app_styles
from extensions.buttons import download_button
from streamlit.logger import get_logger

LOGGER = get_logger(__name__)


def run():

    # set page config
    st.set_page_config(page_title='HDPt', page_icon=':rocket:',
                       initial_sidebar_state='expanded')

    #############################################################
    # HTML Styles
    app_styles.apply(st)

    #############################################################
    # Main page information text
    st.title('Heliospheric Disturbance Propagation Tool (HDPt)')
    st.markdown("""
                ** ⏻ Adjust the disturbance parameters and select
                from different plot options.**
                """)
    st.markdown('---')

    #############################################################
    # Tools

    # Main Sliders
    st.sidebar.markdown('## Disturbance Parameters')
    V0 = st.sidebar.slider('Expansion Speed [km/s]:', 100, 3000, 800)
    a0 = st.sidebar.slider('Expansion Acceleration [m/s2]:', -100, 100, 0)
    alpha = st.sidebar.slider('alpha:', 0., 1., 0.5,
                              step=0.01, format='%1.2f')
    epsilon = st.sidebar.slider('epsilon:', 0.5, 1.5, 1.,
                                step=0.01, format='%1.2f')

    if V0 < 300 or V0 > 2500 or alpha < 0.15 or alpha > 0.85 or epsilon > 1.3 or epsilon < 0.7:
        st.sidebar.warning('Warning: Parameter value is ambiguous.')

    # Explore Section
    st.sidebar.markdown('**Explore Connections:**')
    mode_connect = st.sidebar.selectbox('Select a mode',
                                        ['-', 'Connection to Point',
                                         'Connection at Surface',
                                         'First Connection'])
    if mode_connect == 'Connection to Point':
        r = st.sidebar.slider('Distance from Sun center [Rsun]:',
                              1., 10., 3., step=0.1, format='%1.1f')
        phi = st.sidebar.slider('Helio-Longitude [deg]:',
                                0., 90., 45., step=0.1, format='%1.1f')
        extra_ = {'phi': phi * u.deg, 'r': r * u.R_sun, 'mode': 0}
    elif mode_connect == 'Connection at Surface':
        r = np.logspace(np.log10(1),
                        np.log10(215), 250)
        phi = st.sidebar.slider('Helio-Longitude [deg]:',
                                0., 90., 45., step=0.1, format='%1.1f')
        extra_ = {'phi': phi * u.deg, 'r': r * u.R_sun, 'mode': 1}
    elif mode_connect == 'First Connection':
        r = np.logspace(np.log10(1),
                        np.log10(215), 250)
        phi = st.sidebar.slider('Helio-Longitude [deg]:',
                                0., 90., 45., step=0.1, format='%1.1f')
        extra_ = {'phi': phi * u.deg, 'r': r * u.R_sun, 'mode': 2}
    elif mode_connect == '-':
        extra_ = None

    st.sidebar.markdown('## Plot Options')
    with st.sidebar.expander('Select from different plots to visualize:',
                             expanded=True):
        plt_distprop = st.checkbox('Propagation plot', value=True)
        plt_conntime = st.checkbox('Connection times vs Phi plot', value=False)
        plt_paramete = st.checkbox('Parameters profile plot', value=False)
        plt_cormodel = st.checkbox('Coronal models plot', value=False)

    st.sidebar.markdown('## Coronal Parameters')
    st.sidebar.markdown('')
    with st.sidebar.expander('Magnetic fieldlines configuration:', expanded=False):
        fls_model_ = st.radio('Select between the models', ['Radial', 'Parker spiral', 'Streamer(*)'])
        st.markdown('(*) Mode not implemented yet')
        if fls_model_ == 'Streamer(*)':
            st.warning('This mode is not implemented yet.')
            st.stop()
    with st.sidebar.expander('Density models:', expanded=False):
        density_model_ = st.radio('Select between the models', ['Saito', 'Newkirk', 'Leblanc'])
        nfold_value = st.number_input('N-fold number', 0.5, 5., 1., step=0.1)
    with st.sidebar.expander('Magnetic field models:', expanded=False):
        magnetic_model_ = st.radio('Select between the models', ['1/r^2', '1/r^3'])
        B0_value = st.number_input('B0 at surface [gauss]', 0.5, 10., 2.2, step=0.1)
    with st.sidebar.expander('Solar wind models:', expanded=False):
        st.radio('Select between the models', ['Parker'])
        T0_value = st.number_input('T0 at surface  [MK]',
                                   10.**5, 10.**7, 1.4 * 10**6,
                                   step=10.**5, format='%e')
        usw_model = coronal_models.solarwind_speed_model(np.linspace(1, 215, 215) * u.R_sun, 'Parker', T0_value * u.Kelvin)
        usw_mean = np.mean(usw_model.vsw)
        st.write('Mean SW Speed= {:.2f}'.format(usw_mean))

    #################################################
    # Construct the model
    cmodel_options = {'density_model': {'model': density_model_, 'nfold': nfold_value},
                      'magnetic_model': {'model': magnetic_model_, 'topo': fls_model_,
                                         'B0': B0_value * u.gauss, 'usw': usw_mean},
                      'sw_model': {'model': 'Parker', 'T': T0_value * u.Kelvin},
                      }
    lon_mesh, distance_mesh = np.meshgrid(np.linspace(0, np.pi/2, 4*91)*u.rad,
                                          np.logspace(np.log10(1*u.R_sun.in_units(u.km)),
                                                      np.log10(215*u.R_sun.in_units(u.km)), 250) * u.km)

    cpoints = HDPm.connection_points(lon_mesh, distance_mesh)
    HDPm.coronal_parameters(cpoints, cmodel_options)
    smodel = HDPm.disturbance(a0*(u.m/u.second**2),
                              V0*(u.km/u.second),
                              alpha,
                              epsilon)
    # smodel = smodel.set_parameters(cpoints, cparams)

    plot_buffer = []

    if plt_distprop:
        st.subheader('Propagation Plot')
        left_col, right_col = st.columns(2)
        propa_plt_mode = left_col.selectbox('Select visualization',
                                            ['Disturbance', 'Fieldlines'])
        if propa_plt_mode == 'Disturbance':
            first_cmmode = 'Time'
        elif propa_plt_mode == 'Fieldlines':
            first_cmmode = 'Phi'
        cmmode_ = right_col.selectbox('Select colormap parameter',
                                      [first_cmmode, 'Theta_BN', 'Mach Alfven'])
        left_col, right_col = st.columns(2)
        xlim = left_col.slider('Adjust x-axis limits:', -3, 15, (0, 12))
        ylim = right_col.slider('Adjust y-axis limits:', -15, 15, (-6, 6))

        plt = utils_plot.plot_propagation(smodel, cmodel_options,
                                          app=True, plt_type=propa_plt_mode.lower(),
                                          extra=extra_, xlim=xlim, ylim=ylim, cmmode=cmmode_)

        # Download button
        plot2 = io.BytesIO()
        plt.savefig(plot2, format='png', bbox_inches='tight')
        download_button_str = download_button(plot2.getvalue(),
                                              'HDPmt_propagation_plot.png',
                                              'Download figure as .png file',
                                              pickle_it=False)
        st.markdown(download_button_str, unsafe_allow_html=True)
        plot_buffer.append(('HDPmt_propagation_plot.png', plot2))
        st.markdown('---')

    if plt_conntime:
        st.subheader('Connection times vs Phi plot')

        plt_MA = st.checkbox('Plot MA contour?', value=True)

        plt = utils_plot.plot_phivstime(smodel, cpoints, cmodel_options,
                                        plot_MA=plt_MA,
                                        app=True,
                                        extra=extra_)

        # Download button
        plot2 = io.BytesIO()
        plt.savefig(plot2, format='png', bbox_inches='tight')
        download_button_str = download_button(plot2.getvalue(),
                                              'HDPmt_connection_times_vs_phi.png',
                                              'Download figure as .png file',
                                              pickle_it=False)
        st.markdown(download_button_str, unsafe_allow_html=True)
        plot_buffer.append(('HDPmt_connection_times_vs_phi.png', plot2))
        st.markdown('---')

    if plt_paramete:
        st.subheader('Parameters profile plot')
        left_col, right_col = st.columns(2)
        mode_profile = left_col.selectbox('Select profile mode',
                                          ['Height', 'Time'])
        mode_param = right_col.selectbox('Select shock parameter',
                                         ['Mach Alfven', 'Theta_BN'])
        if mode_param == 'Mach Alfven':
            parameter_ = 'MA'
        elif mode_param == 'Theta_BN':
            parameter_ = 'ThBN'
        plt = utils_plot.plot_parameters_profile(smodel, cmodel_options,
                                                 app=True,
                                                 parameter_mode=parameter_,
                                                 xmode=mode_profile.lower(),
                                                 extra=extra_)
        # Download button png
        plot2 = io.BytesIO()
        plt.savefig(plot2, format='png', bbox_inches='tight')
        col1, col2 = st.columns(2)
        download_button_str = download_button(plot2.getvalue(),
                                              'HDPmt_parameters_profile_plot.png',
                                              'Download figure as .png file',
                                              pickle_it=False)  # TODO Use different filenames
        col1.markdown(download_button_str, unsafe_allow_html=True)
        plot_buffer.append(('HDPmt_parameters_profile_plot.png', plot2))
        # Download button json data
        # TODO Add here an option to download the parameter profile data to a json file
        if extra_ is not None:
            dict_smodel = utils_app.make_smodel_dict(smodel, cmodel_options)
            line = plt.gca().get_lines()[-2]  # TODO: If the plot order change this will result to error
            xd = line.get_xdata()
            yd = line.get_ydata()
            dict_smodel.update({
                'profile': {
                    mode_profile: list(xd),
                    parameter_: list(yd)
                }
            })
            json_buffer = io.BytesIO()

            json_buffer.write(json.dumps(dict_smodel, indent=' ', cls=JsonCustomEncoder).encode())
            download_button_str = download_button(json_buffer.getvalue(),
                                                  'HDPmt_model_parameters_profile.json',
                                                  'Download parameter profile as .json file')
            col2.markdown(download_button_str, unsafe_allow_html=True)
        st.markdown('---')

    if plt_cormodel:
        st.subheader('Coronal models plot')
        utils_plot.plot_coronal_models(smodel, cpoints, cmodel_options, app=True)
        # Download button
        plot2 = io.BytesIO()
        plt.savefig(plot2, format='png', bbox_inches='tight')
        download_button_str = download_button(plot2.getvalue(),
                                              'HDPmt_coronal_models.png',
                                              'Download figure as .png file',
                                              pickle_it=False)
        st.markdown(download_button_str, unsafe_allow_html=True)
        plot_buffer.append(('HDPmt_coronal_models.png', plot2))
        st.markdown('---')

    st.subheader('About HDPmt and this application')
    st.markdown("""
                The _Heliospheric Disturbance Propagation Model & Tool
                (HDPmt)_ is an open-source software package developed to model
                the propagation of heliospheric disturbances. _HDPm_ is a
                simple geometrical model to describe the propagation
                and longitudinal extension of a disturbance in the heliosphere
                and _HDPt_ is a lightweight tool that can be used to perform case studies
                and visualize the modeled results. The implementation of
                the Heliospheric Disturbance Propagation Model is described by
                [Kouloumvakos et al. (2021)](https://www.aanda.org/articles/aa/pdf/2023/01/aa44363-22.pdf).

                **Github**: You can find [here](https://github.com/AthKouloumvakos/HDPmt) the latest version of HDPmt.

                **Citation**: Please cite this software as [Kouloumvakos et al. (2021)](https://ui.adsabs.harvard.edu/abs/2023A%26A...669A..58K/exportcitation)
                """)
    col1, col2 = st.columns((5, 1))
    col1.markdown("The development of the online tool has received funding from the \
                   European Union's Horizon 2020 research and innovation programme \
                   under grant agreement No 101004159 (SERPENTINE).")
    col2.markdown('[<img src="https://serpentine-h2020.eu/wp-content/uploads/2021/02/SERPENTINE_logo_new.png" \
                   height="80">](https://serpentine-h2020.eu)', unsafe_allow_html=True)

    #############################################################
    # I/O sidebar section
    st.sidebar.markdown('## I/O model parameters & plots')

    dict_smodel = utils_app.make_smodel_dict(smodel, cmodel_options)
    json_buffer = io.BytesIO()
    json_buffer.write(json.dumps(dict_smodel, indent=' ', cls=JsonCustomEncoder).encode())
    download_button_str = download_button(json_buffer.getvalue(),
                                          'HDPmt_model_parameters.json',
                                          'Download parameters as .json file')
    st.sidebar.markdown(download_button_str, unsafe_allow_html=True)
    plot_buffer.append(('HDPmt_model_parameters.json', json_buffer))

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
        for file_name, data in plot_buffer:
            zip_file.writestr(file_name, data.getvalue())
    download_button_str = download_button(zip_buffer.getvalue(),
                                          'HDPmt_plots.zip',
                                          'Download param & plots as .zip file',
                                          pickle_it=False)
    st.sidebar.markdown(download_button_str, unsafe_allow_html=True)

    # TODO: Provide parameters as json input
    # json_input = st.sidebar.file_uploader('Upload parameters as .json file',
    #                                      type=['.json'], on_change=None)


if __name__ == '__main__':
    run()
