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


import datetime


def make_smodel_dict(smodel, cmodel_options):
    dict_smodel = {
        'disturbance': {
            'a0': smodel.a0,
            'V0': smodel.V0,
            'alpha': smodel.alpha,
            'epsilon': smodel.epsilon,
        },
        'coronal_model': cmodel_options,
        'misc': {
            'date_run': (datetime.datetime.now()).strftime('%Y-%m-%d %H:%M:%S'),
        },
    }
    return dict_smodel
