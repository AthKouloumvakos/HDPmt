
import astropy.units as u
import datetime

def make_smodel_dict(smodel, cmodel_options):
    dict_smodel ={
    "disturbance": {
        "a0": smodel.a0,
        "V0": smodel.V0,
        "alpha": smodel.alpha,
        "epsilon": smodel.epsilon,
     },
    'coronal_model': cmodel_options,
    "misc": {
        "date_run": (datetime.datetime.now()).strftime('%Y-%m-%d %H:%M:%S'),
     },
    }
    return dict_smodel
