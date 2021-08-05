
import astropy.units as u
import datetime

def make_smodel_dict(smodel):
    dict_smodel ={
    "model": {
        "a0": (smodel.a0.to_value(u.km/u.s**2)).tolist(),
        "V0": (smodel.V0.to_value(u.km/u.s)).tolist(),
        "alpha": smodel.alpha,
        "epsilon": smodel.epsilon,
        "density_model": smodel.coronal_models.density_model.type,
        "NFold": smodel.coronal_models.density_model.NFold,
        "magnetic_model": smodel.coronal_models.magnetic_model.type,
        "B0": (smodel.coronal_models.magnetic_model.B0.to_value(u.gauss)).tolist(),
        "sw_model": smodel.coronal_models.sw_model.type,
        "T0": (smodel.coronal_models.sw_model.T0.to_value(u.megaKelvin)).tolist()
    },
    "misc": {
        "date_run": (datetime.datetime.now()).strftime('%Y-%m-%d %H:%M:%S'),
        "units": {
                  "a0": "km/s2",
                  "V0": "km/s",
                  "alpha": "",
                  "epsilon": "",
                  "B0": "gauss", 
                  "T0": "MK",
                  }
    },
    }
    return dict_smodel
