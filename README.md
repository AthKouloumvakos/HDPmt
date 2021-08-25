# Heliospheric Disturbance Propagation Model & Tool (HDPmt)

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Version](https://img.shields.io/github/v/release/AthKouloumvakos/HDPmt)](https://github.com/AthKouloumvakos/HDPmt/releases)
[![Release Date](https://img.shields.io/github/release-date/AthKouloumvakos/HDPmt)](https://github.com/AthKouloumvakos/HDPmt/releases)
[![Downloads](https://img.shields.io/github/downloads/AthKouloumvakos/HDPmt/total)](https://github.com/AthKouloumvakos/HDPmt)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

_HDPm_ is a simple geometrical model to describe the propagation and longitudinal extension of a disturbance
in the heliosphere and _HDPt_ is a lightweight tool that can be used to visualize the model results and perform
case studies. The implementation of the Heliospheric Disturbance Propagation Model is described by Kouloumvakos et al. (2021).
A preview of this tool is available at [https://athkouloumvakos.github.io/HDPmt](https://athkouloumvakos.github.io/HDPmt).

## 💾 Installation

_HDPmt_ is written in Python >=3.8 and has some package requirements, which are listed in the requirements.txt and environment.yml files. 
To run localy this application we recomend to create its own virtual enviroment in python.

**Recomended (conda)**

Because of a range of dependencies that packages have, the simplest way to work with _HDPmt_ 
is in conda and to create its own environment using the ```conda env create```. 
If you already have conda installed, then with the prompt active you can,
```cd``` the root directory of _HDPmt_ and in your terminal do:

```python
# Create a virtual environment in python using conda
#  see https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
conda env create -f environment.yml
conda info --envs

# Activate the enviroment
conda activate hdpmt
```

**Alternative (pip)**

```python
# You can create a virtual environment in python inside the project folder.
#  see https://docs.python.org/3/library/venv.html
python3 -m venv env

# Activate the enviroment
source env/bin/activate

# install the required packages using pip3
pip3 install -r requirements.txt

# You can deactivate a virtual environment by 
#  typing “deactivate” in your shell.
deactivate
```

Now you can run any part of the HDPmt (see Usage section).

## 🐾 Run localy the HDPmt application
Install the required python packages, and then run the application with streamlit. 
```
# Run the application using streamlit
streamlit run HDPt_app.py
```
The application should now open in the default browser!

## 📙 Usage

Complete documentation of the _HDPmt_ can be found in HDPmt_doc.ipynb.

## Acknowledging or Citing HDPmt

If you use HDPmt for work or research presented in a publication, please cite this package by including in the Acknowledgement section the following: "This research has made use of HDPmt v?.?.? (fill the version), an open-source and free Python package (include zenodo citation here) developed to model the propagation of a heliospheric disturbance.". Include the HDPmt citation in your paper's reference list.

## 📦 Usefull python packages:
        
- [SunPy](https://sunpy.org/): The community-developed, free and open-source solar data analysis environment for Python.
- [AstroPy](https://www.astropy.org/): The Astropy Project is a community effort to develop a single core package for Astronomy in Python.

-----

The development of the online tool has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 101004159 (SERPENTINE). 

[<img src="https://serpentine-h2020.eu/wp-content/uploads/2021/02/SERPENTINE_logo_new.png" height="100">](https://serpentine-h2020.eu)
