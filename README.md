# Heliospheric Disturbance Propagation Model & Tool (HDPmt)

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Version](https://img.shields.io/github/v/release/AthKouloumvakos/HDPmt)](https://github.com/AthKouloumvakos/HDPmt/releases)
[![Release Date](https://img.shields.io/github/release-date/AthKouloumvakos/HDPmt)](https://github.com/AthKouloumvakos/HDPmt/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![flake8](https://github.com/AthKouloumvakos/HDPmt/actions/workflows/flake8.yml/badge.svg)
![pytest](https://github.com/AthKouloumvakos/HDPmt/actions/workflows/pytest.yml/badge.svg)

_HDPm_ is a simple geometrical model to describe the propagation and longitudinal extension of a disturbance
in the heliosphere and _HDPt_ is a lightweight tool that can be used to visualize the model results and perform
case studies. The implementation of the Heliospheric Disturbance Propagation Model is described by Kouloumvakos et al. (2021).
A preview of this tool is available in [![Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://hpdm-tool.streamlitapp.com).

## 💾 Installation

_HDPmt_ is written in Python >=3.8 and has some package requirements, which are listed in the requirements.txt and environment.yml files.
To run locally this application we recommend to create its own virtual environment in python.

**Recommended (conda)**

We create the virtual environment (see [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)) and install the dependent packages by doing the following in the terminal,

```python
conda env create -f environment.yml
conda info --envs
```

Then every time before using _HDPmt_, you have to activate the environment and when finishing your work deactivate it using the following commands,

```python
# Activate the enviroment
conda activate HPDmt

# Here you run _HDPmt_ (see Run locally section)

# When you are done you can deactivate a virtual environment
conda deactivate
```

**Alternative (pip)**

You can create a virtual environment in Python inside the _HDPmt_ project (root) folder using ```pip``` and by doing the following in the terminal,

```python
# Create the virtual environment in _HDPmt_'s root folder
python3 -m venv env

# Activate the environment
source env/bin/activate

# install the required packages using pip3
pip3 install -r requirements.txt

# When you are done you can deactivate a virtual environment
deactivate
```

Now you can run any part of the _HDPmt_ (see Run locally section).

You may also add your _HDPmt_ directory to the environment variable ```PYTHONPATH```. This is useful if you need to run _HDPmt_ tests or when you need to run some of the package modules out of streamlit.

In the terminal use the following and change the \<HPDmtRootDir\> with your path.

```
export PYTHONPATH="${PYTHONPATH}:<HPDmtRootDir>/HPDmt"
```

For a permanent solution, if you're using bash (on a Mac or GNU/Linux distribution), add the above line to your ```~/.bashrc``` file (changing the \<HPDmtRootDir\> with your path first).

## 🐾 Run locally the HDPmt application

Install the required python packages, and then run the application with streamlit. Open a terminal to _HDPmt_ directory and run

```
# Run the application using streamlit
streamlit run HDPt_app.py
```
The application should now open in the default browser!

## 📙 Usage

Complete documentation of the _HDPmt_ can be found in (under construction)

## Acknowledging or Citing HDPmt

If you use HDPmt for work or research presented in a publication, please cite this package by including in the Acknowledgement section the following: "This research has made use of HDPmt v?.?.? (fill the version), an open-source and free Python package (include zenodo citation here) developed to model the propagation of a heliospheric disturbance.". Include the HDPmt citation in your paper's reference list.

## 📦 Usefull python packages:

- [SunPy](https://sunpy.org/): The community-developed, free and open-source solar data analysis environment for Python.
- [AstroPy](https://www.astropy.org/): The Astropy Project is a community effort to develop a single core package for Astronomy in Python.

-----

The development of the online tool has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 101004159 (SERPENTINE).

[<img src="https://serpentine-h2020.eu/wp-content/uploads/2021/02/SERPENTINE_logo_new.png" height="100">](https://serpentine-h2020.eu)
