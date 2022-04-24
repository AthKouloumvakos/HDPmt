import os

from setuptools import find_packages, setup
from setuptools.config import read_configuration

extras = read_configuration('setup.cfg')

# create home directory
if not os.path.isdir(os.path.join(os.environ['HOME'], '.HDPmt')):
    os.mkdir(os.path.join(os.environ['HOME'], '.HDPmt'))

with open('README.md', 'r') as f:
    long_description = f.read()

# get requirements
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='HDPmt',
    use_scm_version={'write_to': os.path.join('HDPmt', '_version.py')},
    # version="0.4.0",
    description='Heliospheric Disturbance Propagation Model & Tool (HDPmt): An open-source software package developed to model the propagation of heliospheric disturbances.',
    url='https://github.com/AthKouloumvakos/HDPmt',
    # long_description=long_description,
    long_description_content_type='text/markdown',

    author='Athanasios Kouloumvakos',
    author_email='athkouloumvakos@gmail.com',
    license='GPL-3.0',
    # license_file="LICENSE.md",

    python_requires='>=3.8',
    install_requires=requirements,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'hdpmt = HDPmt.HDPt_cli:main'
        ]
    },
    extras_require=extras,
)
