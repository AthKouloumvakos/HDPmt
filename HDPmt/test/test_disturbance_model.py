"""
Test the disturbance model
"""

import astropy.units as u
import pytest_check as check
from astropy.tests.helper import assert_quantity_allclose

from HDPmt.HDPm import disturbance


def test_disturbance_parameters_unit():
    """
    Tests if the disturbance parameters have the correct units
    """
    dst = disturbance()
    check.equal(dst.V0.unit, u.km/u.second)
    check.equal(dst.a0.unit, u.km/u.second**2)
    check.equal((dst.alpha * u.dimensionless_unscaled).unit, u.dimensionless_unscaled)
    check.equal((dst.epsilon * u.dimensionless_unscaled).unit, u.dimensionless_unscaled)


def test_disturbance_parameters_values():
    """
    Tests if the disturbance parameters have the correct units
    """
    dst = disturbance()
    check.equal(dst.V0.value, 800)
    check.equal(dst.a0.value, 0)
    check.equal(dst.alpha, 0.5)
    check.equal(dst.epsilon, 1)


def test_disturbance_center_atT0():
    """
    Tests if the disturbance center is at the Sun surface for t0
    """
    dst = disturbance()
    x, y = dst.propagate(time=0*u.second)
    assert_quantity_allclose(x**2 + y**2, (1*u.R_sun)**2, atol=(5e-10*u.R_sun)**2)


def test_disturbance_center_distance():
    """
    Tests if the calculation of the disturbance's center is correct
    """
    dst = disturbance(V0=1*(u.R_sun/u.second))
    dsh = dst.dsh(time=1*u.second)
    check.equal(dsh.to_value(u.R_sun), 1.5)
    dst = disturbance(V0=1*(u.R_sun/u.second), alpha=1)
    dsh = dst.dsh(time=1*u.second)
    check.equal(dsh.to_value(u.R_sun), 2)


def test_disturbance_iscircular():
    """
    Tests if the disturbance is a circle for epsilon eq. 1
    """
    dst = disturbance(V0=1*(u.R_sun/u.second))
    x, y = dst.propagate(time=1*u.second)
    dsh = dst.dsh(time=1*u.second)
    assert_quantity_allclose((x - dsh)**2 + y**2, (1*u.R_sun)**2, atol=(5e-10*u.R_sun)**2)
