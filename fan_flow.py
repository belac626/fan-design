"""Determines fan flow and geometric properties."""
<<<<<<< HEAD
import math as m
import pathlib as pth
from configparser import ConfigParser
from math import degrees as d
from math import radians as r
from typing import ClassVar

from dataclasses import dataclass
=======
>>>>>>> b1d1aabfeebac9ac45a2b06c3a2bedd09093bad0

# pylint: disable=R0902
# - allow more than 7 instance variables
# pylint: disable=R0913
# - allow more than 5 arguments to a function
# pylint: disable=R0914
# - allow more than 15 local variables
# pylint: disable=R0915
# - allow more than 50 statemnts in class or function
# pylint: disable=C0103
# - allow variable names less than 3 characters
# pylint: disable=C0412
# - allow imports to be sorted by isort and not grouped

import math as m
import pathlib as pth
from configparser import ConfigParser
from math import degrees as d
from math import radians as r
from typing import ClassVar

from dataclasses import dataclass

################################
# #Class: ThermoFlow
# Holds stage flow characteristics

# ---Flow Attributes:
# r: reaction coefficient
# phi: loading coefficient
# psi: flow coefficient
# rpm: shaft speed (rpm)
# radius: station in blade (m)
# span: blade length (m)

# u: blade speed (m/s)
# cx: axial air velocity (m/s)
# beta1: relative air inlet angle (deg)
# beta2: relative air exit angle (deg)
# betam: atan(0.5*(tan(beta1) + tan(beta2)) (deg)
# alpha1: absolute air inlet angle (deg)
# alpha2: absolute air exit angle (deg)
# alpham: atan(0.5*(tan(alpha1) + tan(alpha2)) (deg)
# w1: relative air inlet speed (m/s)
# w2: relative air exit speed (m/s)
# c1: absolute air inlet speed (m/s)
# c2: absolute air exit speed (m/s)
# c1: absolute air inlet speed (m/s)
# c2: absolute air exit speed (m/s)
# wt1: tangential relative air inlet speed (m/s)
# wt2: tangential relative air outlet speed (m/s)

# T1: inlet static temperature (K)
# P1: inlet static pressure (bar)
################################


@dataclass
class FanFlow:
    """Set flow properties."""

    reaction: float = 0
    phi: float = 0
    psi: float = 0
    rpm: float = 0
    radius: float = 0
    span: float = 0

    u: float = 0
    beta1: float = 0
    beta2: float = 0
    betam: float = 0
    alpha1: float = 0
    alpha2: float = 0
    alpham: float = 0
    cx: float = 0
    w1: float = 0
    w2: float = 0
    c1: float = 0
    c2: float = 0
    wt1: float = 0
    wt2: float = 0
    ct1: float = 0
    ct2: float = 0

    def get_fan_flow_config(self, filename: str, station: str):
        """Read .ini file."""
        file = pth.Path(f'{pth.Path.cwd()}/Config/{filename}')
        config = ConfigParser()
        config.read(file)  # Fan_Stage.ini

        self.reaction = config.getfloat('flow', 'r')
        self.phi = config.getfloat('flow', 'phi')
        self.rpm = config.getfloat('flow', 'rpm')
        self.alpha1 = config.getfloat('flow', 'alpha1')
        self.span = (config.getfloat('tip', 'radius')
                     - config.getfloat('root', 'radius'))
        if station == 'mean':
            self.radius = m.sqrt((config.getfloat('root', 'radius')**2
                                  + config.getfloat('tip', 'radius')**2)/2)
        else:
            self.radius = config.getfloat(station, 'radius')

    def calc_fan_flow(self, iphi=None):
        """Calculate flow angles."""
        if iphi is not None:
            self.phi = iphi
        self.psi = 2*(1 - self.reaction - self.phi*m.tan(r(self.alpha1)))
        self.u = self.rpm/60*2*m.pi*self.radius
        self.beta1 = d(m.atan((1/self.phi) - m.tan(r(self.alpha1))))
        self.beta2 = d(m.atan((2*self.reaction-1)/self.phi
                              + m.tan(r(self.alpha1))))
        self.betam = d(m.atan((m.tan(r(self.beta1))
                               + m.tan(r(self.beta2)))/2))
        self.alpha2 = d(m.atan((1/self.phi) - m.tan(r(self.beta2))))
        self.alpham = d(m.atan((m.tan(r(self.alpha1))
                                + m.tan(r(self.alpha2)))/2))
        self.cx = self.phi*self.u
        self.w1 = self.cx/m.cos(r(self.beta1))
        self.w2 = self.cx/m.cos(r(self.beta2))
        self.c1 = self.cx/m.cos(r(self.alpha1))
        self.c2 = self.cx/m.cos(r(self.alpha2))
        self.wt1 = self.w1*m.sin(r(self.beta1))
        self.wt2 = self.w1*m.sin(r(self.beta2))
        self.ct1 = self.c1*m.sin(r(self.alpha1))
        self.ct2 = self.c2*m.sin(r(self.alpha2))

################################
# #Class: Stage
# Holds stage characteristics
# #Attributes:
# root: root (Class: ThermoFlow)
# mean: mean (Class: ThermoFlow)
# tip: tip (Class: ThermoFlow)

# rpm: angular velocity of rotor (float) (rpm)
# span: blade height (float) (m)
################################


@dataclass
class Stage:
    """Set stage properties."""

    root: ClassVar = FanFlow()
    mean: ClassVar = FanFlow()
    tip: ClassVar = FanFlow()

    rpm: float = 0
    span: float = 0

    def get_stage_config(self, filename: str):
        """Read .inp files."""
        self.root.get_fan_flow_config(filename=filename, station='root')
        self.mean.get_fan_flow_config(filename=filename, station='mean')
        self.tip.get_fan_flow_config(filename=filename, station='tip')

        file = pth.Path(f'Config/{filename}')
        config = ConfigParser()
        config.read(file)

        self.rpm = config.getfloat('flow', 'rpm')
        self.span = (config.getfloat('tip', 'radius')
                     - config.getfloat('root', 'radius'))

    def calcstage(self):
        """Calculate stage properties (assumes cosntant dcx/dr)."""
        self.mean.calc_fan_flow()

        self.root.calc_fan_flow()
        lphi = 1e-6
        hphi = 2 - 1e-6
        root_phi = self.root.phi
        while abs(self.mean.cx - self.root.cx) > 1e-3:
            if self.root.cx < self.mean.cx:
                lphi = root_phi
            else:
                hphi = root_phi
            root_phi = (hphi + lphi)/2
            self.root.calc_fan_flow(iphi=root_phi)

        self.tip.calc_fan_flow()
        lphi = 1e-6
        hphi = 2 - 1e-6
        tip_phi = self.tip.phi
        while abs(self.mean.cx - self.tip.cx) > 1e-3:
            if self.tip.cx < self.mean.cx:
                lphi = tip_phi
            else:
                hphi = tip_phi
            tip_phi = (hphi + lphi)/2
            self.tip.calc_fan_flow(iphi=tip_phi)
