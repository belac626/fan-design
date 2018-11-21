"""Determines fan flow and geometric properties."""
import math as m
import pathlib as pth
from configparser import ConfigParser
from math import degrees as d
from math import radians as r

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
class FanFlow():
    """Set flow properties."""

    def __init__(self):
        """Instantiate flow properties."""
        self.reaction = 0
        self.phi = 0
        self.psi = 0
        self.rpm = 0
        self.radius = 0
        self.span = 0

        self.u = 0
        self.beta1 = 0
        self.beta2 = 0
        self.betam = 0
        self.alpha1 = 0
        self.alpha2 = 0
        self.alpham = 0
        self.cx = 0
        self.w1 = 0
        self.w2 = 0
        self.c1 = 0
        self.c2 = 0
        self.wt1 = 0
        self.wt2 = 0
        self.ct1 = 0
        self.ct2 = 0

    def GetFanFlowConfig(self, filename:str, station:str):
        """Read .ini file."""
        cfp = ConfigParser()
        file = pth.Path(f'Config/{filename}')
        cfp.read(file)  # Fan_Stage.ini

        self.reaction = cfp.getfloat('flow', 'r')
        self.phi = cfp.getfloat('flow', 'phi')
        self.rpm = cfp.getfloat('flow', 'rpm')
        self.alpha1 = cfp.getfloat('flow', 'alpha1')
        self.span = (cfp.getfloat('tip', 'radius')
                     - cfp.getfloat('root', 'radius'))
        if station == 'mean':
            self.radius = m.sqrt((cfp.getfloat('root', 'radius')**2
                                  + cfp.getfloat('tip', 'radius')**2)/2)
        else:
            self.radius = cfp.getfloat(station, 'radius')

    def CalcFanFlow(self, iphi=None):
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
class Stage():
    """Set stage properties."""

    def __init__(self):
        """Instantiate stage properties."""
        self.root = FanFlow()
        self.mean = FanFlow()
        self.tip = FanFlow()

        self.rpm = 0
        self.span = 0

    def GetStageConfig(self, filename: str):
        """Read .inp files."""
        self.root.GetFanFlowConfig(filename=filename, station='root')
        self.mean.GetFanFlowConfig(filename=filename, station='mean')
        self.tip.GetFanFlowConfig(filename=filename, station='tip')

        cfp = ConfigParser()
        file = pth.Path(f'Config/{filename}')
        cfp.read(file)

        self.rpm = cfp.getfloat('flow', 'rpm')
        self.span = (cfp.getfloat('tip', 'radius')
                     - cfp.getfloat('root', 'radius'))

    def CalcStage(self):
        """Calculate stage properties (assumes cosntant dcx/dr)."""
        self.mean.CalcFanFlow()

        self.root.CalcFanFlow()
        lphi = 1e-6
        hphi = 2 - 1e-6
        rootPhi = self.root.phi
        while (abs(self.mean.cx - self.root.cx) > 1e-3):
            if (self.root.cx < self.mean.cx):
                lphi = rootPhi
            else:
                hphi = rootPhi
            rootPhi = (hphi + lphi)/2
            self.root.CalcFanFlow(iphi=rootPhi)

        self.tip.CalcFanFlow()
        lphi = 1e-6
        hphi = 2 - 1e-6
        tipPhi = self.tip.phi
        while (abs(self.mean.cx - self.tip.cx) > 1e-3):
            if (self.tip.cx < self.mean.cx):
                lphi = tipPhi
            else:
                hphi = tipPhi
            tipPhi = (hphi + lphi)/2
            self.tip.CalcFanFlow(iphi=tipPhi)
