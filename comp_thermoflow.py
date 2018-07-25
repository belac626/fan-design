"""Determines compressor thermodynamic and flow properties."""
import math as m
import os
from configparser import ConfigParser

################################
# #Class: ThermoFlow
# Holds stage flow characteristics


# ---Flow Attributes:
# r: reaction coefficient
# phi: loading coefficient
# psi: flow coefficient
# rpm: shaft speed (rpm)
# radius: station in blade (in)

# u: blade speed (in/s)
# cx: axial air velocity (in/s)
# beta1: relative air inlet angle (deg)
# beta2: relative air exit angle (deg)
# alpha1: absolute air inlet angle (deg)
# alpha2: absolute air exit angle (deg)
# w1: relative air inlet speed (in/s)
# w2: relative air exit speed (in/s)
# c1: absolute air inlet speed (in/s)
# c2: absolute air exit speed (in/s)
# wt1: tangential relative air inlet speed (in/s)
# wt2: tangential relative air outlet speed (in/s)


# ---Thermodynamic Attributes
# M: axial mach number
# T1: inlet static temperature (K)
# T2: exit static temperature (K)
# T01: inlet stagnation temperature (K)
# T02: exit stagnation temperature (K)
# p1: inlet static pressure (bar)
# p2: exit static pressure
# p01: inlet stagnation pressure
# p02: exit stagnation pressure

# rho: density (air) (kg/m^3)
# cp: constant pressure specific heat (air) (J/kg*K)
# cv: constant volume specific heat (air) (J/kg*K)
# gamma: specific heat ratio (air)
# R: ideal gas constant (air) (J/kg*K)

# dT0: change in stagnation temperature
# PR: pressure ratio
# Cp: static pressure rise coefficient
# W: ideal work required
################################
class ThermoFlow():
    """Set flow properties."""

    def __init__(self):
        """Instantiate flow properties."""
        self.r = 0
        self.phi = 0
        self.psi = None
        self.rpm = 0
        self.radius = 0

        self.u = 0
        self.beta1 = 0
        self.beta2 = 0
        self.alpha1 = 0
        self.alpha2 = 0
        self.cx = 0
        self.w1 = 0
        self.w2 = 0
        self.c1 = 0
        self.c2 = 0
        self.wt1 = 0
        self.wt2 = 0
        self.ct1 = 0
        self.ct2 = 0

        self.M = 0

        self.T1 = 0
        self.T2 = 0
        self.T01 = 0
        self.T02 = 0
        self.p1 = 0
        self.p2 = 0
        self.p01 = 0
        self.p02 = 0

        self.rho = 1.2
        self.mu = 1.8*10**(-5)
        self.Re = 0
        self.cp = 1005
        self.cv = 718
        self.gamma = self.cp/self.cv
        self.R = self.cp - self.cv

        self.dT0 = 0
        self.PR = 0
        self.Cp = 0
        self.W = 0
        self.tau = 0

    def GetThermoFlowConfig(self, filename, station):
        """Read .inp file."""
        cfp = ConfigParser()

        os.chdir(r'.\Config')
        cfp.read(filename)

        self.r = cfp.getfloat('flow', 'r')
        self.phi = cfp.getfloat('flow', 'phi')
        self.psi = None if cfp.get('flow', 'psi') == 'None' \
            else cfp.getfloat('flow', 'psi')
        self.rpm = cfp.getfloat('flow', 'rpm')
        self.alpha1 = None if cfp.get('flow', 'alpha1') == 'None' \
            else cfp.getfloat('flow', 'alpha1')

        self.p1 = cfp.getfloat('thermo', 'p1')
        self.T1 = cfp.getfloat('thermo', 'T1')

        if station == 'mean':
            self.radius = (cfp.getfloat('root', 'radius')
                           + cfp.getfloat('tip', 'radius'))/2
        else:
            self.radius = cfp.getfloat(station, 'radius')

        os.chdir('..')

    def CalcThermoFlow(self, iphi=None):
        """Calculate flow angles."""
        # Flow
        if self.psi is None:
            self.psi = 2*(1 - self.r - self.phi*m.tan(m.radians(self.alpha1)))

        if iphi is not None:
            self.phi = iphi

        self.u = self.rpm/60*2*m.pi*self.radius
        self.beta1 = m.degrees(m.atan(self.psi/self.phi
                                      + (self.r - self.psi/2)/self.phi))
        self.beta2 = m.degrees(m.atan((self.r - self.psi/2)/self.phi))
        if self.alpha1 is None:
            self.alpha1 = m.degrees(m.atan((1 - self.phi/2 - self.r)/self.phi))

        self.alpha2 = m.degrees(m.atan(1/self.phi
                                       - m.tan(m.radians(self.beta2))))
        self.cx = self.phi*self.u
        self.w1 = self.cx/m.cos(m.radians(self.beta1))
        self.w2 = self.cx/m.cos(m.radians(self.beta2))
        self.c1 = self.cx/m.cos(m.radians(self.alpha1))
        self.c2 = self.cx/m.cos(m.radians(self.alpha2))
        self.wt1 = self.w1*m.sin(m.radians(self.beta1))
        self.wt2 = self.w1*m.sin(m.radians(self.beta2))
        self.ct1 = self.c1*m.sin(m.radians(self.alpha1))
        self.ct2 = self.c2*m.sin(m.radians(self.alpha2))

        # Thermodynamics (converting all english measurements to metric)
        self.M = self.cx*0.0254/m.sqrt(self.gamma*self.R*self.T1)
        self.dT0 = self.psi*(self.u*0.0254)**2/self.cp

        self.T01 = self.T1*(1 + ((self.gamma - 1)*self.M**2)/2)
        self.p01 = self.p1*(self.T01/self.T1)**(self.gamma/(self.gamma - 1))
        self.T02 = self.T01 + self.dT0
        self.p02 = self.p01*(self.T02/self.T01)**(self.gamma/(self.gamma - 1))
        self.T2 = self.T02/(1 + ((self.gamma - 1)*self.M**2)/2)
        self.p2 = self.p02/(self.T02/self.T2)**(self.gamma/(self.gamma - 1))

        self.PR = (1 + (self.dT0/self.T01))**(self.gamma*1/(self.gamma-1))
        self.Cp = (self.p2 - self.p1)/(self.p01 - self.p1)
        self.W = self.cp*self.dT0
        self.tau = self.W/(self.u/self.radius)
        self.capacity = self.gamma/m.sqrt(self.gamma - 1)*self.M\
            * (1 + (self.gamma - 1)/2*self.M**2)\
            ** (-(1/2)*(self.gamma + 1)/(self.gamma - 1))


################################
# #Class: Stage
# Holds stage characteristics
# #Attributes:
# root: root (Class: ThermoFlow)
# mean: mean (Class: ThermoFlow)
# tip: tip (Class: ThermoFlow)
# rootRadius: location of root flow passage (in)
# tipRadius: location of tip flow passage (in)
################################
class Stage():
    """Set stage properties."""

    def __init__(self):
        """Instantiate stage properties."""
        self.root = ThermoFlow()
        self.mean = ThermoFlow()
        self.tip = ThermoFlow()

        self.rho = 1.177
        self.cp = 1005
        self.cv = 718
        self.gamma = self.cp/self.cv
        self.R = self.cp - self.cv

        self.m = 0
        self.dT0 = 0
        self.PR = 0
        self.Cp = 0
        self.W = 0
        self.tau = 0

    def GetStageConfig(self, filename):
        """Read .inp files."""
        self.root.GetThermoFlowConfig(filename=filename, station='root')
        self.mean.GetThermoFlowConfig(filename=filename, station='mean')
        self.tip.GetThermoFlowConfig(filename=filename, station='tip')

    def CalcStage(self):
        """Calculate stage properties."""
        self.mean.CalcThermoFlow()

        self.root.CalcThermoFlow()
        lphi = 1e-6
        hphi = 2 - 1e-6
        rootPhi = self.root.phi
        while (abs(self.mean.cx - self.root.cx) > 1e-3):
            if (self.root.cx < self.mean.cx):
                lphi = rootPhi
            else:
                hphi = rootPhi
            rootPhi = (hphi + lphi)/2
            self.root.CalcThermoFlow(iphi=rootPhi)

        self.tip.CalcThermoFlow()
        lphi = 1e-6
        hphi = 2 - 1e-6
        tipPhi = self.tip.phi
        while (abs(self.mean.cx - self.tip.cx) > 1e-3):
            if (self.tip.cx < self.mean.cx):
                lphi = tipPhi
            else:
                hphi = tipPhi
            tipPhi = (hphi + lphi)/2
            self.tip.CalcThermoFlow(iphi=tipPhi)

        self.dT0 = self.mean.dT0
        self.PR = self.mean.PR
        self.Cp = self.mean.Cp
        self.W = self.mean.W/1000
        self.tau = self.mean.tau
        self.m = self.rho*self.mean.cx*0.0254 \
            * m.pi*((self.tip.radius*0.0254)**2
                    - (self.root.radius*0.0254)**2)
        print('\n')
        print('Mean Stage Properties')
        print('    m  : {var} kg/s'.format(var=round(self.m, 2)))
        print('    dT0: {var} K'.format(var=round(self.dT0, 2)))
        print('    PR : {var}'.format(var=round(self.PR, 2)))
        print('    Cp : {var}'.format(var=round(self.Cp, 2)))
        print('    W  : {var} kW'.format(var=round(self.W*self.m, 2)))
        print('    torque: {var} Nm'.format(var=round(self.tau*self.m, 2)))
