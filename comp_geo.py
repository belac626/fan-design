"""Determines compressor blade geometry."""
import math as m
import os
from configparser import ConfigParser

import utils as u


################################
# #Class: Airfoil
# Gets and calculates variables needed to parameterize BP3333 airfoil
# #Attributes: (all linear dimensions normalized to chord length)
# xcr: stagger ratio (also very nearly xc) (1 = vle, 0 = vte) (float)
# ycr: max camber altitude ratio (1 = triangular camber, 0 = no camber) (float)
# xc: location of max camber (float)
# yc: max camber ratio (float)
# kc: curvature of camber crest (float)
# xt: location of max thickness (float)
# yt: max thickness ratio (float)
# kt: curvature of thickness crest (float)
# vle: inlet velocity angle (float)
# vte: exit velocity angle (float)
# cle: leading edge camber angle (for solidworks, from chord) (float) (deg)
# cte: trailing edge camber angle (for solidworks, from chord) (float) (deg)
# wte: trailing edge wedge angle (float) (deg) (NACA definition)
# rle: leading edge radius (float) (NACA definition)
################################
class Airfoil():
    """Set Bezier-PARSEC variables."""

    def __init__(self):
        """Instantiate airfoil properties."""
        self.N = 0
        self.AR = 0
        self.DF = 0
        self.DH = 0

        self.xcr = 0
        self.ycr = 0
        self.xc = 0
        self.yc = 0
        self.kc = 0
        self.bc = 0
        self.xt = 0
        self.yt = 0
        self.kt = 0
        self.bt = 0
        self.vle = 0
        self.vte = 0
        self.cle = 0
        self.cte = 0
        self.wte = 0
        self.rle = 0

        self.chord = 0
        self.sigma = 0
        self.camber = 0
        self.incidence = 0
        self.deviation = 0
        self.stagger = 0

    def GetAirfoilConfig(self, filename, station):
        """Read .inp file."""
        cfp = ConfigParser()

        os.chdir(r'.\Config')
        cfp.read(filename)

        self.N = cfp.getfloat('blade', 'N')
        self.AR = cfp.getfloat('blade', 'AR')

        self.xcr = cfp.getfloat(station, 'xcr')
        self.ycr = cfp.getfloat(station,  'ycr')
        self.xt = cfp.getfloat(station, 'xt')
        self.yt = cfp.getfloat(station, 'yt')

        os.chdir('..')

    def CalcAirfoil(self, stage, blade, station):
        """Calculate Bezier-PARSEC variables and airfoil properties."""
        location = {
            ('Rotor', 'root'): (stage.root.beta1,
                                stage.root.beta2,
                                stage.root.w1,
                                stage.root.w2,
                                stage.root.radius),
            ('Rotor', 'mean'): (stage.mean.beta1,
                                stage.mean.beta2,
                                stage.mean.w1,
                                stage.mean.w2,
                                stage.mean.radius),
            ('Rotor', 'tip'): (stage.tip.beta1,
                               stage.tip.beta2,
                               stage.tip.w1,
                               stage.tip.w2,
                               stage.tip.radius),
            ('Stator', 'root'): (stage.root.alpha2,
                                 stage.root.alpha1,
                                 stage.root.c2,
                                 stage.root.c1,
                                 stage.root.radius),
            ('Stator', 'mean'): (stage.mean.alpha2,
                                 stage.mean.alpha1,
                                 stage.mean.c2,
                                 stage.mean.c1,
                                 stage.mean.radius),
            ('Stator', 'tip'): (stage.tip.alpha2,
                                stage.tip.alpha1,
                                stage.tip.c2,
                                stage.tip.c1,
                                stage.tip.radius),
            ('IGV', 'root'): (0,
                              stage.root.alpha1,
                              stage.root.cx,
                              stage.root.c1,
                              stage.root.radius),
            ('IGV', 'mean'): (0,
                              stage.mean.alpha1,
                              stage.mean.cx,
                              stage.mean.c1,
                              stage.mean.radius),
            ('IGV', 'tip'): (0,
                             stage.tip.alpha1,
                             stage.tip.cx,
                             stage.tip.c1,
                             stage.tip.radius),
            }
        self.vle, self.vte, v1, v2, radius = location[(blade, station)]
        self.DH = v2/v1
        self.chord = (stage.tip.radius - stage.root.radius)/self.AR
        self.spacing = (2*m.pi*radius)/self.N
        self.sigma = self.chord/self.spacing
        self.DF = 1 - self.DH\
            + abs(v1*m.sin(m.radians(self.vle))
                  - v2*m.sin(m.radians(self.vte)))/(2*self.sigma*v1)

        # Assumes double circular arc airfoil to estimate stagger
        self.stagger = self.vle*self.xcr + self.vte*(1 - self.xcr)
        self.camber = self.vle - self.vte

        # self.chord = self.spacing*abs(v1*m.sin(m.radians(self.vle))
        #                               - v2*m.sin(m.radians(self.vte)))\
        #     / ((self.DF + self.DH - 1)*(2*v1))

        # Solidworks blade angles
        self.cle = self.vle - self.stagger
        self.cte = self.stagger - self.vte

        # Assume NACA 4 definition of wedge angle and leading edge radius
        self.wte = (2*m.degrees(m.atan(1.16925*(2*self.yt))))
        self.rle = 1.1019*(2*self.yt)**2

        # Calculate natural xc,yc from cle and cte using ycr altitutde ratio
        self.xc = m.tan(m.radians(self.cte))/(m.tan(m.radians(self.cle))
                                              + m.tan(m.radians(self.cte)))
        # self.ycr = (self.xc + 1)/2
        self.yc = self.ycr*self.xc*abs(m.tan(m.radians(self.cle)))

        # Assume NACA 4 curvature
        self.kc, self.kt = u.Curvatures(xc=self.xc, yc=self.yc,
                                        xt=self.xt, yt=self.yt)

        self.bc, self.bt = u.Beziers(xc=self.xc, yc=self.yc, kc=self.kc,
                                     xt=self.xt, yt=self.yt, kt=self.kt,
                                     cle=self.cle, cte=self.cte, rle=self.rle,
                                     blade=blade, station=station)

        k_sh = 1.0
        k_it = (-0.0214
                + 19.17*(2*self.yt)
                - 122.3*(2*self.yt)**2
                + 312.5*(2*self.yt)**3)
        i_010 = ((0.0325 - 0.0674*self.sigma)
                 + (-0.002364 + 0.0913*self.sigma)*self.vle
                 + (.0000164 - 0.000238*self.sigma)*self.vle**2)
        n = ((-0.063 - 0.02274*self.sigma)
             + (-0.0035 + 0.0029*self.sigma)*self.vle
             - (0.0000379 + 0.0000111*self.sigma)**self.vle**2)
        self.incidence = k_sh*k_it*i_010 + n*self.camber
        m_dev = ((0.23*(2*self.xc))**2 + self.vte/500)
        if blade != 'IGV':
            self.deviation = m_dev*self.camber*(1/self.sigma)**0.5
        else:
            self.deviation = m_dev*self.camber*(1/self.sigma)


class Blade():
    """Set blade properties."""

    def __init__(self):
        """Instantiate blade properties."""
        self.root = Airfoil()
        self.mean = Airfoil()
        self.tip = Airfoil()
        self.N = 0

    def GetBladeConfig(self, filename):
        """Read .inp files."""
        self.root.GetAirfoilConfig(filename, 'root')
        self.mean.GetAirfoilConfig(filename, 'mean')
        self.tip.GetAirfoilConfig(filename, 'tip')
        self.N = self.root.N

    def CalcBlade(self, stage, blade):
        """Calculate blade properties."""
        self.root.CalcAirfoil(stage, blade, 'root')
        self.mean.CalcAirfoil(stage, blade, 'mean')
        self.tip.CalcAirfoil(stage, blade, 'tip')

        print('\n')
        print('Mean {blade} Performance'.format(blade=blade))
        print('    DF: {var} (<=0.6)'.format(var=round(self.mean.DF, 2)))
        print('    DH: {var} (>=0.72)'.format(var=round(self.mean.DH, 2)))
        print('    i : {var} deg'.format(var=round(self.mean.incidence, 2)))
        print('    d : {var} deg'.format(var=round(self.mean.deviation, 2)))
