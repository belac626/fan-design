"""Writes geometry vars to Solidworks equation files."""
import os

import numpy as np


################################
# #Class: Solidworks
# Holds Solidworks equation file data
# #Attributes:
# cwd: current working directory (python files)
# root: project directory
# dir: directory of filename
# name: filename that Solidworks reads
################################
class Solidworks():
    """Set Solidworks global variables."""

    def __init__(self):
        """Instantiate file properties."""
        self.cwd = os.getcwd()
        os.chdir('..')
        self.root = os.getcwd()
        os.chdir(r'.\Solidworks\Vars')
        self.dir = os.getcwd()
        os.chdir(self.cwd)
        self.name = ''

    def WriteConfig(self, stage, blade, filename):
        """Write Solidworks equation file."""
        self.name = filename
        os.chdir(self.dir)

        keys = ['Blades', 'ID', 'OD']
        stations = ['root', 'mean', 'tip']
        for s in stations:
            keys.extend(['Chord_' + s, 'Stagger_' + s,
                         'xt_' + s, 'yt_' + s, 'kt_' + s, 'bt_' + s,
                         'xc_' + s, 'yc_' + s, 'kc_' + s, 'bc_' + s,
                         'cle_' + s, 'cte_' + s, 'rle_' + s, 'wte_' + s])

        vars = [blade.z, stage.root.radius*2, stage.tip.radius*2,
                blade.root.chord, blade.root.stagger,
                blade.root.xt, blade.root.yt, blade.root.kt, blade.root.bt,
                blade.root.xc, blade.root.yc, blade.root.kc, blade.root.bc,
                blade.root.cle, blade.root.cte, blade.root.rle, blade.root.wte,
                blade.mean.chord, blade.mean.stagger,
                blade.mean.xt, blade.mean.yt, blade.mean.kt, blade.mean.bt,
                blade.mean.xc, blade.mean.yc, blade.mean.kc, blade.mean.bc,
                blade.mean.cle, blade.mean.cte, blade.mean.rle, blade.mean.wte,
                blade.tip.chord, blade.tip.stagger,
                blade.tip.xt, blade.tip.yt, blade.tip.kt, blade.tip.bt,
                blade.tip.xc, blade.tip.yc, blade.tip.kc, blade.tip.bc,
                blade.tip.cle, blade.tip.cte, blade.tip.rle, blade.tip.wte]
        blade_dict = dict(zip(keys, np.round(vars, 4)))

        with open(self.name, 'w') as SWConfig:
            output = []
            for key, value in blade_dict.items():
                line = f'"{key}"= {value}'  # .format(key=key, value=value)
                output.append(line)
            output = '\n'.join(output)
            SWConfig.write(output)

        os.chdir(self.cwd)
