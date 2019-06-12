"""Used to implement variables in Solidworks Leaf Shredder."""
import os
import pathlib as pth

import utils as ut
import write_file as wf
from fan_flow import Stage
from fan_geo import Blade

CWD = os.getcwd()

STAGE = Stage()
STAGE.get_stage_config(filename='Stage.ini')
STAGE.calcstage()

IGV = Blade()
IGV.get_blade_config(filename='IGV.ini')
IGV.mean.rle = 1.1019*(2*IGV.mean.yt)**2
IGV.mean.kt = ut.thicknesscurvature(IGV.mean.xt, IGV.mean.yt)
IGV.mean.bt = ut.thicknessbezier(IGV.mean.xt, IGV.mean.yt,
                                 IGV.mean.kt, IGV.mean.rle,
                                 blade='IGV', station='mean')
ROTOR = Blade()
ROTOR.get_blade_config(filename='Rotor.ini')
ROTOR.calcblade(stage=STAGE, blade='Rotor')
STATOR = Blade()
STATOR.get_blade_config(filename='Stator.ini')
STATOR.calcblade(stage=STAGE, blade='Stator')

POD_DICT = {
    'Pod chord': 12,
    'Rotor chord': ROTOR.root.chord/0.0254,
    'Stator chord': STATOR.root.chord/0.0254,
    'ID': STAGE.root.radius*2/0.0254,
    'OD': STAGE.tip.radius*2/0.0254,
    'Shaft OD Nom.': 0.5,
    'Shaft OD': '"Shaft OD Nom." + 0.1',
    'Wall Thickness': 0.25,
    'Pre Blades': 3,
    'Blades': 3,
    'Tip Clearance': 0.05,
    # 'Chord': IGV.mean.chord,
    'xt': IGV.mean.xt,
    'yt': IGV.mean.yt,
    'kt': IGV.mean.kt,
    # 'xc': IGV.mean.xc,
    # 'yc': IGV.mean.yc,
    # 'kc': IGV.mean.kc,
    'rle': IGV.mean.rle,
    # 'cle': abs(IGV.mean.cle),
    # 'cte': abs(IGV.mean.cte),
    # 'wte': igb.mean.wte,
    'bt': IGV.mean.bt,
    # 'bc': IGV.mean.bc,
    # 'Stagger': IGV.mean.stagger
}

os.chdir('..')
POD_FILE = pth.Path('Solidworks/Vars/Pod.txt')
with open(POD_FILE, 'w') as pod:
    OUTPUT = []
    for key, value in POD_DICT.items():
        line = f'"{key}"= {value}'
        OUTPUT.append(line)
    OUTPUT = '\n'.join(OUTPUT)
    pod.write(OUTPUT)
os.chdir(CWD)

wf.write_sw_config(stage=STAGE, blade=ROTOR, filename='Rotor.txt')

wf.write_sw_config(stage=STAGE, blade=STATOR, filename='Stator.txt')

wf.DumpVars.dump_stage_csv(filename='Stage_dump.csv', stage=STAGE)
wf.DumpVars.dump_blade_csv(filename='Blade_dump.csv',
                           igv=IGV, rotor=ROTOR, stator=STATOR)
wf.DumpVars.dump_json(filename='Vars_dump.json',
                      stage=STAGE, igv=IGV, rotor=ROTOR, stator=STATOR)
