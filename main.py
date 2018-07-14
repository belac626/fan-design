"""Used to implement variables in Solidworks Leaf Shredder."""
import json
import os

import utils as u
from write_sw import Solidworks

compressor = False

if not compressor:
    from fan_flow import Stage
    from fan_geo import Blade
else:
    from comp_thermoflow import Stage
    from comp_geo import Blade

cwd = os.getcwd()
if compressor is True:
    stage = Stage()
    stage.GetStageConfig(filename='Comp_Stage.ini')
    stage.CalcStage()

    igv = Blade()
    igv.GetBladeConfig(filename='Comp_IGV.ini')
    igv.CalcBlade(stage=stage, blade='IGV')
    rotor = Blade()
    rotor.GetBladeConfig(filename='Comp_Rotor.ini')
    rotor.CalcBlade(stage=stage, blade='Rotor')
    stator = Blade()
    stator.GetBladeConfig(filename='Comp_Stator.ini')
    stator.CalcBlade(stage=stage, blade='Stator')

    pod_dict = {
        'Pod chord': 12,
        'Rotor chord': rotor.root.chord,
        'Stator chord': stator.root.chord,
        'ID': stage.root.radius*2,
        'OD': stage.tip.radius*2,
        'Shaft OD Nom.': 0.5,
        'Shaft OD': '"Shaft OD Nom." + 0.1',
        'Wall Thickness': 0.25,
        'Pre Blades': 3,
        'Blades': 3,
        'Tip Clearance': 0.05,
        'Chord': igv.mean.chord,
        'xt': igv.mean.xt,
        'yt': igv.mean.yt,
        'kt': igv.mean.kt,
        'xc': igv.mean.xc,
        'yc': igv.mean.yc,
        'kc': igv.mean.kc,
        'rle': igv.mean.rle,
        'cle': abs(igv.mean.cle),
        'cte': abs(igv.mean.cte),
        'wte': igv.mean.wte,
        'bt': igv.mean.bt,
        'bc': igv.mean.bc,
        'Stagger': igv.mean.stagger
    }
else:
    stage = Stage()
    stage.GetStageConfig(filename='Fan_Stage.ini')
    stage.CalcStage()

    igv = Blade()
    igv.GetBladeConfig(filename='Fan_IGV.ini')
    igv.mean.rle = 1.1019*(2*igv.mean.yt)**2
    igv.mean.kt = u.ThicknessCurvature(igv.mean.xt, igv.mean.yt)
    igv.mean.bt = u.ThicknessBezier(igv.mean.xt, igv.mean.yt,
                                    igv.mean.kt, igv.mean.rle,
                                    blade='IGV', station='mean')
    rotor = Blade()
    rotor.GetBladeConfig(filename='Fan_Rotor.ini')
    rotor.CalcBlade(stage=stage, blade='Rotor')
    stator = Blade()
    stator.GetBladeConfig(filename='Fan_OSV.ini')
    stator.CalcBlade(stage=stage, blade='OSV')

    pod_dict = {
        'Pod chord': 12,
        'Rotor chord': rotor.root.chord,
        'Stator chord': stator.root.chord,
        'ID': stage.root.radius*2,
        'OD': stage.tip.radius*2,
        'Shaft OD Nom.': 0.5,
        'Shaft OD': '"Shaft OD Nom." + 0.1',
        'Wall Thickness': 0.25,
        'Pre Blades': 3,
        'Blades': 3,
        'Tip Clearance': 0.05,
        # 'Chord': igv.mean.chord,
        'xt': igv.mean.xt,
        'yt': igv.mean.yt,
        'kt': igv.mean.kt,
        # 'xc': igv.mean.xc,
        # 'yc': igv.mean.yc,
        # 'kc': igv.mean.kc,
        'rle': igv.mean.rle,
        # 'cle': abs(igv.mean.cle),
        # 'cte': abs(igv.mean.cte),
        # 'wte': igb.mean.wte,
        'bt': igv.mean.bt,
        # 'bc': igv.mean.bc,
        # 'Stagger': igv.mean.stagger
    }

os.chdir('..')
os.chdir(r'.\Solidworks\Vars')
with open('Pod.txt', 'w') as pod:
    output = []
    for key, value in pod_dict.items():
        line = '"{key}"= {value}'.format(key=key, value=value)
        output.append(line)
    output = '\n'.join(output)
    pod.write(output)
os.chdir(cwd)

rotor_file = Solidworks()
rotor_file.WriteConfig(stage=stage, blade=rotor,
                       filename='Rotor.txt')

stator_file = Solidworks()
stator_file.WriteConfig(stage=stage, blade=stator,
                        filename='Stator.txt')

with open('vars_dump.json', 'w') as dump:
    json.dump(json.loads(json.dumps([[stage.root.__dict__,
                                     stage.mean.__dict__,
                                     stage.tip.__dict__],
                                     [igv.root.__dict__,
                                      igv.mean.__dict__,
                                      igv.tip.__dict__],
                                     [rotor.root.__dict__,
                                     rotor.mean.__dict__,
                                     rotor.tip.__dict__],
                                     [stator.root.__dict__,
                                     stator.mean.__dict__,
                                     stator.tip.__dict__]]),
                         parse_float=lambda x: round(float(x), 3)),
              dump, indent=2)
