"""Used to implement variables in Solidworks Leaf Shredder."""
import os

import utils as u
import write_file as wf

compressor = False

if compressor:
    from comp_thermoflow import Stage
    from comp_geo import Blade
else:
    from fan_flow import Stage
    from fan_geo import Blade

cwd = os.getcwd()
if compressor:
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
        'Rotor chord': rotor.root.chord*0.0254,
        'Stator chord': stator.root.chord*0.0254,
        'ID': stage.root.radius*2*0.0254,
        'OD': stage.tip.radius*2*0.0254,
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
        line = f'"{key}"= {value}'
        output.append(line)
    output = '\n'.join(output)
    pod.write(output)
os.chdir(cwd)

rotor_file = wf.Solidworks()
rotor_file.WriteConfig(stage=stage, blade=rotor, filename='Rotor.txt')

stator_file = wf.Solidworks()
stator_file.WriteConfig(stage=stage, blade=stator, filename='Stator.txt')

wf.DumpVars.Dump_Blade_Csv(filename=r'.\Output\Blade_dump.csv',
                           rotor=rotor, stator=stator)
wf.DumpVars.Dump_Stage_Csv(filename=r'.\Output\Stage_dump.csv',
                           stage=stage)
wf.DumpVars.Dump_Json(filename=r'.\Output\Vars_dump.json',
                      stage=stage, igv=igv, rotor=rotor, stator=stator)
