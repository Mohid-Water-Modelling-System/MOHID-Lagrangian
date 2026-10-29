#!/usr/bin/env python3
# Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com
"""Compile and run small water-property checks against the built Fortran library.

Compile the main project first. No particle cases or production XMLs are run.
The compiler/link settings are taken from the existing CMake build.
"""
from pathlib import Path
import argparse
import re
import shlex
import subprocess
import tempfile
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def run(build):
    link_file = build/'src/MOHIDLagrangianApp/CMakeFiles/MOHIDLagrangian.dir/link.txt'
    flags_file = build/'src/MOHIDLagrangianLib/CMakeFiles/MOHIDLagrangianLib.dir/flags.make'
    if not link_file.exists() or not flags_file.exists():
        raise SystemExit('A configured and compiled CMake Makefiles build is required.')
    link = shlex.split(link_file.read_text())
    compiler = link[0]
    settings = flags_file.read_text()
    flags = shlex.split(re.search(r'^Fortran_FLAGS = (.*)$',settings,re.M)[1])
    includes = shlex.split(re.search(r'^Fortran_INCLUDES = (.*)$',settings,re.M)[1])
    cleaned = []
    skip = False
    for flag in flags:
        if skip:
            skip = False
            continue
        if flag in ['-module','-J']:
            skip = True
        elif not flag.startswith('-J'):
            cleaned.append(flag)
    with tempfile.TemporaryDirectory(prefix='water-properties-') as directory:
        work = Path(directory)
        objects = []
        for name in ['legacy_water_reference','test_water_properties']:
            obj = work/(name+'.o')
            subprocess.run([compiler,*cleaned,*includes,'-I'+str(work),'-I'+str(build/'modules'),
                            '-c',str(HERE/(name+'.f90')),'-o',str(obj)],cwd=work,check=True)
            objects.append(str(obj))
        exe = work/'test_water_properties'
        command = []
        for i,token in enumerate(link):
            if token.endswith('MOHIDLagrangian.f90.o'):
                command.extend(objects)
            elif i>0 and link[i-1]=='-o':
                command.append(str(exe))
            else:
                command.append(token)
        subprocess.run(command,cwd=link_file.parents[2],check=True)
        subprocess.run([str(exe)],cwd=work,check=True)
        invalid = [('density-flag',v) for v in ['-1','2','0.5','NaN','bad']]
        invalid += [('viscosity-flag',v) for v in ['-1','2','0.5','NaN','bad']]
        invalid += [('dynamic-viscosity',v) for v in ['0','-0.001','NaN','Infinity','bad']]
        invalid += [('mean-density',v) for v in ['0','-1','NaN','Infinity','bad']]
        for option,value in invalid:
            result = subprocess.run([str(exe),option,value],cwd=work,capture_output=True,text=True)
            if result.returncode==0:
                raise SystemExit(f'FAIL: accepted {option}={value}')
        print(f'PASS: {len(invalid)} invalid-setting checks')
        # Exercise the real public XML initialization, without running a simulation.
        fixture = REPO/'MOHID-Lagrangian-CodeVerification/CodeVerification-Buoyancy-deriven'
        aliases = ['MeanKVisco', 'MeanKViscosity', 'MeanKinematicViscosity']
        cases = [([], 1.09e-6)]
        cases += [([(name, 2e-6)], 2e-6) for name in aliases]
        cases += [([(aliases[0], 2e-6), (aliases[1], 3e-6)], 3e-6)]
        cases += [([(aliases[2], 4e-6), (aliases[1], 3e-6), (aliases[0], 2e-6)], 4e-6),
                  ([(aliases[0], 2e-6), (aliases[1], 3e-6), (aliases[2], 4e-6)], 4e-6)]
        for index, (entries, expected) in enumerate(cases):
            tree = ET.parse(fixture/'Buoyancy_settling_Def.xml')
            root = tree.getroot()
            constants = root.find('./caseDefinitions/constants')
            for name, value in entries:
                ET.SubElement(constants, name, value=str(value))
            for node in root.findall('.//file'):
                node.set('name', str(fixture/node.get('name')))
            config = work/f'kinematic_alias_{index}.xml'
            tree.write(config, encoding='utf-8', xml_declaration=True)
            result = subprocess.run([str(exe), 'kinematic-xml', str(config), str(expected)],
                                    cwd=work, capture_output=True, text=True)
            if result.returncode != 0:
                raise SystemExit(result.stdout+result.stderr)
        print(f'PASS: {len(cases)} kinematic-viscosity XML compatibility checks')



if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,default=REPO/'build')
    args=parser.parse_args()
    run(args.build.resolve())
