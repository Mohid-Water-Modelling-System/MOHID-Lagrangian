#!/usr/bin/env python3
# Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com
"""Prepare zero-current forcing or analyze an existing buoyancy run."""
from pathlib import Path
import argparse, base64, csv, json, os
import xml.etree.ElementTree as ET
import numpy as np
from netCDF4 import Dataset
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mohid-buoyancy-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CASES = {'settling': 'Buoyancy_settling_Def', 'rising': 'Buoyancy_rising_Def'}
INPUT = HERE/'nc_fields_3D/hydrodynamic/Buoyancy_20240101_0000.nc4'
DURATION = 3600.
CHECKPOINTS = (1800., 3600.)


def prepare_inputs():
    INPUT.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(INPUT, 'w', format='NETCDF4_CLASSIC') as ds:
        for dim, values, units in [
            ('time', np.arange(0., 10801., 3600.), 'seconds since 2024-01-01 00:00:00'),
            ('level', np.linspace(-3000., 0., 13), 'm'),
            ('lat', np.linspace(.9, 1.1, 21), 'degrees_north'),
            ('lon', np.linspace(.9, 1.1, 21), 'degrees_east')]:
            ds.createDimension(dim, len(values))
            v = ds.createVariable(dim, 'f8', (dim,)); v[:] = values; v.units = units
            if dim == 'time': v.calendar = 'standard'
        for key, value, units in [('u', 0., 'm s-1'), ('v', 0., 'm s-1'),
                                  ('w', 0., 'm s-1'), ('bathymetry', -3000., 'm')]:
            v = ds.createVariable(key, 'f8', ('time','level','lat','lon'),
                                  fill_value=9.969209968386869e36, zlib=True)
            v[:] = np.full(v.shape, value); v.units = units
        ds.description = 'All-wet buoyancy verification; zero resolved velocity; bottom at -3000 m.'
    print('Prepared:', INPUT)


def read_vtu(path):
    """Read MOHID inline binary/ascii PointData, across all Pieces.

    Avoid requiring VTK cell connectivity: only particle points/data are needed.
    Supports raw appended data, inline base64, ASCII, and UInt32/UInt64 headers.
    Fails on unsupported compression.
    """
    content=Path(path).read_bytes()
    appended=None
    if b'<AppendedData' in content:
        start=content.index(b'<AppendedData'); tagend=content.index(b'>',start)+1
        end=content.rindex(b'</AppendedData>')
        blob=content[tagend:end].lstrip()
        if not blob.startswith(b'_'): raise ValueError('Missing appended data marker')
        appended=blob[1:]
        if b'encoding="base64"' in content[start:tagend]: appended=base64.b64decode(appended)
        content=content[:tagend]+content[end:]
    root=ET.fromstring(content)
    if root.get('compressor'): raise ValueError('Compressed VTU not supported')
    endian='<' if root.get('byte_order','LittleEndian')=='LittleEndian' else '>'
    header=np.dtype(endian+('u8' if root.get('header_type')=='UInt64' else 'u4'))
    types={'Float64':'f8','Float32':'f4','Int64':'i8','Int32':'i4',
           'UInt64':'u8','UInt32':'u4','UInt8':'u1','Int8':'i1'}
    def array(node):
        dtype=np.dtype(endian+types[node.get('type')]); text=''.join((node.text or '').split())
        if node.get('format')=='ascii': value=np.fromstring(text,sep=' ',dtype=dtype)
        elif node.get('format')=='appended':
            offset=int(node.get('offset')); count=int(np.frombuffer(appended,offset=offset,count=1,dtype=header)[0])
            value=np.frombuffer(appended,offset=offset+header.itemsize,count=count//dtype.itemsize,dtype=dtype)
        elif node.get('format')=='binary':
            # VTK may encode header separately, including base64 padding.
            hchars=4*((header.itemsize+2)//3)
            head=base64.b64decode(text[:hchars])
            nbytes=int(np.frombuffer(head[:header.itemsize],dtype=header)[0])
            raw=base64.b64decode(text)
            if len(raw)<header.itemsize+nbytes:
                raw=base64.b64decode(text[:hchars])+base64.b64decode(text[hchars:])
            value=np.frombuffer(raw[header.itemsize:header.itemsize+nbytes],dtype=dtype)
        else: raise ValueError('Unsupported VTU data encoding')
        nc=int(node.get('NumberOfComponents','1'))
        return value.reshape(-1,nc) if nc>1 else value
    pieces=[]
    for piece in root.findall('.//Piece'):
        if int(piece.get('NumberOfPoints','0'))==0: continue
        data={n.get('Name'):array(n) for n in piece.findall('./PointData/DataArray')}
        data['position']=array(piece.find('./Points/DataArray'))
        pieces.append(data)
    if not pieces: raise ValueError(f'No particle points in {path}')
    data={k:np.concatenate([p[k] for p in pieces]) for k in pieces[0]}
    ids=data['id']; order=np.argsort(ids)
    if len(np.unique(ids))!=len(ids): raise ValueError('Particle IDs are not unique')
    return {k:v[order] for k,v in data.items()}


def specification(work, kind):
    name = CASES[kind]
    """Use the configuration saved with the run, not a subsequently edited source XML."""
    root = ET.parse(work/(name+'.xml')).getroot()
    def constant(key):
        return float(root.find('./caseDefinitions/constants/'+key).get('value'))
    for key in ('UseConstantWaterDensity', 'UseConstantWaterDynamicViscosity'):
        if constant(key) != 1: raise ValueError(key+' must be 1 for constant terminal speed')
    for key in ('DiffusionCoeffHorz','DiffusionCoeffVert','DiffusionCoeffHorzBkg',
                'DiffusionCoeffHorTurbVelLimiter','WindDragCoeff','StokesDriftCoeff','ResuspensionProb'):
        if constant(key) != 0: raise ValueError(key+' must be zero')
    source_type = root.find('.//sourceTypes/types/type')
    if source_type.get('type') != 'paper' or source_type.get('property') != kind+'_sphere':
        raise ValueError('Expected the '+kind+'_sphere paper material')
    material = ET.parse(work/'materialTypes.xml').getroot().find('./paper/'+source_type.get('property'))
    if material is None: raise ValueError('Missing material in the saved materialTypes.xml')
    method = root.find(".//parameter[@key='Integrator']").get('value')
    if int(method) != 1: raise ValueError('This age-zero startup reference is configured for Euler')
    point = root.find('.//sourceDefinitions/source/point')
    initial = np.array([float(point.get(c)) for c in ('x','y','z')])
    rho = constant('MeanDensity'); mu = constant('MeanDynamicViscosity')
    diameter = 2*float(material.find('radius').get('value'))
    particle_rho = float(material.find('density').get('value'))
    gravity = -float(root.find('.//constants/Gravity').get('z'))
    dt = float(root.find('.//timestep').get('dt'))
    if min(rho, mu, diameter, particle_rho, gravity, dt) <= 0: raise ValueError('Nonpositive physical parameter')
    nu = mu/rho
    contrast = abs(1-particle_rho/rho)
    dimensionless_diameter = diameter**3*gravity*contrast/nu**2
    speed = nu/diameter*dimensionless_diameter*(38.1+.93*dimensionless_diameter**(4/7))**(-7/8)
    speed *= 1 if rho >= particle_rho else -1
    if (kind == 'settling' and speed >= 0) or (kind == 'rising' and speed <= 0):
        raise ValueError('Particle density does not give the requested motion direction')
    if not -3000 < initial[2]+speed*DURATION < 0:
        raise ValueError('Analytical trajectory reaches the bottom or surface')
    return dict(motion=kind, rho=rho, mu=mu, nu=nu, diameter=diameter, particle_density=particle_rho,
                gravity=gravity, dt=dt, initial=initial, settling_velocity=speed)


def tex_number(value):
    if value == 0: return '$0$'
    if abs(value) >= 1e-3: return f'${value:.8g}$'
    mantissa, exponent = f'{value:.4e}'.split('e')
    return rf'${mantissa}\times10^{{{int(exponent)}}}$'


def write_table(work, samples, spec):
    nu_tex = tex_number(spec['nu'])[1:-1]
    lines = [r'\begin{table}[htbp]', r'\centering',
        (rf'\caption{{Comparison of analytical and numerical positions for buoyancy-driven {spec["motion"]}. '
         rf'Particle diameter $d={spec["diameter"]*1000:g}~\mathrm{{mm}}$, '
         rf'$\rho_p={spec["particle_density"]:g}~\mathrm{{kg\,m^{{-3}}}}$, '
         rf'$\rho_{{sw}}={spec["rho"]:g}~\mathrm{{kg\,m^{{-3}}}}$, '
         rf'$\nu_{{sw}}={nu_tex}~\mathrm{{m^2\,s^{{-1}}}}$, '
         rf'and $g={spec["gravity"]:g}~\mathrm{{m\,s^{{-2}}}}$. '
         r'All other transport contributions are disabled.}'),
        (r'\label{MOHID-Lagrangian-verrification-buoyancy}' if spec['motion']=='settling'
         else r'\label{MOHID-Lagrangian-verrification-buoyancy-rising}'), r'\footnotesize',
        r'\renewcommand{\arraystretch}{1.15}', r'\resizebox{\textwidth}{!}{%',
        r'\begin{tabular}{c c c c c c c c}', r'\toprule',
        r'Time & Solution & Longitude ($^\circ$) & Latitude ($^\circ$) & $z$ (m) & $Error_{\mathrm{lon}}$ (\%) & $Error_{\mathrm{lat}}$ (\%) & $Error_z$ (\%) \\',
        r'\midrule']
    p=spec['initial']
    lines.append(f'0 h & Initial & {p[0]:.12f} & {p[1]:.12f} & {p[2]:.12f} & & & '+r'\\')
    for t in CHECKPOINTS:
        row=samples[t]; ana=row['analytical']; num=row['position']
        lines += [r'\midrule',
            f'{t/3600:g} h & Analytical & {ana[0]:.12f} & {ana[1]:.12f} & {ana[2]:.12f} & & & '+r'\\',
            f'{t/3600:g} h & Numerical & {num[0]:.12f} & {num[1]:.12f} & {num[2]:.12f} & -- & -- & '+tex_number(row['relative_error_z_percent'])+' '+r'\\']
    lines += [r'\bottomrule',r'\end{tabular}',r'}',r'\par\smallskip\raggedright\footnotesize',
        (r'Time is elapsed buoyancy-driven motion after the initial age-zero resting step. '
         r'$z_{\mathrm{ana}}(t)=z_0+v_3^{\mathrm{buoyancy}}t$; '
         r'$Error_z=100|z_{\mathrm{num}}-z_{\mathrm{ana}}|/|z_{\mathrm{ana}}-z_0|$. '
         r'-- denotes zero analytical horizontal displacement. Errors use unrounded positions.'),
        r'\end{table}']
    (work/'verification_table.tex').write_text('\n'.join(lines)+'\n')


def analyze(work, kind):
    name = CASES[kind]
    if not work.is_dir(): raise ValueError(f'No output folder: {work}. Run the corresponding RunCase_Buoyancy script first.')
    spec=specification(work, kind); rows={}; first_id=None
    files=sorted(work.glob(name+'_[0-9][0-9][0-9][0-9][0-9].vtu'))
    if not files: raise ValueError('No particle VTU files found')
    max_horizontal=0.; max_error=0.; observed_initial=False
    for file in files:
        data=read_vtu(file)
        if len(data['id']) != 1: raise ValueError(f'{file}: expected exactly one particle')
        if first_id is None: first_id=data['id'][0]
        if data['id'][0] != first_id: raise ValueError('Particle identity changed')
        if not all(np.all(np.isfinite(v)) for v in data.values()): raise ValueError('Nonfinite VTU data')
        if np.any(data['state'] != 0): raise ValueError('Particle left the wet water state')
        if np.any(data['velocity'] != 0): raise ValueError('Resolved current is not zero')
        for key,target in [('density',spec['particle_density']),('radius',spec['diameter']/2),
                           ('condition',1.),('degradation_rate',0.),('particulate',1.)]:
            if not np.allclose(data[key],target,rtol=1e-12,atol=1e-14):
                raise ValueError(f'{key} changed or is configured incorrectly')
        age=float(data['age'][0]); pos=data['position'][0]
        # Buoyancy explicitly returns zero at age=0. For Euler the first step
        # advances age, while position remains at the release point.
        t=max(0.,age-spec['dt'])
        if age <= spec['dt']:
            observed_initial=True
            if not np.allclose(pos,spec['initial'],rtol=0,atol=1e-10):
                raise ValueError('Unexpected motion during the initial age-zero resting step')
        if not (.9 < pos[0] < 1.1 and .9 < pos[1] < 1.1 and -3000 < pos[2] <= 0):
            raise ValueError('Particle reached a domain boundary')
        if t > DURATION+1e-8: continue
        ana=spec['initial']+np.array([0.,0.,spec['settling_velocity']*t])
        error=float(abs(pos[2]-ana[2]))
        horizontal=float(np.max(abs(pos[:2]-ana[:2])))
        max_horizontal=max(max_horizontal,horizontal);max_error=max(max_error,error)
        row=dict(time_s=t,age_s=age,position=pos,analytical=ana,error_z_m=error,
                 relative_error_z_percent=100*error/abs(ana[2]-spec['initial'][2]) if t else 0.)
        if t in rows and not np.allclose(pos,rows[t]['position'],rtol=0,atol=1e-10):
            raise ValueError('Position changed without increasing settling time')
        rows[t]=row
    if not observed_initial: raise ValueError('Missing initial snapshot; startup cannot be verified')
    for t in (0., *CHECKPOINTS):
        if t not in rows: raise ValueError(f'Missing exact output at settling time {t:g} s; run to final time with OutputWriteTime=dt')
    records=[rows[t] for t in sorted(rows)]
    with (work/'trajectory.csv').open('w',newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['time_s','age_s','longitude','latitude','z_numerical_m','z_analytical_m','absolute_error_z_m','relative_error_z_percent'])
        for row in records:
            writer.writerow([row['time_s'],row['age_s'],*row['position'],row['analytical'][2],
                             row['error_z_m'],row['relative_error_z_percent']])
    passed=max_horizontal<=1e-10 and max_error<=1e-6
    summary={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in spec.items()}
    summary.update(passed=passed,duration_s=DURATION,max_horizontal_error_deg=max_horizontal,
                   max_vertical_error_m=max_error,absolute_vertical_tolerance_m=1e-6,
                   checkpoints=[{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in rows[t].items()}
                                for t in CHECKPOINTS])
    (work/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    write_table(work,rows,spec)
    with plt.rc_context({'font.family':'serif','font.size':12,'axes.linewidth':1.2,
                         'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42}):
        fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        times=np.array([r['time_s'] for r in records])/3600
        theoretical=np.array([r['analytical'][2] for r in records])
        numerical=np.array([r['position'][2] for r in records])
        axes[0].plot(times,theoretical,'k--',lw=1.8,label='Analytical')
        stride=max(1,len(times)//16)
        axes[0].plot(times,numerical,color='#0072B2',lw=1.4,marker='o',markevery=stride,
                     mfc='white',ms=5,label='Numerical')
        axes[0].set(title='(a) Vertical trajectory',xlabel=spec['motion'].capitalize()+' time (h)',ylabel='$z$ (m)')
        axes[0].legend(frameon=False)
        axes[1].plot(times,np.abs(numerical-theoretical),color='#D55E00',lw=1.8)
        axes[1].set(title='(b) Absolute vertical error',xlabel=spec['motion'].capitalize()+' time (h)',ylabel='Absolute error (m)')
        axes[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
        for ax in axes:
            ax.grid(ls='--',alpha=.35);ax.set_xlim(0,DURATION/3600);ax.set_xticks([0,.25,.5,.75,1]);ax.tick_params(top=True,right=True)
        fig.savefig(work/'Buoyancy_verification.pdf',bbox_inches='tight')
        fig.savefig(work/'Buoyancy_verification.png',dpi=300,bbox_inches='tight');plt.close(fig)
    print(f'{kind}: {"PASS" if passed else "FAIL"}: max vertical error {max_error:.6g} m; horizontal error {max_horizontal:.6g} degrees')
    print('Outputs:',work/'verification_table.tex',work/'Buoyancy_verification.pdf')
    if not passed: raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-inputs',action='store_true',help='Generate the NetCDF only; never run the model')
    parser.add_argument('--case', choices=('both','settling','rising'), default='both')
    parser.add_argument('--input',type=Path,help='Existing output folder; requires --case settling or rising')
    args=parser.parse_args()
    if args.prepare_inputs:
        prepare_inputs()
    else:
        if args.input and args.case=='both': parser.error('--input requires --case settling or rising')
        for kind in (CASES if args.case=='both' else (args.case,)):
            work=args.input.resolve() if args.input else HERE/(CASES[kind]+'_out')
            analyze(work, kind)
