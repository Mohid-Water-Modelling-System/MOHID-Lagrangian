#!/usr/bin/env python3
# Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com
"""Prepare constant forcing or analyze an existing coupled-process run."""
from pathlib import Path
import argparse, base64, csv, json, os
import xml.etree.ElementTree as ET
import numpy as np
from netCDF4 import Dataset
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mohid-coupledprocesses-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
NAME = 'CoupledProcesses_Def'
DURATION = 3600.
CHECKPOINTS = (1800., 3600.)
RADIUS = 6378137.
STOKES = np.array([0.2, 0.4])
CURRENT = np.array([0.2, 0.3, 0.])
WIND = np.array([10., 15.])
HYDRO = HERE/'nc_fields_3D/hydrodynamic/Currents_20240101_0000.nc4'
WAVES = HERE/'nc_fields_3D/waves/Stokes_20240101_0000.nc4'

METEO = HERE/'nc_fields_3D/winds/Wind_20240101_0000.nc4'

def prepare_inputs():
    for path, kind in [(HYDRO,'hydro'), (METEO,'wind'), (WAVES,'waves')]:
        path.parent.mkdir(parents=True, exist_ok=True)
        with Dataset(path,'w',format='NETCDF4_CLASSIC') as ds:
            coords=[('time',np.arange(0.,10801.,3600.),'seconds since 2024-01-01 00:00:00')]
            if kind=='hydro':coords += [('level',np.linspace(-200,0,5),'m')]
            coords += [('lat',np.linspace(.9,1.1,21),'degrees_north'),('lon',np.linspace(.9,1.1,21),'degrees_east')]
            for dim,values,units in coords:
                ds.createDimension(dim,len(values));v=ds.createVariable(dim,'f8',(dim,))
                v[:]=values;v.units=units
                if dim=='time':v.calendar='standard'
            dimensions=('time','level','lat','lon') if kind=='hydro' else ('time','lat','lon')
            fields={'hydro':[('u',.2,'m s-1'),('v',.3,'m s-1'),('w',0.,'m s-1'),('bathymetry',-200.,'m')],
                    'wind':[('u10',10.,'m s-1'),('v10',15.,'m s-1'),('rad',0.,'W m-2')],
                    'waves':[('VSDX',.2,'m s-1'),('VSDY',.4,'m s-1'),('hs',0.,'m'),('ts',10.,'s'),
                             ('wd',0.,'degrees'),('wl',100.,'m')]}[kind]
            for name,value,units in fields:
                v=ds.createVariable(name,'f8',dimensions,fill_value=9.969209968386869e36,zlib=True)
                v[:]=np.full(v.shape,value);v.units=units
            ds.description='Coupled-process verification: spatially and temporally constant '+kind+' fields'
        print('Prepared:',path)


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


def specification(work):
    root=ET.parse(work/(NAME+'.xml')).getroot()
    constants=root.find('./caseDefinitions/constants')
    def value(key): return float(constants.find(key).get('value'))
    for key in ('DiffusionCoeffHorz','DiffusionCoeffVert','DiffusionCoeffHorzBkg',
                'DiffusionCoeffHorTurbVelLimiter','ResuspensionProb'):
        if value(key)!=0:raise ValueError(key+' must be zero')
    coeff=value('StokesDriftCoeff');wind_coeff=value('WindDragCoeff')
    if coeff!=0.03 or wind_coeff!=0.02:raise ValueError('Expected wave coefficient 0.03 and wind coefficient 0.02')
    for key in ('UseConstantWaterDensity','UseConstantWaterDynamicViscosity'):
        if value(key)!=1:raise ValueError(key+' must be 1')
    kind=root.find('.//sourceTypes/types/type')
    if kind.get('type')!='paper' or kind.get('property')!='neutral_sphere':
        raise ValueError('Expected the neutrally buoyant paper material')
    material=ET.parse(work/'materialTypes.xml').getroot().find('./paper/neutral_sphere')
    rho_p=float(material.find('density').get('value'));rho_water=value('MeanDensity')
    if rho_p!=rho_water:raise ValueError('Particle and water density must match to preserve z=0')
    if int(root.find(".//parameter[@key='Integrator']").get('value'))!=1:
        raise ValueError('This verification is configured for Euler')
    point=root.find('.//sourceDefinitions/source/point')
    initial=np.array([float(point.get(c)) for c in ('x','y','z')])
    if not np.array_equal(initial,[1.,1.,0.]): raise ValueError('Expected release at (1 degree,1 degree,0 m)')
    return dict(initial=initial,stokes_coefficient=coeff,windage_coefficient=wind_coeff,
                dt=float(root.find('.//timestep').get('dt')),current=CURRENT.copy(),wind_10m=WIND.copy(),
                prescribed_stokes_velocity=STOKES.copy(),stokes_velocity=coeff*STOKES.copy(),
                windage_velocity=wind_coeff*WIND.copy(),total_velocity=CURRENT[:2]+wind_coeff*WIND+coeff*STOKES,
                particle_density=rho_p,water_density=rho_water,
                particle_radius=float(material.find('radius').get('value')),earth_radius_m=RADIUS)


def exact(time, spec):
    """Closed-form spherical trajectory with changing meridional position."""
    ld=np.longdouble
    pi=np.arccos(ld(-1))
    initial=spec['initial'];east,north=map(ld,spec['total_velocity'])
    phi0=ld(initial[1])*pi/180
    phi=phi0+north*ld(time)/ld(RADIUS)
    if north==0:
        longitude=ld(initial[0])+east*ld(time)/(ld(RADIUS)*np.cos(phi0))*180/pi
    else:
        longitude=ld(initial[0])+east/north*180/pi*(np.arctanh(np.sin(phi))-np.arctanh(np.sin(phi0)))
    return np.array([float(longitude),float(phi*180/pi),initial[2]])


def tex_number(value):
    if value==0: return '$0$'
    mantissa,exponent=f'{value:.4e}'.split('e')
    return rf'${mantissa}\times10^{{{int(exponent)}}}$'


def write_table(work,rows,spec):
    lines=[r'\begin{table}[htbp]',r'\centering',
        (r'\caption{Comparison of analytical and numerical particle positions for coupled '
         r'advection, windage and Stokes drift. Currents $(u,v)=(0.2,0.3)~\mathrm{m\,s^{-1}}$, '
         r'wind $(u_{10},v_{10})=(10,15)~\mathrm{m\,s^{-1}}$ with $C^{\mathrm{wind}}=0.02$, '
         r'and surface Stokes velocity $(0.2,0.4)~\mathrm{m\,s^{-1}}$ with $C^{\mathrm{wave}}=0.03$. '
         r'The initial position is $(1^\circ,1^\circ,0~\mathrm{m})$. '
         r'The paper tracer is neutrally buoyant; diffusion and resuspension are disabled.}'),
        r'\label{MOHID-Lagrangian-verrification-all}', r'\footnotesize',
        r'\renewcommand{\arraystretch}{1.15}',r'\resizebox{\textwidth}{!}{%',
        r'\begin{tabular}{c c c c c c c c}',r'\toprule',
        r'Time & Solution & Longitude ($^\circ$) & Latitude ($^\circ$) & $z$ (m) & $Error_{\mathrm{lon}}$ (\%) & $Error_{\mathrm{lat}}$ (\%) & $Error_z$ (\%) \\',r'\midrule']
    x=spec['initial']
    lines.append(f'0 h & Initial & {x[0]:.13f} & {x[1]:.13f} & {x[2]:.10f} & & & '+r'\\')
    for t in CHECKPOINTS:
        row=rows[t];a=row['analytical'];n=row['position'];rel=row['relative_error_percent']
        lines += [r'\midrule',
            f'{t/3600:g} h & Analytical & {a[0]:.13f} & {a[1]:.13f} & {a[2]:.10f} & & & '+r'\\',
            f'{t/3600:g} h & Numerical & {n[0]:.13f} & {n[1]:.13f} & {n[2]:.10f} & '+
            tex_number(rel[0])+' & '+tex_number(rel[1])+' & -- '+r'\\']
    lines += [r'\bottomrule',r'\end{tabular}',r'}',r'\par\smallskip\raggedright\footnotesize',
        (r'Time is elapsed particle age. The analytical longitude accounts for the '
         r'instantaneous latitude in the spherical metric. '
         r'$Error_j=100|x_j^{\mathrm{num}}-x_j^{\mathrm{ana}}|/|x_j^{\mathrm{ana}}-x_j(0)|$. '
         r'-- denotes zero analytical vertical displacement. Errors use unrounded positions.'),r'\end{table}']
    (work/'verification_table.tex').write_text('\n'.join(lines)+'\n')


def analyze(work):
    if not work.is_dir(): raise ValueError('Run RunCase_CoupledProcesses.sh first: '+str(work))
    spec=specification(work);rows={};first_id=None
    files=sorted(work.glob(NAME+'_[0-9][0-9][0-9][0-9][0-9].vtu'))
    if not files: raise ValueError('No particle VTU files found')
    for file in files:
        data=read_vtu(file)
        if len(data['id'])!=1: raise ValueError('Expected one tracer')
        if first_id is None:first_id=data['id'][0]
        if data['id'][0]!=first_id:raise ValueError('Particle identity changed')
        if not all(np.all(np.isfinite(v)) for v in data.values()):raise ValueError('Nonfinite VTU data')
        if np.any(data['state']!=0):raise ValueError('Particle left the wet water state')
        if np.any(data['particulate']!=1):raise ValueError('Expected a particulate paper tracer')
        for key,target in [('density',spec['particle_density']),('radius',spec['particle_radius']),
                           ('condition',1.),('degradation_rate',0.)]:
            if not np.allclose(data[key],target,rtol=0,atol=1e-12):raise ValueError('Material property changed: '+key)
        # The age-zero VTU is written before the first current interpolation.
        if float(data['age'][0])>0 and not np.allclose(data['velocity'],spec['current'],rtol=0,atol=1e-12):
            raise ValueError('Resolved current differs from prescribed input')
        for key in ('dVelX','dVelY','dVelZ'):
            if np.any(data[key]!=0):raise ValueError('Unexpected stochastic velocity in '+key)
        t=float(data['age'][0]);position=data['position'][0]
        if not (.9<position[0]<1.1 and .9<position[1]<1.1 and -200<position[2]<=0):
            raise ValueError('Particle reached the domain boundary')
        if t>DURATION+1e-8:continue
        analytical=exact(t,spec);delta=position-analytical
        metric=RADIUS*np.pi/180
        horizontal_error=float(np.hypot(delta[0]*metric*np.cos(np.deg2rad(analytical[1])),delta[1]*metric))
        displacement=analytical-spec['initial']
        relative=np.zeros(2)
        if t>0:relative=100*np.abs(delta[:2])/np.abs(displacement[:2])
        row=dict(time_s=t,position=position,analytical=analytical,horizontal_error_m=horizontal_error,
                 vertical_error_m=float(abs(delta[2])),relative_error_percent=relative)
        if t in rows and not np.array_equal(position,rows[t]['position']):
            raise ValueError('Motion without increasing age')
        rows[t]=row
    for t in (0.,*CHECKPOINTS):
        if t not in rows:raise ValueError(f'Missing exact output at particle age {t:g} s')
    if not np.allclose(rows[0.]['position'],spec['initial'],rtol=0,atol=1e-12):
        raise ValueError('Incorrect initial position')
    records=[rows[t] for t in sorted(rows)]
    maximum_horizontal=max(r['horizontal_error_m'] for r in records)
    maximum_vertical=max(r['vertical_error_m'] for r in records)
    maximum_relative=max(float(np.max(rows[t]['relative_error_percent'])) for t in CHECKPOINTS)
    passed=maximum_horizontal<=1e-3 and maximum_vertical<=1e-10 and maximum_relative<=1e-5
    def serializable(row):return {k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in row.items()}
    summary=serializable(spec)
    summary.update(passed=passed,duration_s=DURATION,integrator='Euler',max_horizontal_error_m=maximum_horizontal,
                   max_vertical_error_m=maximum_vertical,max_checkpoint_relative_error_percent=maximum_relative,
                   horizontal_tolerance_m=1e-3,vertical_tolerance_m=1e-10,
                   checkpoint_relative_tolerance_percent=1e-5,
                   checkpoints=[serializable(rows[t]) for t in CHECKPOINTS])
    (work/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (work/'trajectory.csv').open('w',newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['time_s','lon_numerical_deg','lat_numerical_deg','z_numerical_m',
                         'lon_analytical_deg','lat_analytical_deg','z_analytical_m',
                         'horizontal_error_m','vertical_error_m','relative_error_lon_percent','relative_error_lat_percent'])
        for row in records:writer.writerow([row['time_s'],*row['position'],*row['analytical'],
                                            row['horizontal_error_m'],row['vertical_error_m'],*row['relative_error_percent']])
    write_table(work,rows,spec)
    with plt.rc_context({'font.family':'serif','font.size':12,'axes.linewidth':1.2,
                         'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42}):
        fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        times=np.array([r['time_s'] for r in records])/3600
        analytical=np.array([r['analytical'] for r in records]);numerical=np.array([r['position'] for r in records])
        scale=RADIUS*np.pi/180;phi0=np.deg2rad(spec['initial'][1])
        a=(analytical-spec['initial'])*[scale*np.cos(phi0),scale,1]
        n=(numerical-spec['initial'])*[scale*np.cos(phi0),scale,1]
        for j,(color,marker,label) in enumerate([('#0072B2','o','East'),('#D55E00','s','North')]):
            axes[0].plot(times,a[:,j],color=color,lw=1.6,label=label+' analytical')
            axes[0].plot(times,n[:,j],ls='none',color=color,marker=marker,markevery=max(1,len(times)//12),
                         mfc='white',ms=5,label=label+' numerical')
        axes[0].set(title='(a) Horizontal displacement',xlabel='Time (h)',ylabel='Local displacement (m)')
        axes[0].legend(frameon=False,fontsize=10)
        axes[1].plot(times,[r['horizontal_error_m'] for r in records],color='#009E73',lw=1.8)
        axes[1].set(title='(b) Horizontal position error',xlabel='Time (h)',ylabel='Horizontal error (m)')
        axes[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
        for ax in axes:ax.grid(ls='--',alpha=.35);ax.set_xlim(0,1);ax.set_xticks([0,.25,.5,.75,1]);ax.tick_params(top=True,right=True)
        fig.savefig(work/'CoupledProcesses_verification.pdf',bbox_inches='tight')
        fig.savefig(work/'CoupledProcesses_verification.png',dpi=300,bbox_inches='tight');plt.close(fig)
    print(f'{"PASS" if passed else "FAIL"}: maximum horizontal error {maximum_horizontal:.6g} m; vertical error {maximum_vertical:.6g} m')
    print('Outputs:',work/'verification_table.tex',work/'CoupledProcesses_verification.pdf')
    if not passed:raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-inputs',action='store_true',help='Create NetCDF inputs only; never run the model')
    parser.add_argument('--input',type=Path,default=HERE/(NAME+'_out'),help='Existing output directory')
    args=parser.parse_args()
    if args.prepare_inputs:prepare_inputs()
    else:analyze(args.input.resolve())
