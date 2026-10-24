#!/usr/bin/env python3
"""Prepare zero forcing or verify detritus mass, volume and radius from VTU output."""
from pathlib import Path
import argparse, base64, csv, json, os
import xml.etree.ElementTree as ET
import numpy as np
from netCDF4 import Dataset
os.environ.setdefault('MPLCONFIGDIR','/tmp/mohid-degradation-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
NAME='Degradation_Def'
DURATION=3600.
CHECKPOINTS=(1800.,3600.)
INPUT=HERE/'nc_fields_3D/hydrodynamic/Zero_20240101_0000_[actual_time=20240101].nc4'


def prepare_inputs():
    INPUT.parent.mkdir(parents=True,exist_ok=True)
    with Dataset(INPUT,'w',format='NETCDF4_CLASSIC') as ds:
        for dim,values,units in [('time',np.arange(0.,10801.,3600.),'seconds since 2024-01-01 00:00:00'),
                                ('level',np.linspace(-200,0,5),'m'),('lat',np.linspace(.9,1.1,21),'degrees_north'),
                                ('lon',np.linspace(.9,1.1,21),'degrees_east')]:
            ds.createDimension(dim,len(values));v=ds.createVariable(dim,'f8',(dim,))
            v[:]=values;v.units=units
            if dim=='time':v.calendar='standard'
        for name,value,units in [('u',0.,'m s-1'),('v',0.,'m s-1'),('w',0.,'m s-1'),('bathymetry',-200.,'m')]:
            v=ds.createVariable(name,'f8',('time','level','lat','lon'),fill_value=9.969209968386869e36,zlib=True)
            v[:]=np.full(v.shape,value);v.units=units
        ds.description='Stationary detritus degradation verification: zero resolved velocity'
    print('Prepared:',INPUT)


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
    root=ET.parse(work/(NAME+'.xml')).getroot();constants=root.find('./caseDefinitions/constants')
    def constant(key):return float(constants.find(key).get('value'))
    for key in ('DiffusionCoeffHorz','DiffusionCoeffVert','DiffusionCoeffHorzBkg',
                'DiffusionCoeffHorTurbVelLimiter','WindDragCoeff','StokesDriftCoeff','ResuspensionProb'):
        if constant(key)!=0:raise ValueError(key+' must be zero')
    if constant('UseConstantWaterDensity')!=1:raise ValueError('Constant water density is required')
    kind=root.find('.//sourceTypes/types/type')
    if kind.get('type')!='detritus' or kind.get('property')!='degrading_sphere':raise ValueError('Expected degrading_sphere detritus')
    mat=ET.parse(work/'materialTypes.xml').getroot().find('./detritus/degrading_sphere')
    def prop(key):return float(mat.find(f"property[@key='{key}']").get('value'))
    rho=float(mat.find('density').get('value'));v0=float(mat.find('volume').get('value'))
    if rho!=constant('MeanDensity'):raise ValueError('Neutral buoyancy is required')
    if int(root.find(".//parameter[@key='Integrator']").get('value'))!=1:raise ValueError('Configured reference uses Euler')
    params={key:constant(key) for key in ('ToptBMin','ToptBMax','TBMin','TBMax','BK1','BK2','BK3','BK4',
                                         'Lin0_DegradationRate','Exp0_DegradationRate','Exp1_DegradationRate','Exp2_DegradationRate')}
    # XML rate inputs are per day; globals setters divide these three rates by 86400.
    for key in ('Lin0_DegradationRate','Exp1_DegradationRate','Exp2_DegradationRate'):
        params[key]/=86400.
    s1=np.log(params['BK2']*(1-params['BK1'])/(params['BK1']*(1-params['BK2'])))/(params['ToptBMin']-params['TBMin'])
    s2=np.log(params['BK3']*(1-params['BK4'])/(params['BK4']*(1-params['BK3'])))/(params['TBMax']-params['ToptBMax'])
    ya=np.exp(s1*(prop('temp')-params['TBMin']));yb=np.exp(s2*(params['TBMax']-prop('temp')))
    limit=params['BK1']*ya/(1+params['BK1']*(ya-1))*params['BK4']*yb/(1+params['BK4']*(yb-1))
    point=root.find('.//sourceDefinitions/source/point')
    return dict(density=rho,initial_volume=v0,initial_mass=rho*v0,initial_radius=float(mat.find('radius').get('value')),
                temperature=prop('temp'),salinity=prop('salt'),limiting_factor=float(limit),parameters=params,
                minimum_radius=prop('radius_cr_min'),dt=float(root.find('.//timestep').get('dt')),
                initial_position=np.array([float(point.get(c)) for c in ('x','y','z')]))


def exact(time,spec):
    # The constructor starts at M0. Subsequent updates use M0*L(T)*F(age).
    if time==0:return np.array([spec['initial_mass'],spec['initial_volume'],spec['initial_radius']])
    p=spec['parameters']
    fraction=(1-p['Lin0_DegradationRate']*time)*(p['Exp0_DegradationRate']*np.exp(-p['Exp1_DegradationRate']*time)+
             (1-p['Exp0_DegradationRate'])*np.exp(-p['Exp2_DegradationRate']*time))
    minimum_mass=spec['density']*4*np.pi/3*spec['minimum_radius']**3
    mass=max(spec['initial_mass']*spec['limiting_factor']*fraction,minimum_mass)
    volume=mass/spec['density'];radius=np.cbrt(3*volume/(4*np.pi))
    return np.array([mass,volume,radius])


def tex(value):
    if value==0:return '$0$'
    mantissa,exponent=f'{value:.8e}'.split('e')
    return rf'${mantissa}\times10^{{{int(exponent)}}}$'


def tex_error(value):
    # Display threshold in percent; retain the computed values in CSV and JSON.
    return r'$<10^{-12}$' if value < 1e-12 else tex(value)


def table(work,rows,spec):
    lines=[r'\begin{table}[htbp]',r'\centering',
        (rf'\caption{{Comparison of analytical and numerical particle mass, volume and radius '
         rf'for the detritus degradation verification at constant $T={spec["temperature"]:g}^\circ\mathrm{{C}}$, '
         rf'$\rho={spec["density"]:g}~\mathrm{{kg\,m^{{-3}}}}$, and initial radius '
         rf'$r_0={1000*spec["initial_radius"]:g}~\mathrm{{mm}}$. '
         rf'The temperature limiter is unity to numerical precision. '
         rf'All transport contributions are zero. Euler uses $\Delta t={spec["dt"]:g}~\mathrm{{s}}$.}}'),
        r'\label{tab:verification-degradation}',r'\footnotesize',r'\renewcommand{\arraystretch}{1.15}',
        r'\resizebox{\textwidth}{!}{%',r'\begin{tabular}{c c c c c c c c}',r'\toprule',
        r'Time & Solution & Mass (kg) & Volume ($\mathrm{m^3}$) & Radius (m) & $Error_M$ (\%) & $Error_V$ (\%) & $Error_r$ (\%) \\',
        r'\midrule']
    initial=rows[0.]['numerical']
    lines.append('0 h & Initial & '+' & '.join(tex(v) for v in initial)+r' & & & \\')
    for t in CHECKPOINTS:
        row=rows[t];lines.append(r'\midrule')
        lines.append(f'{t/3600:g} h & Analytical & '+' & '.join(tex(v) for v in row['analytical'])+r' & & & \\')
        lines.append(f'{t/3600:g} h & Numerical & '+' & '.join(tex(v) for v in row['numerical'])+' & '+
                     ' & '.join(tex_error(v) for v in row['relative_error_percent'])+r' \\')
    lines += [r'\bottomrule',r'\end{tabular}',r'}',r'\par\smallskip\raggedright\footnotesize',
        (r'Time is particle age. Numerical mass is reconstructed as $M=\rho V$; density is constant. '
         r'Analytical radius is $r(t)=[3V(t)/(4\pi)]^{1/3}$. '
         r'$Error_q=100|q_{\mathrm{num}}-q_{\mathrm{ana}}|/|q_{\mathrm{ana}}|$, for $q=M,V,r$. '
         r'Errors use unrounded values at the same particle age. '
         r'Errors below $10^{-12}\%$ are displayed as $<10^{-12}\%$, including computed zeros; '
         r'this is a reporting threshold, not an estimate of nonzero error.'),r'\end{table}']
    (work/'verification_table.tex').write_text('\n'.join(lines)+'\n')


def analyze(work):
    if not work.is_dir():raise ValueError('Run RunCase_Degradation.sh first')
    spec=specification(work);rows={};first_id=None
    files=sorted(work.glob(NAME+'_[0-9][0-9][0-9][0-9][0-9].vtu'))
    if not files:raise ValueError('No VTU particle output')
    max_position=0.;rho=spec['density']
    for file in files:
        d=read_vtu(file)
        if len(d['id'])!=1:raise ValueError('Expected one particle')
        if first_id is None:first_id=d['id'][0]
        if d['id'][0]!=first_id:raise ValueError('Particle identity changed')
        if not all(np.all(np.isfinite(v)) for v in d.values()):raise ValueError('Nonfinite output')
        if np.any(d['state']!=0):raise ValueError('Particle is not in wet water')
        for key,target in [('particulate',1.),('density',rho),('initial_volume',spec['initial_volume']),
                           ('temp',spec['temperature']),('salt',spec['salinity']),('radius_cr_min',spec['minimum_radius'])]:
            if not np.allclose(d[key],target,rtol=1e-12,atol=0):raise ValueError('Unexpected change in '+key)
        for key in ('velocity','dVelX','dVelY','dVelZ'):
            if np.any(d[key]!=0):raise ValueError('Nonzero transport in '+key)
        max_position=max(max_position,float(np.max(abs(d['position'][0]-spec['initial_position']))))
        t=float(d['age'][0])
        if t>DURATION+1e-8:continue
        numerical=np.array([rho*float(d['volume'][0]),float(d['volume'][0]),float(d['radius'][0])])
        analytical=exact(t,spec)
        row=dict(time_s=t,numerical=numerical,analytical=analytical,
                 relative_error_percent=100*np.abs(numerical-analytical)/np.abs(analytical))
        if t in rows and not np.array_equal(numerical,rows[t]['numerical']):raise ValueError('Properties changed without aging')
        rows[t]=row
    for t in (0.,*CHECKPOINTS):
        if t not in rows:raise ValueError(f'Missing output at age {t:g}')
    records=[rows[t] for t in sorted(rows)]
    errors=np.array([r['relative_error_percent']/100 for r in records]);max_errors=errors.max(axis=0)
    checks=dict(stationary=max_position<=1e-12,mass=bool(max_errors[0]<=1e-10),
                volume=bool(max_errors[1]<=1e-10),radius_same_time=bool(max_errors[2]<=1e-10))
    passed=all(checks.values())
    def convert(d):return {k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in d.items()}
    summary=convert(spec);summary.update(passed=passed,checks=checks,max_relative_errors=max_errors.tolist(),
        relative_tolerance=1e-10,max_coordinate_drift=max_position,duration_s=DURATION,
        checkpoints=[convert(rows[t]) for t in CHECKPOINTS])
    (work/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (work/'degradation.csv').open('w',newline='') as handle:
        writer=csv.writer(handle);writer.writerow(['time_s','mass_numerical_kg','volume_numerical_m3','radius_numerical_m',
            'mass_analytical_kg','volume_analytical_m3','radius_analytical_m','error_mass_percent','error_volume_percent',
            'error_radius_percent'])
        for row in records:writer.writerow([row['time_s'],*row['numerical'],*row['analytical'],
                                           *row['relative_error_percent']])
    table(work,rows,spec)
    with plt.rc_context({'font.family':'serif','font.size':11,'axes.linewidth':1.2,'pdf.fonttype':42,
                         'xtick.direction':'in','ytick.direction':'in'}):
        fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
        times=np.array([r['time_s'] for r in records])/3600
        num=np.array([r['numerical'] for r in records]);ana=np.array([r['analytical'] for r in records])
        for j,(ax,label,unit) in enumerate(zip(axes.flat,['Mass','Volume','Radius'],['kg',r'$\mathrm{m^3}$','m'])):
            ax.plot(times,ana[:,j],'k--',lw=1.6,label='Analytical, same time')
            ax.plot(times,num[:,j],color='#0072B2',marker='o',markevery=max(1,len(times)//12),mfc='white',ms=4,lw=1,label='Numerical')
            ax.set(title=f'({"abc"[j]}) {label}',xlabel='Time (h)',ylabel=f'{label} ({unit})')
            ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True);ax.legend(frameon=False,fontsize=9)
        ax=axes[1,1]
        for j,label in enumerate(['Mass','Volume','Radius']):ax.plot(times,100*errors[:,j],label=label)
        ax.set(title='(d) Same-time relative errors',xlabel='Time (h)',ylabel='Relative error (%)');ax.legend(frameon=False)
        for ax in axes.flat:ax.grid(ls='--',alpha=.35);ax.set_xlim(0,1);ax.set_xticks([0,.25,.5,.75,1]);ax.tick_params(top=True,right=True)
        fig.savefig(work/'Degradation_verification.pdf',bbox_inches='tight')
        fig.savefig(work/'Degradation_verification.png',dpi=300,bbox_inches='tight');plt.close(fig)
    for key,value in checks.items():print(f'{key}: {"PASS" if value else "FAIL"}')
    print('Reports:',work/'verification_table.tex',work/'summary.json')
    if not passed:raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-inputs',action='store_true')
    parser.add_argument('--input',type=Path,default=HERE/(NAME+'_out'))
    args=parser.parse_args()
    if args.prepare_inputs:prepare_inputs()
    else:analyze(args.input.resolve())
