#!/usr/bin/env python3
# Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com
"""Prepare shared NetCDF forcing or analyze completed advection runs."""
from pathlib import Path
import argparse, base64, csv, json, os, re
from datetime import datetime, timedelta
import xml.etree.ElementTree as ET
import numpy as np
from netCDF4 import Dataset
os.environ.setdefault('MPLCONFIGDIR','/tmp/mohid-advection-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

HERE=Path(__file__).resolve().parent
# The project compiles default reals at 64 bits (-r8 / -fdefault-real-8).
MODEL_PI=float(np.pi)
DEGREE=6378137.*MODEL_PI/180.
A=5e-5; B=2.5e-5; GAMMA=1e-3; T=3600.; ORIGIN_TIME=240.
INITIAL=np.array([100.,100.,-100.])
METHODS={'Eluer1':1,'Eluer2':2,'RK4':3}
STEPS=(120,60,30,15)
INPUT=HERE/'nc_fields_3D/hydrodynamic/Adv_20240101_0000.nc4'

def exact(t):
    return INITIAL*np.exp(np.asarray(t)[...,None]*(1+GAMMA*np.asarray(t)[...,None]/2)*np.array([A,B,-A-B]))

def prepare_inputs(refinement=1):
    """Same forcing for all dt: linear in time/X/Z; latitude resolves metric factor."""
    INPUT.parent.mkdir(parents=True,exist_ok=True)
    x=np.linspace(-1000,1000,21); y=np.linspace(-1000,1000,200*refinement+1)
    z=np.linspace(-400,0,9); times=np.arange(0.,4081.,60.)
    with Dataset(INPUT,'w',format='NETCDF4_CLASSIC') as ds:
        for dim,values,unit in [('time',times,'seconds since 2024-01-01 00:00:00'),
                                ('level',z,'m'),('lat',y/DEGREE,'degrees_north'),
                                ('lon',x/DEGREE,'degrees_east')]:
            ds.createDimension(dim,len(values));v=ds.createVariable(dim,'f8',(dim,))
            v[:]=values;v.units=unit
            if dim=='time': v.calendar='standard'
        fields={}
        for key,unit in [('u','m s-1'),('v','m s-1'),('w','m s-1'),('bathymetry','m')]:
            fields[key]=ds.createVariable(key,'f8',('time','level','lat','lon'),
                    fill_value=9.969209968386869e36,zlib=True)
            fields[key].units=unit
        shape=(len(z),len(y),len(x))
        for i,time in enumerate(times):
            f=1+GAMMA*(time-ORIGIN_TIME)
            fields['u'][i]=np.broadcast_to(A*x[None,None,:]*np.cos(y[None,:,None]/6378137.)*f,shape)
            fields['v'][i]=np.broadcast_to(B*y[None,:,None]*f,shape)
            fields['w'][i]=np.broadcast_to(-(A+B)*z[:,None,None]*f,shape)
            fields['bathymetry'][i]=-400.
        ds.description='Advection only. X=degree*lon, Y=degree*lat, Z=level; physical t=seconds since midnight minus 240. u includes cos(latitude) to cancel model metric conversion.'
        ds.model_pi=MODEL_PI;ds.latitude_refinement=refinement
    print(INPUT)

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

def case_settings(path):
    """The case filename defines the method and time step."""
    match=re.fullmatch(r'Advetion_Verification_(Eluer1|Eluer2|RK4)_(015|030|060|120)',path.stem)
    if not match:
        raise ValueError(f'{path.name}: use a standard method/time-step filename')
    method,step=match.groups();dt=int(step)
    epoch=datetime(2024,1,1)
    return method,dt,{
        ('parameter','key','Start','value'):(epoch+timedelta(seconds=ORIGIN_TIME-dt)).strftime('%Y %m %d %H %M %S'),
        ('parameter','key','End','value'):(epoch+timedelta(seconds=ORIGIN_TIME+T+dt)).strftime('%Y %m %d %H %M %S'),
        ('parameter','key','BaseDateTime','value'):epoch.strftime('%Y %m %d %H %M %S'),
        ('parameter','key','Integrator','value'):str(METHODS[method]),
        ('parameter','key','OutputWriteTime','value'):str(dt),
        ('timestep',None,None,'dt'):str(dt),
        ('rate_seconds',None,None,'value'):str(dt),
        ('rate_trcPerEmission',None,None,'value'):'0.5',
        ('active',None,None,'start'):'0',
        ('active',None,None,'end'):f'{1.5*dt:g}',
    }


def check_case(path):
    method,dt,settings=case_settings(path)
    root=ET.parse(path).getroot()
    for (tag,key,value,attr),expected in settings.items():
        query=f".//{tag}[@{key}='{value}']" if key else './/'+tag
        nodes=root.findall(query)
        if len(nodes)!=1:
            raise ValueError(f'{path}: expected exactly one {query}')
        actual=nodes[0].get(attr)
        equal=(actual==expected if tag=='parameter' and value in ['Start','End','BaseDateTime']
               else actual is not None and float(actual)==float(expected))
        if not equal:
            raise ValueError(f'{path}: {value or tag}.{attr}={actual}; expected {expected}. '
                             'Run python analyze_results.py --sync-cases, then rerun this case. '
                             'Existing outputs must be regenerated, not relabelled.')
    for key in ['DiffusionCoeffHorz','DiffusionCoeffVert','DiffusionCoeffHorzBkg','DiffusionCoeffHorTurbVelLimiter','WindDragCoeff','StokesDriftCoeff']:
        if float(root.find('caseDefinitions/constants/'+key).get('value'))!=0:
            raise ValueError(f'{path}: {key} must be zero')
    return method,dt


def sync_cases():
    """Repair dependent settings in source XMLs while preserving all comments."""
    updates=[]
    for method in METHODS:
        for dt in STEPS:
            path=HERE/f'Advetion_Verification_{method}_{dt:03d}.xml'
            _,_,settings=case_settings(path);text=path.read_text()
            for (tag,key,value,attr),expected in settings.items():
                pattern=r'<'+tag+r'\b[^>]*>'
                matches=[m for m in re.finditer(pattern,text)
                         if key is None or re.search(r'\b'+key+r'="'+re.escape(value)+'"',m[0])]
                if len(matches)!=1:raise ValueError(f'{path}: cannot identify {value or tag}')
                match=matches[0]
                changed,count=re.subn(r'(\b'+attr+r'=")[^"]*(")',
                                      lambda m:m[1]+expected+m[2],match[0])
                if count!=1:raise ValueError(f'{path}: missing or duplicate {attr}')
                text=text[:match.start()]+changed+text[match.end():]
            ET.fromstring(text)
            updates.append((path,text))
    # Validate every source exists before writing any source; never touch _out XMLs.
    for path,text in updates:
        if path.read_text()!=text:path.write_text(text)
        check_case(path)
    print('Synchronized 12 source XMLs from their filenames; existing output folders were not changed.')


def analyze_case(method,dt):
    name=f'Advetion_Verification_{method}_{dt:03d}'
    work=HERE/(name+'_out')
    check_case(work/(name+'.xml'))
    trace=[];cohort=None;previous=None
    for file in sorted(work.glob(name+'_[0-9][0-9][0-9][0-9][0-9].vtu')):
        data=read_vtu(file)
        if len(data['id'])!=1: raise ValueError(f'{file}: expected exactly one particle')
        if cohort is None: cohort=data['id'].copy()
        if not np.array_equal(data['id'],cohort): raise ValueError('Particle ID changed')
        if data['state'][0]!=0 or data['particulate'][0]!=0: raise ValueError('Particle must remain passive and in water')
        p=data['position'][0]*np.array([DEGREE,DEGREE,1.]);t=float(data['age'][0])
        if not np.all(np.isfinite(p)): raise ValueError('Nonfinite position')
        if np.any(np.abs(p[:2])>=999) or not -399<p[2]<-1: raise ValueError('Particle reached boundary')
        if previous is not None and abs(t-previous[0])<1e-8:
            if not np.allclose(p,previous[1],rtol=0,atol=1e-10): raise ValueError('Motion without increasing age')
            continue
        if previous is not None and abs(t-previous[0]-dt)>1e-7: raise ValueError('Missing particle output')
        if previous is None and (abs(t)>1e-8 or not np.allclose(p,INITIAL,rtol=0,atol=1e-10)):
            raise ValueError('Initial position/age mismatch')
        analytical=exact(t);err=p-analytical
        trace.append(dict(time_s=t,X=p[0],Y=p[1],Z=p[2],X_exact=analytical[0],Y_exact=analytical[1],Z_exact=analytical[2],error_m=np.linalg.norm(err)))
        previous=(t,p)
    if len(trace)!=round(T/dt)+1 or abs(trace[-1]['time_s']-T)>1e-7: raise ValueError(f'{name}: physical duration is not 3600 seconds')
    with (work/'trajectory.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(trace[0]));w.writeheader();w.writerows(trace)
    write_trajectory_table(work,method,dt,trace)
    e=previous[1]-exact(T)
    result=dict(method=method,dt_s=dt,error_m=float(np.linalg.norm(e)),error_X_m=float(e[0]),error_Y_m=float(e[1]),error_Z_m=float(e[2]),
                relative_X_percent=float(100*abs(e[0])/abs(exact(T)[0]-INITIAL[0])),
                relative_Y_percent=float(100*abs(e[1])/abs(exact(T)[1]-INITIAL[1])),
                relative_Z_percent=float(100*abs(e[2])/abs(exact(T)[2]-INITIAL[2])),
                observed_order=None,particles=1,duration_s=T)
    (work/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def write_trajectory_table(work, method, dt, trace):
    """Compare the actual VTU trajectory with the reference at common ages."""
    labels={'Eluer1':'Euler','Eluer2':'second-order Euler predictor-corrector','RK4':'fourth-order Runge--Kutta'}
    def percent(value):
        text=f'{value:.6g}'
        if 'e' in text:
            mantissa,exponent=text.split('e')
            text=rf'{mantissa}\times10^{{{int(exponent)}}}'
        return '$'+text+'$'
    def coordinates(p):
        return [f'{p[0]/DEGREE:.17f}', f'{p[1]/DEGREE:.17f}', f'{p[2]:.12f}']
    tex=[r'\begin{table}[htbp]',r'\centering',
         r'\caption{Comparison of analytical and numerical particle positions for the nonuniform, time-dependent three-dimensional advection verification test using '+
         labels[method]+rf' with $\Delta t={dt}~\mathrm{{s}}$. All other transport processes are disabled.}}',
         rf'\label{{tab:verification_advection_trajectory_{method}_{dt:03d}}}',
         r'\footnotesize',r'\renewcommand{\arraystretch}{1.15}',
         r'\resizebox{\textwidth}{!}{%',r'\begin{tabular}{c c c c c c c c}',r'\toprule',
         r'Time & Solution & Longitude ($^\circ$) & Latitude ($^\circ$) & $z$ (m) & $Error_{\mathrm{lon}}$ (\%) & $Error_{\mathrm{lat}}$ (\%) & $Error_z$ (\%) \\',r'\midrule']
    for target in [0.,T/2,T]:
        selected=[r for r in trace if np.isclose(float(r['time_s']),target,rtol=0,atol=1e-7)]
        if len(selected)!=1:raise ValueError(f'{work.name}: missing trajectory-table age {target:g} s')
        row=selected[0]
        numerical=np.array([float(row[k]) for k in ['X','Y','Z']])
        analytical=np.array([float(row[k]) for k in ['X_exact','Y_exact','Z_exact']])
        time=f'{target/3600:g} h'
        if target==0:
            cells=[time,'Initial',*coordinates(numerical),'','','']
            tex.append(' & '.join(cells)+r' \\')
        else:
            tex.append(' & '.join([time,'Analytical',*coordinates(analytical),'','',''])+r' \\')
            displacement=analytical-INITIAL
            errors=[percent(100*abs(numerical[j]-analytical[j])/abs(displacement[j]))
                    if displacement[j]!=0 else '--' for j in range(3)]
            tex.append(' & '.join([time,'Numerical',*coordinates(numerical),*errors])+r' \\')
        if target!=T:tex.append(r'\midrule')
    tex.extend([r'\bottomrule',r'\end{tabular}',r'}',r'\par\smallskip\raggedright\footnotesize',
                r'Time is elapsed particle age, not the XML calendar timestamp. '
                r'$Error_j=100|x_j^{\mathrm{num}}-x_j^{\mathrm{ana}}|/|x_j^{\mathrm{ana}}-x_j(0)|$. '
                r'Errors use unrounded positions; -- denotes a zero analytical displacement.',r'\end{table}'])
    (work/'trajectory_table.tex').write_text('\n'.join(tex)+'\n')


def write_convergence_table(results):
    """Write the manuscript table from measured convergence results."""
    labels={'Eluer1':'Euler','Eluer2':'Second-order multi-step Euler','RK4':'Fourth-order Runge--Kutta'}
    tex=[r'\begin{table}[htbp]',r'\centering',
         r'\caption{Temporal convergence of the three MOHID-Lagrangian integration schemes for the nonuniform, time-dependent three-dimensional advection test at $T=3600~\mathrm{s}$. The final-time position error $E_{\Delta t}$ and observed order $p_{\Delta t}$ are defined in Equations~\eqref{eq:verification_advection_error_norm} and \eqref{eq:verification_advection_observed_order}.}',
         r'\label{tab:verification_advection_convergence}',r'\footnotesize',r'\renewcommand{\arraystretch}{1.15}',r'\begin{tabular}{l c c c}',r'\toprule',
         r'\textbf{Integration scheme} & $\boldsymbol{\Delta t}$ (s) & $\boldsymbol{E_{\Delta t}}$ (m) & \textbf{Observed order}, $p_{\Delta t}$ \\',r'\midrule']
    for method in METHODS:
        group=[r for r in results if r['method']==method]
        for i,r in enumerate(group):
            order='--' if r['observed_order'] is None else f"{r['observed_order']:.3f}"
            label=(rf'\multirow{{{len(group)}}}{{*}}{{{labels[method]}}}' if i==0 else '')
            error=f"{r['error_m']:.6g}"
            if 'e' in error:
                mantissa,exponent=error.split('e')
                error=rf'{mantissa}\times10^{{{int(exponent)}}}'
            tex.append(f"{label} & {r['dt_s']} & ${error}$ & {order}"+r' \\')
        if method!='RK4': tex.append(r'\midrule')
    tex.extend([r'\bottomrule',r'\end{tabular}',r'\par\smallskip\raggedright\footnotesize',
                r'$E_{\Delta t}=\|\mathbf X_{num}(T)-\mathbf X_{ana}(T)\|_2$; '
                r'$p_{\Delta t}=\log_2(E_{\Delta t}/E_{\Delta t/2})$. '
                r'Orders describe the current implementation; no theoretical order is imposed.',r'\end{table}'])
    (HERE/'convergence_table.tex').write_text('\n'.join(tex)+'\n')


def collect_trajectory_tables():
    tables=[]
    for method in METHODS:
        for dt in STEPS:
            work=HERE/f'Advetion_Verification_{method}_{dt:03d}_out'
            tables.append((work/'trajectory_table.tex').read_text())
    (HERE/'trajectory_tables.tex').write_text('\n'.join(tables))


def tables_only():
    """Format existing CSV results without running the model or reading VTUs."""
    traces=[]
    for method in METHODS:
        for dt in STEPS:
            work=HERE/f'Advetion_Verification_{method}_{dt:03d}_out'
            with (work/'trajectory.csv').open() as f:
                traces.append((work,method,dt,list(csv.DictReader(f))))
    with (HERE/'convergence.csv').open() as f:
        results=list(csv.DictReader(f))
    for row in results:
        row['dt_s']=int(row['dt_s']);row['error_m']=float(row['error_m'])
        row['observed_order']=float(row['observed_order']) if row['observed_order'] else None
    for work,method,dt,trace in traces:write_trajectory_table(work,method,dt,trace)
    write_convergence_table(results)
    collect_trajectory_tables()
    print('Tables generated from existing CSV results: trajectory_tables.tex and convergence_table.tex')


def analyze():
    results=[]
    # All cases are required; missing cases must not silently produce a partial convergence table.
    for method in METHODS:
        group=[analyze_case(method,dt) for dt in STEPS]
        for coarse,fine in zip(group,group[1:]):
            if fine['error_m']>0 and coarse['error_m']>0:
                coarse['observed_order']=float(np.log2(coarse['error_m']/fine['error_m']))
        for r in group:
            name=f"Advetion_Verification_{method}_{r['dt_s']:03d}"
            (HERE/(name+'_out')/'summary.json').write_text(json.dumps(r,indent=2)+'\n')
        results.extend(group)
    with (HERE/'convergence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    labels={'Eluer1':'Euler','Eluer2':'Second-order multi-step Euler','RK4':'Fourth-order Runge--Kutta'}
    write_convergence_table(results)
    collect_trajectory_tables()
    with plt.rc_context({'font.family':'serif','font.serif':['STIXGeneral'],'mathtext.fontset':'stix','font.size':12,'pdf.fonttype':42}):
        fig,ax=plt.subplots(figsize=(7.4,5.1))
        for method,color,marker in zip(METHODS,['#0072B2','#D55E00','#009E73'],['o','s','^']):
            group=[r for r in results if r['method']==method]
            ax.loglog([r['dt_s'] for r in group],[r['error_m'] for r in group],color=color,marker=marker,
                      markerfacecolor='white',markeredgewidth=1.6,linewidth=1.8,label=labels[method])
        ax.set(xlabel=r'Time step $\Delta t$ (s)',ylabel=r'Final position error $E_{\Delta t}$ (m)')
        ax.set_xticks(sorted(STEPS),[str(s) for s in sorted(STEPS)])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(which='both',direction='in',top=True,right=True)
        ax.grid(which='major',linestyle='--',color='.85');ax.legend(fontsize=10)
        fig.tight_layout()
        for ext in ['png','pdf']:fig.savefig(HERE/('Advection_convergence.'+ext),dpi=600,bbox_inches='tight')
        plt.close(fig)
    for r in results: print(r['method'],r['dt_s'],f"E={r['error_m']:.6g} m",'p='+str(r['observed_order']))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    action=parser.add_mutually_exclusive_group()
    action.add_argument('--tables-only',action='store_true',help='Generate LaTeX tables from existing CSVs only')
    action.add_argument('--prepare-inputs',action='store_true',help='Generate the shared NetCDF only')
    action.add_argument('--sync-cases',action='store_true',help='Set dependent source XML settings from the method and dt in each filename')
    action.add_argument('--check-case',type=Path,help='Validate one XML before running the model')
    parser.add_argument('--latitude-refinement',type=int,default=1,choices=[1,2],help='Refine metric interpolation grid when preparing input')
    args=parser.parse_args()
    try:
        if args.tables_only:tables_only()
        elif args.prepare_inputs:prepare_inputs(args.latitude_refinement)
        elif args.sync_cases:sync_cases()
        elif args.check_case:
            check_case(args.check_case)
            print(f'{args.check_case.name}: configuration OK')
        else:analyze()
    except (ValueError, FileNotFoundError, ET.ParseError) as error:
        parser.exit(1,f'Error: {error}\n')
