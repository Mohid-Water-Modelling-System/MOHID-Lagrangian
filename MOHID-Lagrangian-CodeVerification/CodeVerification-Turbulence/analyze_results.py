#!/usr/bin/env python3
"""Read MOHID particle VTU outputs and produce statistical tables and figures.

Run the two RunCase_*.sh scripts first. This program does not run the model or
change the case XMLs. --prepare-inputs only regenerates the synthetic NetCDF.
"""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import argparse
import base64
import csv
import hashlib
import json
import os
import platform
import sys
import xml.etree.ElementTree as ET
import numpy as np
from netCDF4 import Dataset
os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ.get("TMPDIR", "/tmp"))/"mohid-verification-matplotlib"))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, ScalarFormatter

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
INPUT_DIR = HERE/'nc_fields_3D/hydrodynamic'
INPUT_FILE = INPUT_DIR/'MOHID_Hydrodynamic_ZeroCurrents_20240101_0000_[actual_time=20240101].nc4'
RADIUS = 6378137.0  # Same radius as utilities.f90 m2geo_comp.
DEGREE = RADIUS * np.pi / 180
# Equatorial center minimizes the local Cartesian longitude approximation error.
CENTER = np.array([0., 0., -100.])
HALF_DEG = 10000 / DEGREE
START = datetime(2024, 1, 1)
CASES = {'rw': ('RandomWalk', 10.), 'rf': ('RandomFlight', 5.)}
FIELDS = ('age', 'dVelX', 'dVelY', 'TPathHor', 'VelStandardDeviation', 'particulate')

def prepare_inputs():
    INPUT_DIR.mkdir(parents=True,exist_ok=True)
    path = INPUT_FILE
    with Dataset(path, 'w', format='NETCDF4_CLASSIC') as ds:
        for dim, values, units in (
            ('time', np.arange(0., 10801., 3600.), 'seconds since 2024-01-01 00:00:00'),
            ('level', [-200., -150., -100., -50., 0.], 'm'),
            ('lat', np.linspace(-HALF_DEG, HALF_DEG, 21), 'degrees_north'),
            ('lon', np.linspace(-HALF_DEG, HALF_DEG, 21), 'degrees_east')):
            ds.createDimension(dim, len(values))
            v = ds.createVariable(dim, 'f8', (dim,)); v[:] = values; v.units = units
            if dim=='time': v.calendar='standard'
        for name, value, units in [('u',0.,'m s-1'),('v',0.,'m s-1'),
                                    ('w',0.,'m s-1'),('bathymetry',-200.,'m')]:
            v=ds.createVariable(name,'f8',('time','level','lat','lon'),
                                fill_value=9.969209968386869e36,zlib=True)
            v[:] = np.full(v.shape,value); v.units=units
        ds.description='All-wet 20 km by 20 km, 200 m deep verification domain; zero u/v/w.'
    return path

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


def plot_msd(kind, work, times, measured, analytical):
    """Publication-style MSD plot; all data points are retained in each curve."""
    name,_=CASES[kind]
    style={
        'font.family':'serif', 'font.serif':['STIXGeneral','DejaVu Serif'],
        'mathtext.fontset':'stix', 'font.size':12,
        'axes.labelsize':14, 'axes.titlesize':15, 'axes.titleweight':'bold',
        'axes.linewidth':1.1, 'xtick.labelsize':12, 'ytick.labelsize':12,
        'xtick.direction':'in', 'ytick.direction':'in',
        'xtick.top':True, 'ytick.right':True,
        'pdf.fonttype':42, 'ps.fonttype':42,
        'savefig.facecolor':'white',
    }
    with plt.rc_context(style):
        fig,ax=plt.subplots(figsize=(7.4,5.1))
        fig.subplots_adjust(left=.15,right=.97,bottom=.15,top=.80)
        marker_indices=np.unique(np.linspace(0,len(times)-1,min(13,len(times)),dtype=int))
        numerical,=ax.plot(times,measured,color='#0072B2',linewidth=2.0,
            marker='o',markevery=marker_indices,markersize=6,
            markerfacecolor='white',markeredgewidth=1.5,
            label='MOHID-Lagrangian',zorder=3)
        reference,=ax.plot(times,analytical,color='black',linewidth=1.7,
            linestyle=(0,(5,3)),label='Analytical solution',zorder=4)
        ax.set_title('Random walk' if kind=='rw' else 'Random flight',pad=12)
        ax.set_xlabel(r'Time, $t$ (s)')
        ax.set_ylabel(r'Horizontal MSD ($\mathrm{m}^2$)')
        span=float(times[-1]-times[0])
        peak=max(float(np.max(measured)),float(np.max(analytical)),1e-12)
        ax.set_xlim(float(times[0])-.025*span,float(times[-1])+.025*span)
        ax.set_ylim(-.025*peak,1.06*peak)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=7))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=7,steps=[1,2,2.5,5,10]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        formatter=ScalarFormatter(useOffset=False); formatter.set_scientific(False)
        ax.yaxis.set_major_formatter(formatter)
        ax.tick_params(which='major',length=5,width=1.0,pad=6)
        ax.tick_params(which='minor',length=2.5,width=.8)
        ax.set_axisbelow(True)
        ax.grid(which='major',color='#D3D3D3',linestyle='--',linewidth=.75,alpha=.85)
        fig.legend(handles=[numerical,reference],loc='upper center',
                   bbox_to_anchor=(.55,.98),ncol=2,frameon=False,
                   handlelength=2.7,columnspacing=2.0,fontsize=12)
        for ext in ('png','pdf'):
            fig.savefig(work/(name+'_MSD_verification.'+ext),dpi=600,bbox_inches='tight',pad_inches=.08)
        plt.close(fig)


def write_verification_table(work, checks, duration=3600.):
    """Write the manuscript table; zero references include one standard error."""
    flight='velocity_mean_1' in checks
    kind='random_flight' if flight else 'random_walk'
    if flight:
        caption=(r'Statistical verification of the piecewise-constant random-flight '
                 r'turbulence model for constant $\sigma_v$ and $T_{eddy}$. '
                 r'Analytical references and numerical results are shown.')
        rows=[
            (r'Mean velocity, $\overline{v_1^{turb}}$', 'velocity_mean_1', r'm\,s^{-1}'),
            (r'Mean velocity, $\overline{v_2^{turb}}$', 'velocity_mean_2', r'm\,s^{-1}'),
            (r'Velocity variance, $\operatorname{Var}(v_1^{turb})$', 'velocity_variance_1', r'm^2\,s^{-2}'),
            (r'Velocity variance, $\operatorname{Var}(v_2^{turb})$', 'velocity_variance_2', r'm^2\,s^{-2}'),
            (r'Cross-covariance, $\operatorname{Cov}(v_1^{turb},v_2^{turb})$', 'velocity_covariance', r'm^2\,s^{-2}'),
            (r'RMS horizontal velocity magnitude', 'rms_horizontal_velocity', r'm\,s^{-1}'),
            (r'Eddy renewal interval', 'eddy_renewal_timing', 's'),
            (r'Recovered horizontal diffusivity, $\widehat D_h^{RF}$', 'recovered_horizontal_diffusivity', r'm^2\,s^{-1}')]
        note=(r'Velocity statistics refer to the turbulent component only. '
              r'The renewal interval is verified against the output timer and velocity changes. '
              r'The recovered diffusivity is estimated from the late-time MSD slope. ')
    else:
        caption=(r'Random-walk statistical verification under zero resolved advection '
                 r'and prescribed spatially uniform diffusivities. Analytical reference '
                 rf'values are evaluated at $T={duration:g}~\mathrm{{s}}$.')
        rows=[
            (r'Mean displacement, $\overline{\Delta X_1}$', 'mean_displacement_1', 'm'),
            (r'Mean displacement, $\overline{\Delta X_2}$', 'mean_displacement_2', 'm'),
            (r'Horizontal variance, $\operatorname{Var}(\Delta X_1)$', 'displacement_variance_1', r'm^2'),
            (r'Horizontal variance, $\operatorname{Var}(\Delta X_2)$', 'displacement_variance_2', r'm^2'),
            (r'Vertical variance, $\operatorname{Var}(\Delta X_3)$', 'displacement_variance_3', r'm^2'),
            (r'Recovered horizontal diffusivity, $\widehat D_h$', 'recovered_horizontal_diffusivity', r'm^2\,s^{-1}')]
        note=''
    def number(value):
        formatted=f'{value:.6g}'
        if 'e' in formatted:
            mantissa,exponent=formatted.split('e')
            return rf'{mantissa}\times 10^{{{int(exponent)}}}'
        if abs(value)>=1000:
            return f'{value:,.2f}'.rstrip('0').rstrip('.').replace(',',r'\,')
        return formatted
    latex=[r'\begin{table}[htbp]', r'\centering',
           r'\caption{'+caption+'}', r'\label{tab:verification_'+kind+'}',
           r'\footnotesize', r'\renewcommand{\arraystretch}{1.15}',
           r'\begin{tabular}{@{}p{5.4cm}p{2.5cm}p{3.5cm}p{1.6cm}@{}}',
           r'\toprule',
           r'\textbf{Statistical quantity} & \textbf{Analytical reference} & '
           r'\textbf{Numerical result} & \textbf{Relative error (\%)} \\', r'\midrule']
    for label,key,unit in rows:
        row=checks[key]
        if key=='eddy_renewal_timing':
            reference=float(row['reference_s'])
            numerical=rf'${number(reference)}~\mathrm{{{unit}}}$ (verified)' if row['passed'] else 'Not verified'
            error='0.00' if row['passed'] else '--'
        else:
            reference=float(row['reference']); estimate=float(row['estimate'])
            value=number(estimate)
            if reference==0:
                uncertainty=number(float(row['standard_error']))
                numerical=(r'\shortstack[l]{'+f'${value}$'+r'\\'+
                           rf'$\pm {uncertainty}~\mathrm{{{unit}}}$'+'}')
            else:
                numerical=rf'${value}~\mathrm{{{unit}}}$'
            error='--' if reference==0 else f'{100*abs(estimate-reference)/abs(reference):.2f}'
        analytical=rf'${number(reference)}~\mathrm{{{unit}}}$'
        latex.append(f'{label} & {analytical} & {numerical} & {error}'+r' \\')
    latex.extend([r'\bottomrule', r'\end{tabular}', r'\vspace{2pt}',
                  r'\par\raggedright\footnotesize\emph{Note:} '+note+
                  r'Relative error is $100|x_{num}-x_{ref}|/|x_{ref}|$. '
                  r'For zero theoretical references, -- denotes an undefined relative error; '
                  r'numerical results are reported as estimate $\pm$ one standard error. '+
                  (r'For mean turbulent velocities, $SE_j=\sqrt{\operatorname{Var}(v_j^{turb})/N}$, '
                   if flight else
                   r'For mean displacements, $SE_j=\sqrt{\operatorname{Var}(\Delta X_j)/N}$, ')+
                  r'where $N$ is the number of independent particles.',
                  r'\end{table}'])
    (work/'verification_table.tex').write_text('\n'.join(latex)+'\n')


def measure_check(checks,name,estimate,target,se):
    tolerance=max(6*float(se),1e-10*max(1.,abs(target)))
    checks[name]={'estimate':float(estimate),'reference':float(target),
                  'standard_error':float(se),'tolerance':tolerance,
                  'passed':bool(abs(estimate-target)<=tolerance)}


def analyze(kind,work,particles,duration):
    name,dt=CASES[kind]; rows=[]; checks={}; first_ids=None; previous_vel=None
    files=[]; previous_positions=None
    previous_t=None; min_clearance=np.inf; rf_timing_ok=True; stationary_ok=True
    last=None; curve_max_standard_errors=0.; curve_ok=True
    for file in sorted(work.glob(name+'_Def_[0-9][0-9][0-9][0-9][0-9].vtu')):
        data=read_vtu(file); pos=data['position']; n=len(pos)
        if n!=particles: raise ValueError(f'{file}: expected {particles}, got {n} particles')
        if first_ids is None: first_ids=data['id'].copy()
        if not np.array_equal(data['id'],first_ids): raise ValueError('Particle cohort changed')
        if not all(np.all(np.isfinite(a)) for a in data.values()): raise ValueError('Nonfinite VTU values')
        if np.any(data['state']!=0): raise ValueError('Beaching/land mask is not water')
        if np.any(data['particulate']!=0): raise ValueError('Particles must be nonparticulate')
        if np.any(data['velocity']!=0): raise ValueError('Resolved velocity is nonzero')
        age=data['age']; t=float(age[0])
        if not np.allclose(age,t,rtol=0,atol=1e-8): raise ValueError('Noninstantaneous release')
        if previous_t is not None and t==previous_t:
            if not np.array_equal(pos,previous_positions): raise ValueError('Motion without increasing particle age')
            continue
        previous_positions=pos.copy()
        files.append(file)
        clearance=np.min(np.column_stack(((HALF_DEG-np.abs(pos[:,0]))*DEGREE,
                        (HALF_DEG-np.abs(pos[:,1]))*DEGREE,pos[:,2]+200,-pos[:,2])))
        min_clearance=min(min_clearance,float(clearance))
        if clearance<=1.: raise ValueError('A particle reached a domain boundary')
        # Equatorial tangent-plane coordinates; finite-latitude metric error < 2e-6 here.
        displacement=(pos-CENTER)*np.array([DEGREE,DEGREE,1.])
        means=displacement.mean(axis=0); variance=displacement.var(axis=0,ddof=1)
        msd=np.mean(np.sum(displacement[:,:2]**2,axis=1))
        if kind=='rw': theory=8*t
        else:
            m=int(t//60); r=t-m*60; theory=.01*(m*60**2+r*r)
            velocity=np.column_stack((data['dVelX'],data['dVelY']))
            expected_timer=(int(round(t/dt))-1)%12*dt+dt
            if t>0: rf_timing_ok &= bool(np.allclose(data['TPathHor'],expected_timer,atol=1e-8,rtol=0))
            if previous_vel is not None and t>0 and expected_timer!=dt:
                stationary_ok &= bool(np.array_equal(velocity,previous_vel))
            if previous_vel is not None and expected_timer==dt:
                rf_timing_ok &= bool(np.any(velocity!=previous_vel))
            if t>0 and not np.allclose(data['VelStandardDeviation'],.1,atol=1e-12,rtol=0):
                raise ValueError('Random-flight sigma is not 0.1')
            previous_vel=velocity.copy()
        square=np.sum(displacement[:,:2]**2,axis=1)
        curve_se=square.std(ddof=1)/np.sqrt(n)
        curve_tolerance=max(6*curve_se,1e-10)
        curve_ok &= bool(abs(msd-theory)<=curve_tolerance)
        if curve_se>0: curve_max_standard_errors=max(curve_max_standard_errors,abs(msd-theory)/curve_se)
        if previous_t is not None and abs(t-previous_t-dt)>1e-8:
            raise ValueError('Missing or irregular output step')
        previous_t=t
        rows.append(dict(time_s=t,mean_x=means[0],mean_y=means[1],mean_z=means[2],
                         var_x=variance[0],var_y=variance[1],var_z=variance[2],
                         msd_h=msd,msd_h_analytical=theory))
        last=(data,displacement)
    if not rows or abs(rows[-1]['time_s']-duration)>1e-8: raise ValueError('Final particle age differs from requested duration')
    if len(rows)!=int(duration/dt)+1: raise ValueError('Unexpected number of VTU outputs')
    checks['msd_time_evolution']={'passed':curve_ok,'max_standard_errors':float(curve_max_standard_errors),
                                  'bound_standard_errors':6}
    data,disp=last; n=particles
    means=disp.mean(axis=0); centered=disp-means; var=disp.var(axis=0,ddof=1)
    for j in range(3 if kind=='rw' else 2):
        measure_check(checks,f'mean_displacement_{j+1}',means[j],0,np.sqrt(var[j]/n))
        target=2*[2,2,1e-5][j]*duration if kind=='rw' else rows[-1]['msd_h_analytical']/2
        se=np.sqrt(max(0.,np.mean(centered[:,j]**4)-var[j]**2)/n)
        measure_check(checks,f'displacement_variance_{j+1}',var[j],target,se)
    squared=np.sum(disp[:,:2]**2,axis=1)
    measure_check(checks,'final_horizontal_msd',squared.mean(),rows[-1]['msd_h_analytical'],squared.std(ddof=1)/np.sqrt(n))
    times=np.array([r['time_s'] for r in rows]); msds=np.array([r['msd_h'] for r in rows])
    # OLS slope on final half of the trajectory. Pointwise curve values are correlated:
    # use independent particle paths, not the OLS iid-residual error, for uncertainty.
    selected=np.flatnonzero(times>=duration/2); ts=times[selected]; centered_t=ts-ts.mean()
    weights=centered_t/np.sum(centered_t**2)/4
    per_particle_D=np.zeros(n)
    for index,weight in zip(selected,weights):
        d=read_vtu(files[index])
        x=(d['position'][:,:2]-CENTER[:2])*DEGREE
        per_particle_D+=weight*np.sum(x*x,axis=1)
    # Discrete synchronized RF has a periodic MSD modulation. Fit the exact reference
    # with the same times/window, rather than assuming an exactly linear finite window.
    refcurve=np.array([r['msd_h_analytical'] for r in rows])
    fit_reference=float(np.dot(weights,refcurve[selected]))
    measure_check(checks,'recovered_horizontal_diffusivity',per_particle_D.mean(),fit_reference,
                  per_particle_D.std(ddof=1)/np.sqrt(n))
    if kind=='rf':
        v=np.column_stack((data['dVelX'],data['dVelY'])); vc=v-v.mean(axis=0)
        for j in range(2):
            measure_check(checks,f'velocity_mean_{j+1}',v[:,j].mean(),0,v[:,j].std(ddof=1)/np.sqrt(n))
            vv=v[:,j].var(ddof=1)
            measure_check(checks,f'velocity_variance_{j+1}',vv,.005,np.sqrt(max(0.,np.mean(vc[:,j]**4)-vv**2)/n))
        measure_check(checks,'velocity_covariance',np.mean(vc[:,0]*vc[:,1]),0,
                      np.std(vc[:,0]*vc[:,1],ddof=1)/np.sqrt(n))
        speed2=np.sum(v*v,axis=1)
        measure_check(checks,'rms_horizontal_velocity',np.sqrt(speed2.mean()),.1,
                      speed2.std(ddof=1)/np.sqrt(n)/(2*np.sqrt(speed2.mean())))
        checks['eddy_renewal_timing']={'passed':rf_timing_ok,'reference_s':60}
        checks['velocity_retained_between_renewals']={'passed':stationary_ok}
        checks['vertical_displacement_zero']={'passed':bool(np.all(disp[:,2]==0))}
    table=work/'statistics.csv'
    with table.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    plot_msd(kind,work,times,msds,refcurve)
    write_verification_table(work,checks,duration)
    report={'case':kind,'passed':all(r['passed'] for r in checks.values()),'particles':n,
            'duration_s':duration,'dt_s':dt,'minimum_boundary_clearance_m':min_clearance,
            'long_time_diffusivity_reference':2. if kind=='rw' else .15,
            'finite_window_fit_reference':fit_reference,'checks':checks,
            'uncertainty':'Independent particles; 6 empirical standard errors. Slope SE from full particle paths, not correlated-time OLS residuals.'}
    (work/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f"{name}: {'PASS' if report['passed'] else 'FAIL'}, MSD={msds[-1]:.5g}, D={per_particle_D.mean():.5g}",flush=True)
    return report



def read_specification(kind, work):
    """Read the actual executed XML; never rewrite it or assume CLI overrides."""
    name,dt=CASES[kind]
    root=ET.parse(work/(name+'_Def.xml')).getroot()
    parameters={n.get('key'):n.get('value') for n in root.findall('execution/parameters/parameter')}
    if parameters['Integrator']!='1': raise ValueError('Verification requires Euler integration')
    if float(root.find('caseDefinitions/simulation/timestep').get('dt'))!=dt:
        raise ValueError('Time step differs from the analytical test specification')
    if int(root.find('caseDefinitions/simulation/DiffusionMethod').get('value'))!=(1 if kind=='rw' else 2):
        raise ValueError('Wrong DiffusionMethod for this test')
    values=({'DiffusionCoeffHorz':2.,'DiffusionCoeffHorzBkg':2.,
        'DiffusionCoeffHorTurbVelLimiter':0.,'DiffusionCoeffVert':1e-5} if kind=='rw' else
        {'VarVelH':.1,'VarVelHX':0.,'ConstantMixingLength':6.})
    for tag,value in values.items():
        if float(root.find('caseDefinitions/constants/'+tag).get('value'))!=value:
            raise ValueError('XML differs from analytical verification specification: '+tag)
    start=datetime.strptime(parameters['Start'],'%Y %m %d %H %M %S')
    end=datetime.strptime(parameters['End'],'%Y %m %d %H %M %S')
    duration=(end-start).total_seconds()-2*dt
    if duration<2*dt or duration%dt: raise ValueError('Invalid duration including solver guard steps')
    source=root.find('caseDefinitions/sourceDefinitions/source')
    per_step=float(source.find('rate_trcPerEmission').get('value'))*dt/float(source.find('rate_seconds').get('value'))
    particles=int(np.floor(1+per_step))
    if particles<100: raise ValueError('Statistical verification needs at least 100 particles')
    return particles,duration,parameters


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=('both','rw','rf'),default='both')
    p.add_argument('--input',type=Path,help='One existing output directory; use with --case rw or rf')
    p.add_argument('--prepare-inputs',action='store_true',help='Create zero-current NetCDF only; do not run or analyze')
    p.add_argument('--figures-only',action='store_true',help='Redraw figures from existing statistics.csv files')
    args=p.parse_args()
    if args.prepare_inputs:
        print('Prepared:',prepare_inputs()); return
    if args.input and args.case=='both': p.error('--input requires --case rw or rf')
    kinds=('rw','rf') if args.case=='both' else (args.case,)
    failed=False
    for kind in kinds:
        name,_=CASES[kind]
        work=args.input.resolve() if args.input else HERE/(name+'_Def_out')
        if not work.is_dir(): p.error(f'No output folder {work}; run RunCase_{name}.sh first')
        if args.figures_only:
            data=np.genfromtxt(work/'statistics.csv',delimiter=',',names=True)
            plot_msd(kind,work,data['time_s'],data['msd_h'],data['msd_h_analytical'])
            print('Updated figure:',work/(name+'_MSD_verification.pdf'),flush=True)
            continue
        particles,duration,parameters=read_specification(kind,work)
        model_log=work/'model.log'
        if model_log.exists() and 'Simulation ended' not in model_log.read_text(errors='replace'):
            raise RuntimeError(f'Model did not complete; inspect {model_log}')
        report=analyze(kind,work,particles,duration)
        metadata={'analyzed_utc':datetime.now(timezone.utc).isoformat(),
                  'python':sys.version,'platform':platform.platform(),
                  'threads':parameters.get('Threads'),'num_blocks':parameters.get('NumBlocks'),
                  'particles':particles,'duration_s':duration,
                  'rng':'Production intrinsic RNG; no seed/replay changes',
                  'executed_xml_sha256':hashlib.sha256((work/(name+'_Def.xml')).read_bytes()).hexdigest()}
        local_input=work/'nc_fields_3D/hydrodynamic'/INPUT_FILE.name
        if local_input.exists(): metadata['input_sha256']=hashlib.sha256(local_input.read_bytes()).hexdigest()
        fingerprint=work/'executable.sha256'
        if fingerprint.exists(): metadata['executable_sha256_record']=fingerprint.read_text().strip()
        (work/'analysis_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
        failed |= not report['passed']
        print('Table and figures:',work,flush=True)
    if failed: raise SystemExit('Statistical verification failed; inspect the output summary.json')

if __name__=='__main__': main()
