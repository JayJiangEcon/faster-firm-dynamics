#!/usr/bin/env python3
"""Build, evaluate, or search with the version 0.2.0 CUDA pipeline."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

from cuda_pipeline import ROOT, VERSION, Pipeline, digest, write_json


def build(args):
    output = args.output.resolve()
    if output.exists():
        raise ValueError('Choose a new build directory')
    fc = shutil.which(args.fc)
    nvcc = None if args.inputs_only else shutil.which(args.nvcc)
    if not fc or (not args.inputs_only and not nvcc):
        raise ValueError('GNU Fortran and, for a full build, nvcc must be available')
    if not str(args.arch).isdigit():
        raise ValueError('CUDA architecture must be numeric, e.g. 80')
    manifest = json.loads((ROOT/'results/integrated-manifest.json').read_text())
    for name, expected in manifest['models'][args.model]['files'].items():
        if digest(ROOT/name) != expected:
            raise ValueError(f'CPU source differs from the recorded package: {name}')
    output.mkdir(parents=True)
    cpp, fdir = output/'cuda', output/'fortran'
    cpp.mkdir(); fdir.mkdir()
    source_names = ['src/fortran/base_lib.f90', f'src/fortran/{args.model}/diag_firms_solve.f90',
                    'src/fortran/prepare_cuda_inputs.F90', 'src/fortran/finalize_numeric.f90']
    for name in source_names:
        shutil.copyfile(ROOT/name, fdir/Path(name).name)
    # Compile only the parameter declarations and two input-construction
    # routines. The preparer binary contains no policy solver or objective.
    source = (fdir/'diag_firms_solve.f90').read_text()
    params = source[:source.index('end module params')+len('end module params')]
    def routine(start, end):
        begin = source.index(start)
        finish = source.index(end, begin)+len(end)
        return source[begin:finish]
    input_module = params+'\nmodule solve\ncontains\n'+routine(
        'subroutine discretize()', 'end subroutine discretize')+'\n'+routine(
        'double precision function outputnetlabor(', 'end function outputnetlabor')+'\nend module solve\n'
    (fdir/'input_model.f90').write_text(input_module)
    prefix = args.model.lower()
    native = 're_native_parallel.cu' if args.model == 'RE' else 'de_native_simulation_parallel.cu'
    cuda_names = [native, f'{prefix}_action_block_reduce.cu', f'{prefix}_simulation_parallel.cu']
    for name in cuda_names:
        shutil.copyfile(ROOT/'src/cuda'/name, cpp/name)
        source_names.append('src/cuda/'+name)
    flags = shlex.split(args.flags)
    env = dict(os.environ)
    env['PATH'] = str(Path(fc).parent)+os.pathsep+env.get('PATH', '')
    if sys.platform == 'darwin' and shutil.which('xcrun'):
        env.setdefault('SDKROOT', subprocess.check_output(['xcrun','--show-sdk-path'], text=True).strip())
    command = [fc, *flags, '-cpp', *(['-DMODEL_DE'] if args.model == 'DE' else []),
               'base_lib.f90', 'input_model.f90', 'prepare_cuda_inputs.F90',
               '-o', str(output/'prepare-inputs')]
    with (output/'fortran-build.log').open('w') as log:
        subprocess.run(command, cwd=fdir, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    binaries = {'prepare-inputs': digest(output/'prepare-inputs')}
    finisher = 'finalize-numeric.dylib' if sys.platform == 'darwin' else 'finalize-numeric.so'
    link_flags = ['-dynamiclib'] if sys.platform == 'darwin' else ['-shared','-fPIC']
    with (output/'finisher-build.log').open('w') as log:
        subprocess.run([fc,*flags,*link_flags,'finalize_numeric.f90','-o',str(output/finisher)],
                       cwd=fdir,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    binaries[finisher] = digest(output/finisher)

    cuda_flags = ['-O3','-std=c++17','--default-stream','per-thread','--fmad=false',
                  '--prec-div=true','--prec-sqrt=true','-Xcompiler=-pthread',
                  '-Xcompiler=-ffp-contract=off','-Xcompiler=-fno-fast-math',
                  '-Xcompiler=-march=native',f'-arch=sm_{args.arch}']
    if nvcc:
        with (output/'cuda-build.log').open('w') as log:
            subprocess.run([nvcc,*cuda_flags,native,'-o',str(output/'evaluate-cuda')],
                           cwd=cpp, stdout=log, stderr=subprocess.STDOUT, check=True)
        binaries['evaluate-cuda'] = digest(output/'evaluate-cuda')
    record = {'version': VERSION, 'kind': 'cuda-inputs-only' if args.inputs_only else 'cuda-pipeline',
              'model': args.model, 'numerical_finisher': finisher, 'input_module_sha256': digest(fdir/'input_model.f90'), 'source_sha256': {name: digest(ROOT/name) for name in source_names},
              'binary_sha256': binaries, 'Fortran_compiler': subprocess.check_output([fc,'--version'],text=True),
              'Fortran_flags': flags, 'CUDA_flags': cuda_flags,
              'CUDA_compiler': subprocess.check_output([nvcc,'--version'],text=True) if nvcc else None}
    write_json(output/'build.json', record)
    print(output/'build.json')


def check_search_spec(spec, model):
    import numpy as np
    n = {'RE': 8, 'DE': 9}[model]
    bounds = np.asarray(spec['bounds'], dtype=float)
    sigma = np.asarray(spec['sigma_pi_bounds'], dtype=float)
    if spec['model'] != model or bounds.shape != (n,2) or not np.isfinite(bounds).all() or not (bounds[:,0]<bounds[:,1]).all():
        raise ValueError(f'Expected {n} ordered finite bounds for {model}')
    if sigma.shape != (2,) or not np.isfinite(sigma).all() or not 0<=sigma[0]<sigma[1]:
        raise ValueError('Invalid conditional sigma bounds')
    offset = int(model=='DE')
    if bounds[offset,0]<=-1 or bounds[offset,1]>=1 or bounds[1+offset,0]<=0 or bounds[5+offset,0]<=0:
        raise ValueError('Bounds must keep |rho_z| < 1, sigma_z > 0 and sigma_A > 0')
    for key, default, minimum in [('seed',None,0),('max_evaluations',None,1),('population_multiplier',5,1)]:
        value=spec.get(key,default)
        if type(value) is not int or value<minimum:
            raise ValueError(f'Invalid {key}')
    for key,default in [('tol',1e-7),('atol',0)]:
        value=spec.get(key,default)
        if not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
            raise ValueError(f'Invalid {key}')
    if 'initial_population' in spec:
        initial=np.asarray(spec['initial_population'],dtype=float)
        if initial.ndim!=2 or initial.shape[1]!=n or len(initial)<5 or not np.isfinite(initial).all():
            raise ValueError('Initial population must have at least five finite complete policy vectors')
        if not ((initial>=bounds[:,0]) & (initial<=bounds[:,1])).all():
            raise ValueError('Initial population must lie inside the declared bounds')
    return bounds, sigma


def run_search(engine, spec, destination):
    import numpy as np
    import scipy
    from scipy.optimize import differential_evolution
    bounds, sigma = check_search_spec(spec, engine.model)
    budget = spec['max_evaluations']
    record = {'version': VERSION, 'model': engine.model, 'spec': spec,
              'algorithm': 'scipy.optimize.differential_evolution', 'scipy_version': scipy.__version__,
              'strategy': 'best1bin', 'updating': 'deferred', 'vectorized': True,
              'workers': 1, 'polish': False, 'attempted_evaluations': 0,
              'completed_evaluations': 0, 'unique_policy_vectors': 0, 'generations_completed': 0,
              'best': None, 'status': 'running', 'numerical_scope': 'profiled moments 2 through 17'}
    unique=set(); start=time.perf_counter(); destination=Path(destination)
    def save():
        record['search_seconds']=time.perf_counter()-start
        write_json(destination,record)
    class BudgetReached(Exception):
        pass
    def objective(columns):
        values=np.asarray(columns,dtype=float)
        if values.ndim!=2 or values.shape[0]!=len(bounds):
            raise ValueError('Expected a vectorized candidate matrix')
        points=values.T
        remaining=budget-record['attempted_evaluations']
        if remaining<=0:raise BudgetReached()
        accepted=points[:remaining]
        before=record['attempted_evaluations']
        record['attempted_evaluations']+=len(accepted);save()
        rows=engine.evaluate(accepted,sigma_bounds=sigma)
        if len(rows)!=len(accepted):raise RuntimeError('Evaluator lost candidate results')
        result=[]
        for index,(point,row) in enumerate(zip(accepted,rows)):
            record['completed_evaluations']+=1
            unique.add(tuple(float(v).hex() for v in point));record['unique_policy_vectors']=len(unique)
            q=row['objective'] if row['status']=='ok' else math.inf
            result.append(q)
            if math.isfinite(q) and (record['best'] is None or q<record['best']['objective']):
                record['best']={'evaluation':before+index,'objective':q,
                                'parameters':row['parameters'],'parameter_hex':row['parameter_hex']}
        save()
        if len(accepted)<len(points):raise BudgetReached()
        return np.asarray(result)
    def callback(x,convergence):
        record['generations_completed']+=1;save()
    save()
    try:
        result=differential_evolution(objective,bounds,seed=spec['seed'],
            popsize=spec.get('population_multiplier',5),maxiter=budget,
            init=np.asarray(spec['initial_population']) if 'initial_population' in spec else 'latinhypercube',
            polish=False,workers=1,vectorized=True,updating='deferred',callback=callback,
            tol=spec.get('tol',1e-7),atol=spec.get('atol',0))
        record.update(status='finished',success=bool(result.success),message=str(result.message))
    except BudgetReached:
        record.update(status='evaluation_budget_reached',success=False)
    except BaseException as error:
        record.update(status='failed',success=False,message=str(error));save();raise
    save()
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('build');p.add_argument('model',choices=['RE','DE'])
    p.add_argument('--output',type=Path,required=True);p.add_argument('--fc',default=os.environ.get('FC','gfortran'))
    p.add_argument('--nvcc',default=os.environ.get('NVCC','nvcc'));p.add_argument('--arch',default='80')
    p.add_argument('--flags',default='-O3 -funroll-loops -fopenmp -m64')
    p.add_argument('--inputs-only',action='store_true',help='Build only the CPU input constructor')
    for name in ['evaluate','search']:
        p=sub.add_parser(name);p.add_argument('--build',type=Path,required=True)
        p.add_argument('--output',type=Path,required=True);p.add_argument('--inputs-dir',type=Path,default=ROOT/'data/fixed')
        p.add_argument('--batch-size',type=int,default=32);p.add_argument('--simulation-workers',type=int,default=4)
        p.add_argument('--timeout',type=float,default=300)
        p.add_argument('--remote-host');p.add_argument('--remote-port',type=int,default=22)
        p.add_argument('--remote-key');p.add_argument('--remote-known-hosts');p.add_argument('--remote-build')
        if name=='evaluate':
            p.add_argument('--parameters',type=Path,required=True)
            p.add_argument('--sigma-bounds',type=float,nargs=2)
        else:p.add_argument('--spec',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='build':build(args);return
    build_record=json.loads((args.build/'build.json').read_text())
    if args.command=='search':
        spec=json.loads(args.spec.read_text());check_search_spec(spec,build_record['model'])
    else:
        points=json.loads(args.parameters.read_text())
        if points['model']!=build_record['model']:raise ValueError('Parameter model does not match build')
        points=points.get('points',[points.get('parameters')])
    remote=None
    if args.remote_host:
        if not args.remote_build:raise ValueError('--remote-build is required with --remote-host')
        remote={'host':args.remote_host,'port':args.remote_port,'key':args.remote_key,
                'known_hosts':args.remote_known_hosts,'build':args.remote_build}
    start=time.perf_counter()
    with Pipeline(args.build,args.output,args.inputs_dir,args.simulation_workers,args.batch_size,args.timeout,remote) as engine:
        if args.command=='evaluate':
            result={'model':engine.model,'results':engine.evaluate(points,args.sigma_bounds)}
        else:result=run_search(engine,spec,args.output/'search.json')
        result['total_pipeline_seconds']=time.perf_counter()-start
        result['startup_seconds']=engine.startup_seconds
        write_json(args.output/('results.json' if args.command=='evaluate' else 'search.json'),result)
    print(args.output/('results.json' if args.command=='evaluate' else 'search.json'))

if __name__=='__main__':
    main()
