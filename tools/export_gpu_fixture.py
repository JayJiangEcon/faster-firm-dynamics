#!/usr/bin/env python3
"""Build the packaged CPU solver with binary exports and prepare one CUDA job."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parameters', type=Path, required=True)
    p.add_argument('--build', type=Path, required=True, help='New fixture-export build directory')
    p.add_argument('--output', type=Path, required=True, help='New fixture output directory')
    p.add_argument('--fc', default=os.environ.get('FC', 'gfortran'))
    p.add_argument('--flags', default='-O3 -funroll-loops -fopenmp -m64')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--prepare-only', action='store_true', help='Instrument source without compiling or evaluating')
    a = p.parse_args()
    point = json.loads(a.parameters.read_text())
    model = point['model']
    if model not in ('RE', 'DE'):
        p.error('Expected model RE or DE')
    build = a.build.resolve()
    output = a.output.resolve()
    if output.exists():
        p.error('Choose a new output directory')
    subprocess.run([sys.executable, str(ROOT/'build.py'), model, '--output', str(build), '--prepare-only'], check=True)
    from gpu_fixture_dumps import transform
    source = build/'diag_firms_solve.f90'
    source.write_text(transform(source.read_text()))
    record_path = build/'build.json'
    record = json.loads(record_path.read_text())
    record['generated_solver_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
    record['recorded_source_match'] = False
    record['instrumentation'] = 'Binary exports after CPU policy and simulation execution'
    record_path.write_text(json.dumps(record, indent=2)+'\n')
    if a.prepare_only:
        print(record_path)
        return
    fc = shutil.which(a.fc)
    if not fc:
        p.error(f'Compiler not found: {a.fc}')
    flags = shlex.split(a.flags)
    env = dict(os.environ)
    env['PATH'] = str(Path(fc).parent)+os.pathsep+env.get('PATH','')
    if sys.platform == 'darwin' and shutil.which('xcrun'):
        env.setdefault('SDKROOT', subprocess.check_output(['xcrun','--show-sdk-path'], text=True).strip())
    binary = build/'fixture-export'
    subprocess.run([fc,*flags,'base_lib.f90','diag_firms_solve.f90','evaluate.f90','-o',str(binary)], cwd=build, env=env, check=True)
    record.update(compiled=True, compiler=fc, flags=flags, binary=binary.name,
                  binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
    record_path.write_text(json.dumps(record, indent=2)+'\n')
    subprocess.run([sys.executable,str(ROOT/'run.py'),'--build',str(build),'--parameters',str(a.parameters.resolve()),'--output',str(output),'--threads',str(a.threads)], check=True)
    x = point['parameters']
    offset = 1 if model == 'DE' else 0
    options = [('eta-k',x[2+offset]),('phi',x[3+offset]),('gamma',x[4+offset]),('eta-f',x[6+offset]),('eta-d',x[7+offset]),('sigma-pi',x[8+offset])]
    # The C++ batch reader accepts std::quoted strings for paths with spaces.
    line = '--fixture '+json.dumps(str(output))+' '+ ' '.join('--'+k+' '+repr(float(v)) for k,v in options)
    (output/'job.txt').write_text(line+'\n')
    files = {f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(output.glob('RCC_*.bin'))}
    (output/'fixture.json').write_text(json.dumps({'model':model,'parameter_hex':[float(v).hex() for v in x],'files':files,'build_record_sha256':hashlib.sha256(record_path.read_bytes()).hexdigest()},indent=2)+'\n')
    print(output/'job.txt')

if __name__ == '__main__':
    main()
