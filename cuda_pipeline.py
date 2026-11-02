"""Batch input construction and a persistent, reference-free CUDA evaluator."""
from __future__ import annotations

import hashlib
import ctypes
import json
import math
import os
from pathlib import Path
import queue
import shlex
import tarfile
import tempfile
import subprocess
import threading
import time

ROOT = Path(__file__).resolve().parent
VERSION = '0.2.0'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


class StreamProcess:
    """One line per request; JSON completion records delimit each response."""
    def __init__(self, command, cwd, log, timeout=300, env=None):
        self.timeout = timeout
        self.records = queue.Queue()
        self.log = Path(log)
        self.errors = self.log.with_suffix('.stderr.log').open('w')
        self.process = subprocess.Popen(command, cwd=cwd, env=env, stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, stderr=self.errors,
                                        text=True, bufsize=1)
        self.reader = threading.Thread(target=self._read, daemon=True)
        self.reader.start()
        try:
            self.until('ready')
        except BaseException:
            self.close()
            raise

    def _read(self):
        try:
            with self.log.open('w') as log:
                for line in self.process.stdout:
                    log.write(line)
                    log.flush()
                    if line.lstrip().startswith('{'):
                        try:
                            self.records.put(json.loads(line))
                        except ValueError as error:
                            self.records.put(error)
        finally:
            self.records.put(None)

    def until(self, family):
        result = []
        deadline = time.monotonic() + self.timeout
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(f'Worker timed out; see {self.log}')
            try:
                item = self.records.get(timeout=remaining)
            except queue.Empty:
                raise TimeoutError(f'Worker timed out; see {self.log}') from None
            if item is None:
                raise RuntimeError(f'Worker exited before {family}; see {self.log}')
            if isinstance(item, Exception):
                raise RuntimeError(f'Invalid worker JSON; see {self.log}') from item
            result.append(item)
            if item.get('family') == family:
                if item.get('status', 0):
                    raise RuntimeError(f'Worker failed; see {self.log}')
                return result

    def request(self, path):
        text = str(path)
        if '\n' in text or '\r' in text:
            raise ValueError('Newlines in request paths are unsupported')
        self.process.stdin.write(text + '\n')
        self.process.stdin.flush()
        return self.until('batch_complete')

    def close(self):
        if self.process.poll() is None:
            try:
                self.process.stdin.write('QUIT\n')
                self.process.stdin.flush()
                self.process.wait(timeout=5)
            except (OSError, subprocess.TimeoutExpired):
                self.process.kill()
                self.process.wait()
        self.reader.join(timeout=5)
        self.process.stdin.close()
        self.process.stdout.close()
        self.errors.close()


def validate_parameters(model, parameters):
    values = [float(x) for x in parameters]
    if len(values) != {'RE': 9, 'DE': 10}[model] or not all(map(math.isfinite, values)):
        raise ValueError('Expected a complete finite parameter vector')
    offset = int(model == 'DE')
    if abs(values[offset]) >= 1 or values[1+offset] <= 0 or values[5+offset] <= 0:
        raise ValueError('Require |rho_z| < 1, sigma_z > 0 and sigma_A > 0')
    return values


class InputPreparer:
    def __init__(self, model, binary, inputs, directory, timeout=300):
        self.model, self.inputs = model, Path(inputs).resolve()
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True)
        self.shared = self.directory / 'shared'
        self.shared.mkdir()
        env = {k: v for k, v in os.environ.items() if not k.startswith('RCC_')}
        env.update(OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE')
        self.worker = StreamProcess([str(Path(binary).resolve()), '--serve', str(self.shared)],
                                    self.inputs, self.directory/'prepare.stdout.log', timeout, env)
        self.sequence = 0

    def prepare(self, points, directory, profile=False):
        directory = Path(directory).resolve()
        directory.mkdir(parents=True)
        manifest, jobs, bindings = [], [], []
        for i, point in enumerate(points):
            values = validate_parameters(self.model, point)
            case = directory / f'candidate-{i:04d}'
            case.mkdir()
            if '\n' in str(case) or '\r' in str(case):
                raise ValueError('Newlines in paths are unsupported')
            case.joinpath('RCC_pifirmshocks.bin').symlink_to(self.shared/'RCC_pifirmshocks.bin')
            manifest.append('"' + str(case).replace('"', '""') + '" ' + ' '.join(map(repr, values)))
            offset = int(self.model == 'DE')
            options = [('eta-k', values[2+offset]), ('phi', values[3+offset]),
                       ('gamma', values[4+offset]), ('eta-f', values[6+offset]),
                       ('eta-d', values[7+offset]), ('sigma-pi', values[8+offset])]
            jobs.append('--production ' + ('--sigma-profile ' if profile else '') +
                        '--fixture ' + json.dumps(str(case), ensure_ascii=False) + ' ' +
                        ' '.join('--'+key+' '+repr(value) for key, value in options))
            bindings.append({'parameters': values, 'parameter_hex': [v.hex() for v in values]})
        manifest_path = directory/'prepare.txt'
        manifest_path.write_text('\n'.join(manifest)+'\n')
        started = time.perf_counter()
        self.worker.request(manifest_path)
        elapsed = time.perf_counter()-started
        for i, binding in enumerate(bindings):
            case = directory/f'candidate-{i:04d}'
            binding['input_sha256'] = {p.name: digest(p) for p in sorted(case.glob('*.bin'))}
        path = directory/'jobs.txt'
        path.write_text('\n'.join(jobs)+'\n')
        write_json(directory/'inputs.json', {'model': self.model, 'full_CPU_solves': 0,
                    'preparation_seconds': elapsed, 'candidates': bindings})
        self.sequence += 1
        return path, elapsed

    def close(self):
        self.worker.close()


class RemoteGPU:
    """SSH transport for local preprocessing and a persistent remote CUDA worker."""
    def __init__(self, config, output, workers, timeout):
        self.timeout = timeout
        self.ssh = ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15']
        if config.get('key'):
            self.ssh += ['-i', str(Path(config['key']).expanduser().resolve())]
        if config.get('known_hosts'):
            self.ssh += ['-o', 'UserKnownHostsFile='+str(Path(config['known_hosts']).resolve())]
        self.ssh += ['-p', str(config.get('port',22)), config['host']]
        build = str(config['build']).rstrip('/')
        self.record = json.loads(self.command(['cat', build+'/build.json']))
        if self.record.get('kind') != 'cuda-pipeline' or self.record.get('version') != VERSION:
            raise ValueError('Expected a version 0.2.0 remote CUDA build')
        binary = build+'/evaluate-cuda'
        actual = self.command(['sha256sum', binary]).split()[0]
        if actual != self.record['binary_sha256']['evaluate-cuda']:
            raise ValueError('Remote CUDA binary identity changed')
        self.directory = self.command(['mktemp', '-d', '/tmp/ffd-cuda-XXXXXXXX']).strip()
        self.command(['mkdir', self.directory+'/shared'])
        remote_command = 'cd '+shlex.quote(self.directory)+' && exec '+shlex.join([
            binary,'--serve','--simulation-workers',str(workers)])
        self.worker = StreamProcess(self.ssh+[remote_command], output, output/'gpu.stdout.log', timeout)

    def command(self, args, data=None):
        completed = subprocess.run(self.ssh+[shlex.join(args)],input=data,
                                   stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=self.timeout)
        if completed.returncode:
            raise RuntimeError('Remote command failed: '+completed.stderr.decode(errors='replace'))
        return completed.stdout.decode()

    def shared(self, path):
        with tempfile.TemporaryFile() as archive:
            with tarfile.open(fileobj=archive, mode='w') as bundle:
                bundle.add(path, arcname='shared/RCC_pifirmshocks.bin')
            archive.seek(0)
            self.command(['tar','-xf','-','-C',self.directory], archive.read())

    def request(self, jobs):
        # Only model inputs and a manifest cross the transport boundary.
        # Reference solutions and previous optimizer outputs are never uploaded.
        local = Path(jobs).parent
        relative = local.name
        remote = self.directory+'/'+relative
        with tempfile.TemporaryFile() as archive:
            with tarfile.open(fileobj=archive, mode='w') as bundle:
                for case in sorted(local.glob('candidate-*')):
                    for item in sorted(case.glob('*.bin')):
                        name=relative+'/'+case.name+'/'+item.name
                        if item.name=='RCC_pifirmshocks.bin':
                            info=tarfile.TarInfo(name);info.type=tarfile.SYMTYPE
                            info.linkname='../../shared/RCC_pifirmshocks.bin';bundle.addfile(info)
                        else:bundle.add(item,arcname=name)
                import io
                content=Path(jobs).read_text().replace(str(local),remote).encode()
                info=tarfile.TarInfo(relative+'/jobs.txt');info.size=len(content)
                bundle.addfile(info,io.BytesIO(content))
            archive.seek(0)
            self.command(['tar','-xf','-','-C',self.directory],archive.read())
        return self.worker.request(remote+'/jobs.txt')

    def close(self):
        self.worker.close()
        # Keep uploaded input manifests for inspection until the remote host is
        # released; no remote deletion is implicit in closing the worker.


class Pipeline:
    def __init__(self, build, output, inputs=ROOT/'data/fixed', workers=4,
                 batch_size=32, timeout=300, remote=None):
        import numpy as np
        self.np = np
        self.build = Path(build).resolve()
        self.output = Path(output).resolve()
        if self.output.exists():
            raise ValueError('Choose a new output directory')
        self.record = json.loads((self.build/'build.json').read_text())
        self.model = self.record['model']
        if self.record.get('version') != VERSION or self.record.get('kind') not in (
                ('cuda-pipeline','cuda-inputs-only') if remote else ('cuda-pipeline',)):
            raise ValueError('Expected a version 0.2.0 CUDA pipeline build')
        for name, expected in self.record['binary_sha256'].items():
            if digest(self.build/name) != expected:
                raise ValueError(f'Changed build binary: {name}')
        if not 1 <= batch_size <= 64 or not 1 <= workers <= 64 or not math.isfinite(timeout) or timeout <= 0:
            raise ValueError('Use batch size/workers in 1..64 and a positive timeout')
        inputs = Path(inputs).resolve()
        manifest = json.loads((ROOT/'results/integrated-manifest.json').read_text())
        for name, expected in manifest['runtime_input_sha256'].items():
            if digest(inputs/name) != expected:
                raise ValueError(f'Fixed input does not match the recorded contract: {name}')
        library = self.record.get('numerical_finisher')
        if not library:
            raise ValueError('Rebuild the local inputs with the version 0.2.0 numerical finisher')
        self.finalizer = ctypes.CDLL(str(self.build/library))
        array = np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS')
        self.finalizer.macro_variance.argtypes = [array,ctypes.POINTER(ctypes.c_double)]
        self.finalizer.macro_variance.restype = None
        self.finalizer.direct_objective.argtypes = [array,array,array,ctypes.POINTER(ctypes.c_double)]
        self.finalizer.direct_objective.restype = None
        self.target = np.loadtxt(inputs/'_GMMdatamom.txt', delimiter=',')
        self.weight = np.asfortranarray(np.loadtxt(inputs/'_GMMdatawgt.txt', delimiter=','))
        self.output.mkdir(parents=True)
        self.batch_size, self.sequence = batch_size, 0
        self.preparer = self.gpu = None
        start = time.perf_counter()
        try:
            self.preparer = InputPreparer(self.model, self.build/'prepare-inputs', inputs,
                                         self.output/'preparer', timeout)
            if remote:
                self.gpu = RemoteGPU(remote, self.output, workers, timeout)
                if self.gpu.record['model'] != self.model:
                    raise ValueError('Remote and local model identities differ')
                self.gpu.shared(self.preparer.shared/'RCC_pifirmshocks.bin')
            else:
                env = {k: v for k, v in os.environ.items() if not k.startswith('RCC_')}
                self.gpu = StreamProcess([str(self.build/'evaluate-cuda'), '--serve',
                                          '--simulation-workers', str(workers)], self.output,
                                         self.output/'gpu.stdout.log', timeout, env)
        except BaseException:
            self.close()
            raise
        self.startup_seconds = time.perf_counter()-start
        write_json(self.output/'pipeline.json', {'version': VERSION, 'model': self.model,
                  'startup_seconds': self.startup_seconds, 'batch_size': batch_size,
                  'simulation_workers': workers, 'build_sha256': digest(self.build/'build.json'),
                  'fixed_input_sha256': manifest['runtime_input_sha256'],
                  'driver_sha256': {name: digest(ROOT/name) for name in
                                    ('gpu.py','cuda_pipeline.py','sigma_profile.py')},
                  'CPU_policy_solves': 0, 'reference_outputs_required': False,
                  'transport': 'ssh' if remote else 'local',
                  'remote_build': self.gpu.record if remote else None,
                  'remote_input_directory': self.gpu.directory if remote else None})

    def evaluate(self, points, sigma_bounds=None):
        from sigma_profile import minimize_basis
        np = self.np
        n = {'RE': 9, 'DE': 10}[self.model]
        profile = sigma_bounds is not None
        values = np.asarray(points, dtype=float)
        if values.ndim != 2 or values.shape[1] != n-int(profile) or not np.isfinite(values).all():
            raise ValueError('Candidate matrix has the wrong dimension or nonfinite entries')
        if profile:
            bounds = np.asarray(sigma_bounds, dtype=float)
            if bounds.shape != (2,) or not np.isfinite(bounds).all() or not 0 <= bounds[0] < bounds[1]:
                raise ValueError('Expected finite ordered nonnegative sigma bounds')
        rows = []
        for offset in range(0, len(values), self.batch_size):
            chunk = values[offset:offset+self.batch_size]
            complete = [list(x)+([0.0] if profile else []) for x in chunk]
            directory = self.output/f'batch-{self.sequence:06d}'
            self.sequence += 1
            start = time.perf_counter()
            jobs, prep_seconds = self.preparer.prepare(complete, directory, profile)
            gpu_start = time.perf_counter()
            reports = self.gpu.request(jobs)
            gpu_seconds = time.perf_counter()-gpu_start
            result = [r for r in reports if r.get('family') == 'evaluation']
            if len(result) != len(chunk) or [r['job'] for r in result] != list(range(len(chunk))):
                raise RuntimeError('GPU response count or candidate order changed')
            current = []
            for x, observation in zip(complete, result):
                if observation['profile'] != profile:
                    raise RuntimeError('GPU returned the wrong evaluation mode')
                row = {'parameters': x, 'parameter_hex': [float(v).hex() for v in x]}
                try:
                    macro = np.asarray(observation['macro_output'],dtype=float)
                    if macro.shape != (302,) or not np.isfinite(macro[48:300]).all() or not (macro[48:300]>0).all():
                        raise ValueError('Invalid macro output path')
                    variance = ctypes.c_double()
                    self.finalizer.macro_variance(macro,ctypes.byref(variance))
                    if profile:
                        m0, mp, mm = [np.asarray(observation[k], dtype=float) for k in ('m0','mp','mm')]
                        for basis in (m0,mp,mm): basis[16] = variance.value
                        c1, c2 = 0.5*(mp-mm), 0.5*(mp+mm)-m0
                        selected = minimize_basis(m0, c1, c2, self.target, self.weight, sigma_bounds)
                        row.update(selected)
                        row['parameters'] = x[:-1]+[selected['sigma_pi']]
                        row['parameter_hex'] = [float(v).hex() for v in row['parameters']]
                        row['basis'] = {key: [float(v) if math.isfinite(v) else None for v in vector]
                                        for key, vector in zip(('m0','mp','mm','c1','c2'),(m0,mp,mm,c1,c2))}
                    else:
                        m = np.asarray(observation['moments'], dtype=float)
                        if m.shape != (17,) or not np.isfinite(m[1:]).all():
                            raise ValueError('Nonfinite objective moments')
                        m[16] = variance.value
                        objective = ctypes.c_double()
                        self.finalizer.direct_objective(m,self.target,self.weight,ctypes.byref(objective))
                        q = objective.value
                        if not math.isfinite(q):
                            raise ValueError('Nonfinite objective')
                        row.update(objective=q, moments=[float(v) if math.isfinite(v) else None for v in m])
                    row['status'] = 'ok'
                except ValueError as error:
                    # Keep the full diagnostic output, but only selected moments
                    # determine whether the configured objective is finite.
                    row.update(status='invalid_objective', objective=None, error=str(error))
                current.append(row)
            write_json(directory/'results.json', {'model': self.model, 'profile': profile,
                       'preparation_seconds': prep_seconds, 'GPU_request_seconds': gpu_seconds,
                       'pipeline_seconds': time.perf_counter()-start, 'results': current,
                       'GPU_reports': reports})
            rows.extend(current)
        return rows

    def close(self):
        if self.gpu is not None:
            self.gpu.close()
        if self.preparer is not None:
            self.preparer.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
