"""Scheduler and protocol contracts; these tests do not validate model economics."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from cuda_pipeline import StreamProcess
from gpu import run_search


class FakeEngine:
    model='RE'
    def __init__(self):
        self.points=[]
    def evaluate(self, points, sigma_bounds=None):
        rows=[]
        for x in points:
            self.points.append(np.array(x))
            parameters=[float(v) for v in x]+[0.25]
            rows.append({'objective':float(np.sum(x*x)), 'parameters':parameters,
                         'parameter_hex':[v.hex() for v in parameters],'status':'ok'})
        return rows


class SearchContracts(unittest.TestCase):
    def spec(self,budget):
        initial=np.random.default_rng(3).uniform(.1,.8,(5,8))
        return {'model':'RE','bounds':[[.1,.8]]*8,'sigma_pi_bounds':[0,1],
                'seed':9,'max_evaluations':budget,'initial_population':initial.tolist(),
                'tol':0,'atol':0}
    def test_batch_matches_serial_deferred_trajectory(self):
        spec=self.spec(15);engine=FakeEngine()
        with tempfile.TemporaryDirectory() as d:
            result=run_search(engine,spec,Path(d)/'search.json')
        serial=[]
        def objective(x):
            serial.append(x.copy());return float(np.sum(x*x))
        reference=differential_evolution(objective,spec['bounds'],seed=spec['seed'],
            init=np.array(spec['initial_population']),maxiter=2,polish=False,
            updating='deferred',tol=0,atol=0)
        np.testing.assert_array_equal(engine.points,serial)
        self.assertEqual(result['completed_evaluations'],15)
        self.assertEqual(result['best']['objective'],float(reference.fun))
    def test_budget_stops_at_exact_prefix(self):
        engine=FakeEngine()
        with tempfile.TemporaryDirectory() as d:
            result=run_search(engine,self.spec(7),Path(d)/'search.json')
        self.assertEqual(len(engine.points),7)
        self.assertEqual(result['attempted_evaluations'],7)
        self.assertEqual(result['completed_evaluations'],7)
        self.assertEqual(result['status'],'evaluation_budget_reached')
        self.assertEqual(result['generations_completed'],0)
    def test_failed_batch_is_not_reported_as_completed(self):
        class Broken(FakeEngine):
            def evaluate(self,*args,**kwargs):raise RuntimeError('worker failure')
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'search.json'
            with self.assertRaisesRegex(RuntimeError,'worker failure'):
                run_search(Broken(),self.spec(10),path)
            result=json.loads(path.read_text())
        self.assertEqual(result['status'],'failed')
        self.assertEqual(result['completed_evaluations'],0)
        self.assertEqual(result['attempted_evaluations'],5)
        self.assertIsNone(result['best'])


class ProtocolContracts(unittest.TestCase):
    def test_multiple_requests_share_one_process(self):
        source='''import sys,json
print(json.dumps({'family':'ready'}),flush=True)
for i,line in enumerate(sys.stdin):
 if line.strip()=='QUIT':break
 print(json.dumps({'family':'evaluation','sequence':i}),flush=True)
 print(json.dumps({'family':'batch_complete','status':0}),flush=True)
'''
        with tempfile.TemporaryDirectory() as d:
            worker=StreamProcess([sys.executable,'-u','-c',source],d,Path(d)/'worker.log',2)
            try:
                self.assertEqual(worker.request('a')[0]['sequence'],0)
                self.assertEqual(worker.request('b')[0]['sequence'],1)
            finally:worker.close()
    def test_unexpected_exit_does_not_hang(self):
        with tempfile.TemporaryDirectory() as d:
            worker=StreamProcess([sys.executable,'-u','-c',"print('{\"family\":\"ready\"}',flush=True)"],d,Path(d)/'worker.log',2)
            try:
                with self.assertRaises((RuntimeError,BrokenPipeError)):
                    worker.request('a')
            finally:worker.close()


if __name__=='__main__':unittest.main()
