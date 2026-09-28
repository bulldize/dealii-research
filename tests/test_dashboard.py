import importlib.util
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch
from http.server import ThreadingHTTPServer
from urllib.request import urlopen, Request
from urllib.error import HTTPError

spec = importlib.util.spec_from_file_location('dashboard', Path(__file__).resolve().parents[1]/'scripts/dashboard.py')
d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

class DashboardTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
        self.cases = [dict(name='test', suite='convergence', steps=10)]
        self.dest = self.root/'convergence/test'
    def tearDown(self): self.tmp.cleanup()
    def state(self): return d.snapshot(self.root, self.cases)['cases'][0]['state']
    @patch.object(d, 'processes', return_value=[])
    def test_stopped_output_is_not_running(self, _):
        self.assertEqual(self.state(), 'pending'); self.dest.mkdir(parents=True)
        self.assertEqual(self.state(), 'incomplete')
        (self.dest/'status.json').write_text('{"status":"failed"}')
        self.assertEqual(self.state(), 'failed')
        (self.dest/'summary.json').write_text('{"completed":false}')
        self.assertEqual(self.state(), 'failed')
        (self.dest/'summary.json').write_text('{"completed":true}')
        self.assertEqual(self.state(), 'completed')
    def test_live_case_and_partial_csv(self):
        self.dest.mkdir(parents=True)
        with patch.object(d,'processes',return_value=[(123,['/app/giesekus','config',str(self.dest)])]):
            self.assertEqual(self.state(), 'running')
        history=self.dest/'history.csv';history.write_text('step,time\n1,0.01\n2,0.')
        self.assertEqual(d.last_row(history), {'step':'1','time':'0.01'})
    def test_auth_and_log_allowlist(self):
        (self.root/'pipeline.log').write_text('<script>untrusted log</script>')
        server=ThreadingHTTPServer(('127.0.0.1',0),d.make_handler(self.root,self.cases,'secret-token'))
        worker=threading.Thread(target=server.serve_forever,daemon=True);worker.start()
        base='http://127.0.0.1:'+str(server.server_port)
        try:
            with self.assertRaises(HTTPError) as err:urlopen(base+'/api/status')
            self.assertEqual(err.exception.code,401);err.exception.close()
            def fetch(path):return urlopen(Request(base+path,headers={'Authorization':'Bearer secret-token'}))
            with fetch('/api/status') as r:self.assertEqual(len(json.load(r)['cases']),1)
            with fetch('/api/log?name=pipeline') as r:self.assertIn(b'<script>',r.read())
            with self.assertRaises(HTTPError) as err:fetch('/api/log?name=../../etc/passwd')
            self.assertEqual(err.exception.code,404);err.exception.close()
        finally:server.shutdown();server.server_close();worker.join()

if __name__=='__main__':unittest.main()
