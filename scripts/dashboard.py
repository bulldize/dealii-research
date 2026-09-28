#!/usr/bin/env python3
"""Read-only paper experiment dashboard. Python 3 standard library only."""
import argparse
import csv
import hmac
import json
import os
from pathlib import Path
import secrets
import shutil
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlsplit

REPO = Path(__file__).resolve().parents[1]
LOGS = ('pipeline', 'dealii-configure', 'dealii-build', 'dealii-install',
        'research-build', 'calibration', 'convergence', 'contraction', 'plotting')


def read_json(path, default=None):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def tail(path, size=12000):
    try:
        with path.open('rb') as f:
            f.seek(0, 2)
            f.seek(max(0, f.tell() - size))
            return f.read().decode('utf-8', errors='replace')
    except OSError:
        return ''


def last_row(path):
    try:
        with path.open() as f:
            header = next(csv.reader([f.readline()]))
        # Ignore any last line still being written by the simulation.
        lines = tail(path, 8192).split('\n')[:-1]
        for line in reversed(lines):
            row = next(csv.reader([line]))
            if len(row) == len(header) and row[0].isdigit():
                return dict(zip(header, row))
    except (OSError, ValueError, StopIteration):
        pass
    return {}


def processes():
    result = []
    for path in Path('/proc').glob('[0-9]*/cmdline'):
        try:
            args = [s.decode(errors='replace') for s in path.read_bytes().split(b'\0') if s]
            if args:
                result.append((int(path.parent.name), args))
        except OSError:
            pass
    return result


def snapshot(root, manifest):
    procs = processes()
    cases = []
    for case in manifest:
        dest = root / case['suite'] / case['name']
        summary = read_json(dest / 'summary.json', {})
        status = read_json(dest / 'status.json', {})
        active = any(Path(args[0]).name == 'giesekus' and str(dest) in args for _, args in procs)
        state = ('completed' if summary.get('completed') is True else
                 'running' if active else 'failed' if status.get('status') == 'failed' else
                 'incomplete' if dest.exists() else 'pending')
        history = last_row(dest / 'history.csv')
        cases.append(dict(name=case['name'], suite=case['suite'], state=state,
                          steps=case['steps'], history=history))
    try:
        phase = (root / 'phase.txt').read_text().strip()
        updated = (root / 'phase.txt').stat().st_mtime
    except OSError:
        phase, updated = '尚无阶段记录', None
    pipeline_alive = False
    try:
        pid = int((root / 'pipeline.pid').read_text())
        pipeline_alive = any(p == pid and any(Path(a).name == 'server_pipeline.sh' for a in args)
                             for p, args in procs)
    except (OSError, ValueError):
        pass
    return dict(time=time.time(), phase=phase, phase_updated=updated,
                pipeline_alive=pipeline_alive, process_detection=Path('/proc/self').exists(),
                cases=cases, calibration=read_json(root / 'calibration/settings.json', {}))


class Resources:
    def __init__(self):
        self.previous = None
        self.cpu = None
        self.lock = threading.Lock()

    def sample(self, root):
        with self.lock:
            try:
                values = list(map(int, Path('/proc/stat').read_text().splitlines()[0].split()[1:9]))
                total, idle = sum(values), values[3] + values[4]
                if self.previous and total > self.previous[0]:
                    self.cpu = round(100 * (1 - (idle-self.previous[1])/(total-self.previous[0])), 1)
                self.previous = total, idle
            except (OSError, ValueError, IndexError):
                pass
            memory = {}
            try:
                fields = {r.split(':')[0]: int(r.split()[1])*1024
                          for r in Path('/proc/meminfo').read_text().splitlines()}
                memory = dict(total=fields['MemTotal'], used=fields['MemTotal']-fields['MemAvailable'],
                              swap_used=fields['SwapTotal']-fields['SwapFree'])
            except (OSError, ValueError, KeyError):
                pass
            disk = shutil.disk_usage(root if root.exists() else REPO)
            return dict(cpu_percent=self.cpu, logical_cpus=os.cpu_count(), memory=memory,
                        disk_free=disk.free, gpu='本程序不使用 GPU')


def make_handler(root, manifest, token):
    resources = Resources()
    names = {c['name']: c for c in manifest}

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *_):
            pass  # Do not write tokens from query strings into access logs.

        def do_GET(self):
            url = urlsplit(self.path)
            query = parse_qs(url.query)
            supplied = self.headers.get('Authorization', '').removeprefix('Bearer ') or query.get('token', [''])[0]
            if token and not hmac.compare_digest(supplied, token):
                self.send_data(401, b'Authentication required', 'text/plain'); return
            if url.path == '/':
                self.send_data(200, (REPO / 'web/dashboard.html').read_bytes(), 'text/html'); return
            if url.path == '/api/status':
                result = snapshot(root, manifest)
                result['resources'] = resources.sample(root)
                self.send_data(200, json.dumps(result, allow_nan=False).encode(), 'application/json'); return
            if url.path == '/api/log':
                name = query.get('name', ['pipeline'])[0]
                if name in LOGS:
                    path = root / (name + '.log')
                elif name in names:
                    case = names[name]; path = root / case['suite'] / name / 'run.log'
                else:
                    self.send_data(404, b'Unknown log', 'text/plain'); return
                self.send_data(200, tail(path, 24000).encode(), 'text/plain'); return
            self.send_data(404, b'Not found', 'text/plain')

        def send_data(self, status, data, kind):
            self.send_response(status)
            self.send_header('Content-Type', kind+'; charset=utf-8')
            self.send_header('Content-Length', str(len(data)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.send_header('Referrer-Policy', 'no-referrer')
            self.end_headers()
            self.wfile.write(data)
    return Handler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=REPO/'results/full-paper')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', type=int, default=8766)
    parser.add_argument('--token-file', type=Path)
    args = parser.parse_args()
    token = ''
    if args.host not in ('127.0.0.1', 'localhost') and not args.token_file:
        parser.error('Binding to the network requires --token-file (created automatically if absent)')
    if args.token_file:
        if not args.token_file.exists():
            args.token_file.parent.mkdir(parents=True, exist_ok=True)
            fd = os.open(args.token_file, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            with os.fdopen(fd, 'w') as f:
                f.write(secrets.token_urlsafe(32)+'\n')
        token = args.token_file.read_text().strip()
        if len(token) < 24:
            parser.error('Token must contain at least 24 characters')
    manifest = read_json(REPO/'configs/paper/manifest.json')
    if not manifest:
        parser.error('Missing or invalid experiment manifest')
    server = ThreadingHTTPServer((args.host, args.port), make_handler(args.root.resolve(), manifest, token))
    print(f'Dashboard listening on {args.host}:{args.port}; root={args.root.resolve()}', flush=True)
    server.serve_forever()


if __name__ == '__main__':
    main()
