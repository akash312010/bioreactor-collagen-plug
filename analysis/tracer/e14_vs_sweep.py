"""Fluent run E14 (channel G1, Pe = 69) against the independent-solver diffusivity sweep.

The scaled tracer problem depends only on Pe, so the sweep in G1 at D/D0 = 2 (dsweep_2.json) is equivalent to
Pe = 69 at the baseline diffusivity, the Peclet number of E14 (E8 at half its inlet velocity). The baseline
sweep (dsweep_1.json) corresponds to E8 itself. Replacement times are in flush times t_ideal.
Run:  python analysis/tracer/e14_vs_sweep.py
"""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
E = json.loads((HERE/'E_all.json').read_text())

def fluent(run):
    o = np.loadtxt(ROOT/'data'/'tracer_reports'/run/'outlet_massweighted.out', skiprows=3)
    v = np.loadtxt(ROOT/'data'/'tracer_reports'/run/'volume_average.out', skiprows=3)
    tid = E[run]['t_ideal_ms']*1e-3
    def t99(d):
        i = int(np.argmax(d[:, 1] >= 0.99))
        return float(np.interp(0.99, d[i-1:i+1, 1], d[i-1:i+1, 2])/tid)
    return t99(o), t99(v)

sweep = {f: json.loads((HERE/f'dsweep_{f}.json').read_text()) for f in ('1', '2', '0.5')}
rows = [('E8 (Fluent), Pe 138', *fluent('E8')),
        ('sweep D/D0=1 (Pe 138)', sweep['1']['mass_t99'], sweep['1']['vol_t99']),
        ('E14 (Fluent), Pe 69', *fluent('E14')),
        ('sweep D/D0=2 (Pe 69)', sweep['2']['mass_t99'], sweep['2']['vol_t99']),
        ('sweep D/D0=0.5 (Pe 276)', sweep['0.5']['mass_t99'], sweep['0.5']['vol_t99'])]
print(f"{'case':<26}{'t99,o':>8}{'t99,V':>8}{'ratio':>8}")
for name, o, v in rows:
    ratio = v/o if v == v else float('nan')
    print(f"{name:<26}{o:8.2f}{v:8.2f}{ratio:8.2f}" if v == v else f"{name:<26}{o:8.2f}{'>12':>8}{'>3.6':>8}")
