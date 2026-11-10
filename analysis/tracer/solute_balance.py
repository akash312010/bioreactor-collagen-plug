"""Solute balance of the tracer simulations (main-text Eq. 8, Supplemental Table S5).

For each run in data/tracer_reports/ it checks  t_ideal dc_V/dt = 1 - c_o  (after half a flush time),
fits the late decay time of the interior residual r = 1 - c_V between residuals 0.001 and 0.03
(exchange number Lambda_e = tau_e/t_ideal), and takes the median ratio r/(1 - c_o) in the same window,
which the balance predicts to equal Lambda_e. Also prints t99,V/t99,o.
Run:  python analysis/tracer/solute_balance.py      -> solute_balance.json
"""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
E = json.loads((HERE/'E_all.json').read_text())
RUNS = ['E11', 'E12', 'E1', 'E2', 'E3', 'E4', 'E6', 'E7', 'E14', 'E8', 'E9', 'E10']
GEOM = {'p28': 'G1', 'p35': 'G3', 'p45': 'G5'}

def load(run, name):
    d = np.loadtxt(ROOT/'data'/'tracer_reports'/run/name, skiprows=3)   # time step, value, flow time
    return d[:, 2], d[:, 1]

def cross(x, c, lev):
    i = int(np.argmax(c >= lev))
    if c[i] < lev:
        return float('nan')
    return float(np.interp(lev, c[i-1:i+1], x[i-1:i+1]))

out = {}
for run in RUNS:
    t, co = load(run, 'outlet_massweighted.out')
    t2, cv = load(run, 'volume_average.out')
    assert np.allclose(t, t2), run
    tid = E[run]['t_ideal_ms']*1e-3
    x = t/tid
    y = 1-cv
    m = (y > 1e-3) & (y < 3e-2) & (t > t[np.argmax(y < 0.5)])
    tau = -1/np.polyfit(t[m], np.log(y[m]), 1)[0]
    dev = np.abs(np.gradient(cv, t)*tid-(1-co))[x > 0.5].max()
    out[run] = dict(channel=GEOM[E[run]['geom']], Pe=E[run]['Pe_w'], t_ideal_ms=tid*1e3,
                    Lambda_e=tau/tid, residual_ratio_median=float(np.median(y[m]/(1-co[m]))),
                    balance_max_deviation=float(dev), t99_out=cross(x, co, .99), t99_V=cross(x, cv, .99))
    out[run]['ratio99'] = out[run]['t99_V']/out[run]['t99_out']
(HERE/'solute_balance.json').write_text(json.dumps(out, indent=1))
print('Run  Ch.  Lambda_e  r/(1-c_o)  t99V/t99o  max|balance|')
for r, v in out.items():
    print(f"{r:<4} {v['channel']}   {v['Lambda_e']:6.2f}    {v['residual_ratio_median']:6.2f}     {v['ratio99']:5.2f}     {v['balance_max_deviation']:.1e}")
