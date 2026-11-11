"""Interior replacement reconstructed from the outlet history alone (main-text Eq. 19, Supplemental Table S6).

Solute conservation gives the interior residual from the outlet record,
    r(t) = (1/t_ideal) * integral_t^infinity [1 - c_o(t')] dt'.
The integral over the recorded interval is evaluated by the trapezoidal rule. Beyond the end of the record
the outlet residual is continued as a single exponential fitted to its last factor-of-three decrease, which
adds the area r_o(t_end)*tau_o. The reconstructed 99% time is compared with the directly computed one, with
and without that continuation. The script also reports the ratio of interior to outlet replacement time at
90%, 95% and 99%, which shows that the interior lag of the pocket runs is a late-time effect.
Run:  python analysis/tracer/outlet_reconstruction.py      -> outlet_reconstruction.json
"""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_labels import PAPER_LABEL

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
    _, cv = load(run, 'volume_average.out')
    tid = E[run]['t_ideal_ms']*1e-3
    x = t/tid
    ro = 1-co
    seg = 0.5*(ro[1:]+ro[:-1])*np.diff(t)
    recorded = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])       # integral from t to the end of the record
    m = (ro <= 3*ro[-1]) & (ro > 0)
    tau_o = -1/np.polyfit(t[m], np.log(ro[m]), 1)[0] if (m.sum() > 5 and ro[-1] > 1e-9) else 0.0
    tail = ro[-1]*tau_o
    direct = cross(x, cv, 0.99)
    rec = cross(x, 1-(recorded+tail)/tid, 0.99)
    rec_no_tail = cross(x, 1-recorded/tid, 0.99)
    out[run] = dict(paper_label=PAPER_LABEL[run], channel=GEOM[E[run]['geom']], t99_V_direct=direct, t99_V_from_outlet=rec,
                    rel_difference=rec/direct-1, t99_V_from_outlet_no_continuation=rec_no_tail,
                    rel_difference_no_continuation=rec_no_tail/direct-1,
                    record_length_flush_times=float(x[-1]), outlet_residual_at_end=float(ro[-1]),
                    continuation_area_over_t_ideal=float(tail/tid),
                    ratio_V_over_o={str(lev): cross(x, cv, lev)/cross(x, co, lev) for lev in (0.9, 0.95, 0.99)})
(HERE/'outlet_reconstruction.json').write_text(json.dumps(out, indent=1))
print('Paper (Fluent run)  Ch.  t99V direct  from outlet   diff    no continuation   tV/to at 90%  95%   99%')
for r, v in out.items():
    q = v['ratio_V_over_o']
    print(f"{v['paper_label']:<4}({r:<4})       {v['channel']}   {v['t99_V_direct']:8.2f}   {v['t99_V_from_outlet']:8.2f}   {100*v['rel_difference']:+5.1f}%"
          f"   {100*v['rel_difference_no_continuation']:+6.1f}%          {q['0.9']:.2f}  {q['0.95']:.2f}  {q['0.99']:.2f}")
