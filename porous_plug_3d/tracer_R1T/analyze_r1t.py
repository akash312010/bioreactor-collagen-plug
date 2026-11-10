"""R1T: 3D tracer on the frozen flow of the porous-plug run R1 (mesh M1, G3 permeability, 4.99 mL/min).

Reads the Fluent report series in reports/ (mass-weighted outlet concentration and plug volume average,
unit step at the inlet at t = 0) and compares the plug-scale replacement with the axisymmetric advective
calculation (analysis/porous/plug_transport_curves.json). Writes R1T_metrics.json and R1T_curves.npz.

Times are referred to the plug by subtracting the mean transit times of the inlet conduit (interior and
outlet) and of the outlet conduit (outlet only), and are scaled by tbar = phi*pi*R^2*delta/Q.
Run:  python porous_plug_3d/tracer_R1T/analyze_r1t.py
"""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REP = HERE/'reports'

def load(name):
    d = np.loadtxt(REP/name, skiprows=3)          # columns: time step, value, flow time
    return d[:, 2], d[:, 1]

def surface(name, key):
    for line in (REP/name).read_text().splitlines():
        f = line.split()
        if len(f) == 2 and f[0] == key:
            return float(f[1])
    raise KeyError(key)

def cross(x, c, lev):
    i = int(np.argmax(c >= lev))
    return float(np.interp(lev, c[i-1:i+1], x[i-1:i+1]))

def r1t_metrics():
    t, co = load('outlet_massweighted.out')
    t2, cv = load('plug_volume_average.out')
    assert np.allclose(t, t2)
    A_in, A_out = surface('area.txt', 'inlet'), surface('area.txt', 'outlet')
    u_in = surface('velocity_magnitude_area_weighted.txt', 'inlet')
    Q = A_in*u_in
    R, H, phi = 16.664e-3, 26.873e-3, 0.27
    tbar = phi*np.pi*R**2*H/Q
    t_in, t_out = A_in*31.68e-3/Q, A_out*1.777e-3/Q
    I = float(np.trapezoid(1-co, t))
    res = dict(Q_mL_min=Q*6e7, tbar_s=tbar, t_conduit_in_s=t_in, t_conduit_out_s=t_out,
               outlet_residual_integral_s=I, plug_residence_from_integral_s=I-t_in-t_out,
               plug_residence_rel_err=(I-t_in-t_out)/tbar-1)
    x_o = (t-t_in-t_out)/tbar
    x_v = (t-t_in)/tbar
    for lev in (0.5, 0.9, 0.99):
        res[f'out_t{int(lev*100)}'] = cross(x_o, co, lev)
        res[f'vol_t{int(lev*100)}'] = cross(x_v, cv, lev)
    res['ratio99'] = res['vol_t99']/res['out_t99']
    res['outlet_first_exceeds_1e-3'] = cross(x_o, co, 1e-3)
    for name, x, c in (('out', x_o, co), ('vol', x_v, cv)):
        y = 1-c
        m = (y < 0.03) & (y > 1e-3) & (x > 1)
        res[f'tail_time_constant_{name}'] = float(-1/np.polyfit(x[m], np.log(y[m]), 1)[0])
    # plug-scale solute balance, main-text Eq. (8): tbar dcV/dt = 1 - c_port_out once the inlet port is fresh
    dcv = np.gradient(cv, t)*tbar
    cpo = np.interp(t+t_out, t, co)
    m = t > 30
    res['balance_max_deviation_after_30s'] = float(np.max(np.abs(dcv[m]-(1-cpo[m]))))
    y = 1-cv
    mm = (y < 0.03) & (y > 1e-3)
    res['tail_residual_ratio_median'] = float(np.median(y[mm]/(1-cpo[mm])))
    ax = ROOT/'analysis'/'porous'/'plug_transport_curves.json'
    if ax.exists():
        A = json.loads(ax.read_text())['summary']
        res['axisymmetric_advective'] = dict(outlet=A['outlet'], interior=A['interior'], first_arrival=A['first_arrival'])
    return res, dict(x_o=x_o, x_v=x_v, co=co, cv=cv)

if __name__ == '__main__':
    res, curves = r1t_metrics()
    (HERE/'R1T_metrics.json').write_text(json.dumps(res, indent=1))
    np.savez(HERE/'R1T_curves.npz', **curves)
    for k, v in res.items():
        print(k, v if isinstance(v, dict) else round(v, 4))
