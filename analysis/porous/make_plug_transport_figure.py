"""Figure: flow distribution and replacement at the plug scale (port-fed Darcy plug).
(a) Fraction of the plug volume whose superficial speed exceeds a multiple of the uniform value Q/A.
(b) Unreplaced fraction at the outlet and in the interior against time in mean residence times.
(b) also shows the 3D Fluent tracer run R1T (symbols) when porous_plug_3d/tracer_R1T/R1T_curves.npz
(from porous_plug_3d/tracer_R1T/analyze_r1t.py) is present.
Run:  python analysis/porous/make_plug_transport_figure.py [h_mm]
      -> figures/Figure_Plug_Transport_V2.pdf/.png and analysis/porous/plug_transport_curves.json
Panels (c-e) of the paper's Fig. 7 are added by porous_plug_3d/tracer_R1T/compose_figure7.py.
Arial is used when installed; otherwise the metric-identical Liberation Sans."""
import sys, json, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plug_transport import summary
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
OUT = ROOT/'figures'; OUT.mkdir(exist_ok=True)
R1T = ROOT/'porous_plug_3d'/'tracer_R1T'/'R1T_curves.npz'

CAND = [Path('C:/Windows/Fonts/arial.ttf'), Path('/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'),
        Path('/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf')]
FP = next(p for p in CAND if p.is_file()); font_manager.fontManager.addfont(str(FP))
FONT = font_manager.FontProperties(fname=str(FP)).get_name()
plt.rcParams.update({'font.family': FONT, 'font.size': 12, 'axes.labelsize': 12, 'xtick.labelsize': 12,
    'ytick.labelsize': 12, 'legend.fontsize': 12, 'mathtext.fontset': 'custom', 'mathtext.rm': FONT,
    'mathtext.it': FONT+':italic', 'axes.linewidth': .7, 'lines.linewidth': 1.3, 'xtick.direction': 'in',
    'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True, 'figure.facecolor': 'white',
    'savefig.facecolor': 'white', 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
ORANGE, BLACK, GREY = '#bc5618', '#171717', '#777777'
CM = 1/2.54

h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.025
res, f, (e, dpsi, T, t, co, cv) = summary(h)
sp = f['speed'].ravel(); w = (f['vol']/f['vol'].sum()).ravel()
o = np.argsort(sp); s_sorted = sp[o]; above = 1-np.cumsum(w[o])
s = np.geomspace(0.05, 60, 300); frac = np.interp(s, s_sorted, above)

fig, ax = plt.subplots(1, 2, figsize=(17*CM, 7.6*CM), layout='constrained')
ax[0].loglog(s, frac, color=BLACK)
ax[0].axvline(1, color=GREY, lw=.8, ls='--')
ax[0].set(xlabel=r'Superficial speed / $(Q/A)$', ylabel='Plug volume fraction\nabove this speed', xlim=(0.05, 60), ylim=(1e-4, 1.5))
ax[1].semilogy(t, 1-co, color=BLACK, label='Outlet')
ax[1].semilogy(t, 1-cv, '--', color=ORANGE, label='Interior')
tu = np.linspace(0, 1, 200)
ax[1].semilogy(tu, 1-tu, ':', color=GREY, label='Uniform feeding')
if R1T.exists():
    F = np.load(R1T)
    for x, c, mk, col, lab in ((F['x_o'], F['co'], 'o', BLACK, '3D, outlet'), (F['x_v'], F['cv'], 's', ORANGE, '3D, interior')):
        xs = np.arange(0.25, 4.5, 0.25); ys = np.interp(xs, x, 1-c)
        ax[1].semilogy(xs, ys, mk, mec=col, mfc='white', ms=4.5, mew=1.0, ls='none', label=lab)
ax[1].axhline(0.01, color=GREY, lw=.8, ls='--')
ax[1].set(xlabel=r'Time / mean residence time, $t/\bar{t}$', ylabel='Unreplaced fraction', xlim=(0, 4.5), ylim=(3e-3, 1.3))
from matplotlib.legend_handler import HandlerTuple
hl, lb = ax[1].get_legend_handles_labels(); d = dict(zip(lb, hl))
if '3D, outlet' in d:     # lines: axisymmetric advective calculation; symbols: 3D simulation
    hs = [(d['Outlet'], d['3D, outlet']), (d['Interior'], d['3D, interior']), d['Uniform feeding']]
    ax[1].legend(hs, ['Outlet', 'Interior', 'Uniform feeding'], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)},
                 frameon=False, loc='upper right', handlelength=2.6, borderaxespad=0.2, labelspacing=0.3)
else:
    ax[1].legend(frameon=False, loc='upper right', handlelength=1.6)
for a, L in zip(ax, ('(a)', '(b)')):
    a.text(-0.25, 1.04, L, transform=a.transAxes, ha='left', va='bottom', fontsize=12)
for ext in ('pdf', 'png'): fig.savefig(OUT/f'Figure_Plug_Transport_V2.{ext}', dpi=300)
json.dump(dict(summary=res, speed_multiple=s.tolist(), volume_fraction_above=frac.tolist(),
               t_over_tbar=t[::10].tolist(), outlet_unreplaced=(1-co)[::10].tolist(), interior_unreplaced=(1-cv)[::10].tolist()),
          open(HERE/'plug_transport_curves.json', 'w'), indent=1)
for x in (1.18, 7.1): print('fraction above', x, float(np.interp(x, s_sorted, above)))
print(json.dumps(res['outlet']), json.dumps(res['interior']), res['first_arrival'])
