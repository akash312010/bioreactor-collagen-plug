"""Single-column (8.2 cm) charts at final placed Arial12 for the manuscript.

Uses only existing results (headline_numbers.json, device CSV, curvi_test.json,
fcurve.json). Panel letters only; no panel titles. Black/burnt-orange palette.
Arial is used when installed (author's Windows machine); otherwise the
metric-identical Liberation Sans is used for this preview build.
Run:  python make_singlecol_figures.py <ROOT>   (ROOT = Paper_Final folder)
"""
import sys, json, csv
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path(sys.argv[1]); OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT/'figures'
OUT.mkdir(parents=True, exist_ok=True)
CM = 1/2.54; W = 8.2*CM
from figure_style import FAMILY as FONT
from run_labels import PAPER_LABEL
plt.rcParams.update({'font.family': FONT, 'font.size': 12, 'axes.labelsize': 12,
    'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 12,
    'mathtext.fontset': 'custom', 'mathtext.rm': FONT, 'mathtext.it': FONT+':italic',
    'axes.linewidth': .7, 'lines.linewidth': 1.2, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True, 'figure.facecolor': 'white',
    'savefig.facecolor': 'white', 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
ORANGE, BLACK, GREY = '#bc5618', '#171717', '#777777'
J = lambda p: json.loads((ROOT/p).read_text())
H = J('headline_numbers.json')

def save(fig, name):
    for ext in ('pdf', 'svg'): fig.savefig(OUT/f'{name}.{ext}')
    fig.savefig(OUT/f'{name}.png', dpi=250); plt.close(fig)
def letter(ax, s, x=-0.2): ax.text(x, 1.02, s, transform=ax.transAxes, ha='left', va='bottom', fontsize=12)

# ---- Fig. 3: pressure budget, device loss vs plug demand -------------------
rows = list(csv.DictReader(open(ROOT/'data/device/Part1_data_reconciled.csv')))
Qd = sorted({float(r['Q_ml_min']) for r in rows})
dev = {m: [float(r['dP_Pa']) for r in rows if r['mesh'] == m] for m in ('0.5M', '1M', '2M')}
lo = np.min([dev[m] for m in dev], axis=0); hi = np.max([dev[m] for m in dev], axis=0)
cap = H['capacity']; F = cap['F_finest']
s = np.array(list(cap['slope_uniform_mLmin_per_kPa'].values())); s3 = cap['slope_uniform_mLmin_per_kPa']['G3']
Q = np.linspace(5, 30, 60)
fig, ax = plt.subplots(figsize=(W, 7*CM), layout='constrained')
for fac, col, ls, lab in [(F, ORANGE, '--', 'Port-fed plug'), (1, BLACK, '-', 'Uniformly fed plug')]:
    ax.fill_between(Q, 1e3*Q*fac/s.max(), 1e3*Q*fac/s.min(), color=col, alpha=.12, lw=0)
    ax.semilogy(Q, 1e3*Q*fac/s3, ls, color=col, label=lab)
ax.errorbar(Qd, dev['1M'], yerr=[np.array(dev['1M'])-lo, hi-np.array(dev['1M'])], fmt='o',
            color=BLACK, mfc='white', ms=4, capsize=2, lw=.8, label='Open device')
ax.set(xlabel=r'Flow rate, $Q$ (mL min$^{-1}$)', ylabel=r'Pressure difference (Pa)', xlim=(3, 32), ylim=(.2, 5e5))
ax.legend(frameon=False, loc='center right', bbox_to_anchor=(1, .42))
save(fig, 'Fig_Pressure_Budget')

# ---- Fig. 5: theory vs CFD, reference value + range in one panel ------------
C = J('analysis/theory/curvi_test.json'); G = [f'G{i}' for i in range(1, 6)]; x = np.arange(5)
cuts = ['4_8', '6_12', '8_16', '10_20']
def val(g, sig, cut, key):
    d = C[g][f'sig{sig}'][cut]
    return 1/d['cfd_over_lub'] if key == 'leading' else d[key]
fig, ax = plt.subplots(figsize=(W, 9.0*CM), layout='constrained')
for key, lab, col, mk, sh, filled in [('leading', 'Leading order', GREY, 'o', -.2, False),
                                       ('no_curv_over_cfd', 'Offset/curvature omitted', BLACK, 's', 0, False),
                                       ('theory_over_cfd', 'Full second order', ORANGE, 'D', .2, True)]:
    ref = [val(g, 4, '6_12', key) for g in G]
    rng = np.array([[val(g, sg, c, key) for sg in (3, 4) for c in cuts] for g in G])
    ax.vlines(x+sh, rng.min(1), rng.max(1), color=col, lw=1.2)
    ax.plot(x+sh, ref, mk, color=col, mfc=col if filled else 'white', ms=5, label=lab)
ax.axhline(1, color=GREY, lw=.8, ls='--'); ax.set_xticks(x, G)
ax.set(ylabel=r'$\Delta p_{\rm model}/\Delta p_{\rm CFD}$', xlabel='Channel', ylim=(.62, 1.32))
ax.legend(frameon=False, loc='lower left', bbox_to_anchor=(-.27, 1.01), fontsize=12, borderaxespad=0)
save(fig, 'Fig_Theory_vs_CFD')

# ---- Fig. 6: capacity and spreading factor, stacked --------------------------
fig, axs = plt.subplots(2, 1, figsize=(W, 15*CM), layout='constrained')
P = np.geomspace(.5, 100, 200)
for fac, col, ls, lab in [(1, BLACK, '-', 'Uniform feeding'), (F, ORANGE, '--', 'Port feeding')]:
    axs[0].fill_between(P, P*s.min()/fac, P*s.max()/fac, color=col, alpha=.12, lw=0)
    axs[0].loglog(P, P*s3/fac, ls, color=col, label=lab)
axs[0].set(xlabel='Pressure across plug,\n'+r'$\Delta P_c$ (kPa)', ylabel=r'$Q_{\max}$ (mL min$^{-1}$)', xlim=(.5, 100), ylim=(.03, 300))
axs[0].legend(frameon=False, loc='upper left'); letter(axs[0], '(a)')
fc = J('analysis/porous/fcurve.json'); a = np.geomspace(.03, 1, 250); R, Ht = 16.664, 26.873
axs[1].loglog(a, 1+np.pi*R/(2*Ht)*(1-a)**1.5/a, color=ORANGE, label='Correlation')
axs[1].loglog(fc['a_over_R'], fc['F'], 'o', mfc='white', mec=BLACK, ms=4.5, label='Darcy calculation')
axs[1].plot((fc['device']['a1']+fc['device']['a2'])/(2*R), F, '*', ms=10, color=ORANGE, mec=BLACK, mew=.5, label='Device ports')
axs[1].set(xlabel=r'Port-to-plug radius ratio, $a/R$', ylabel=r'Spreading factor, $F$', xlim=(.025, 1.1), ylim=(.9, 100))
axs[1].legend(frameon=False, loc='upper right', markerfirst=False); letter(axs[1], '(b)')
save(fig, 'Fig_Plug_Capacity')

# ---- Fig. 7: parity plot, outlet vs internal t99 ----------------------------
tr = H['tracer']['runs']
sty = {'G3': ('o', BLACK, 'white'), 'G5': ('s', BLACK, BLACK), 'G1': ('D', ORANGE, ORANGE)}
fig, ax = plt.subplots(figsize=(W, 8.7*CM), layout='constrained')
lim = (1.2, 10)
ax.plot(lim, lim, '--', color=GREY, lw=.8)
done = set()
for r, d in tr.items():
    mk, ec, fc_ = sty[d['geom']]
    ax.loglog(d['t99_out'], d['t99_V'], mk, mec=ec, mfc=fc_, ms=5.5,
              label=d['geom'] if d['geom'] not in done else None); done.add(d['geom'])
for r, (tx, ty) in {'E8': (1.86, 8.2), 'E14': (1.72, 4.9), 'E7': (1.85, 2.85), 'E10': (3.25, 4.6)}.items():
    ax.text(tx, ty, PAPER_LABEL[r], fontsize=12)
ax.text(1.3, 6.6, 'Interior lags outlet', fontsize=12, color=GREY)
ax.text(4.2, 1.5, 'Interior first', fontsize=12, color=GREY)
ax.set(xlim=lim, ylim=lim, xlabel=r'Outlet $t_{99,o}/t_{\rm ideal}$', ylabel=r'Internal $t_{99,V}/t_{\rm ideal}$')
ax.set_aspect('equal',adjustable='box')
for axis in (ax.xaxis, ax.yaxis):
    axis.set_major_formatter(matplotlib.ticker.ScalarFormatter()); axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks([2, 3, 5, 10]); ax.set_yticks([2, 3, 5, 10])
handles = [ax.plot([], [], m, mec=e, mfc=f, ms=5.5, ls='none')[0] for g, (m, e, f) in sty.items()]
ax.legend(handles, ['G3', 'G5', 'G1'], frameon=False, loc='lower center', bbox_to_anchor=(.5,1.02), ncol=3, columnspacing=.6, handletextpad=.4, handlelength=.8, borderaxespad=0)
save(fig, 'Fig_Replacement_Parity')
print('font:', FONT, '-> wrote 4 figures to', OUT)
