"""Author's final-size figure specification (Arial 12 pt at 17 cm).
Uses Arial when installed (author's Windows machine); otherwise the metric-identical
Liberation Sans so that layout is unchanged."""
from pathlib import Path
import matplotlib as mpl
from matplotlib import font_manager
CM = 1 / 2.54
WIDTH = 17 * CM
CANDIDATES = [Path('C:/Windows/Fonts/arial.ttf'), Path('/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'),
              Path('/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf')]
FONT = next(p for p in CANDIDATES if p.is_file())
font_manager.fontManager.addfont(str(FONT))
FAMILY = font_manager.FontProperties(fname=str(FONT)).get_name()
STYLE = {
    'font.family': FAMILY, 'font.size': 12,
    'axes.labelsize': 12, 'axes.titlesize': 12,
    'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 12,
    'mathtext.fontset': 'custom', 'mathtext.rm': FAMILY,
    'mathtext.it': FAMILY + ':italic', 'mathtext.bf': FAMILY + ':bold',
    'axes.linewidth': .7, 'lines.linewidth': 1.2,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
}
mpl.rcParams.update(STYLE)
