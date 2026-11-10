"""Author's final-size figure standard: Arial 12 pt, black text, vector PDFs."""
from pathlib import Path
import json
import os

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
from matplotlib.text import Text
from matplotlib.image import AxesImage, FigureImage

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "Manuscript/figures/2026-10-21_arial12"
WIDTH_CM = 17.0
WIDTH_IN = WIDTH_CM / 2.54
FONT_SIZE = 12

# Windows Arial is installed on the author's machine. Other machines must also
# supply Arial; do not silently deliver a different font under this style name.
for filename in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
    path = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts" / filename
    if path.exists():
        font_manager.fontManager.addfont(str(path))
font_manager.findfont("Arial", fallback_to_default=False)

STYLE = {
    "font.family": "Arial", "font.size": FONT_SIZE,
    "axes.labelsize": FONT_SIZE, "axes.titlesize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE, "figure.titlesize": FONT_SIZE,
    "text.color": "black", "axes.labelcolor": "black",
    "axes.titlecolor": "black", "xtick.color": "black", "ytick.color": "black",
    "axes.edgecolor": "black", "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 600,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial", "mathtext.default": "it",
    "axes.spines.top": False, "axes.spines.right": False,
}
matplotlib.rcParams.update(STYLE)


def enforce_style(fig):
    """Override legacy per-annotation font/colour settings before export."""
    for artist in fig.findobj(match=Text):
        artist.set_color("black")
        artist.set_fontfamily("Arial")
        artist.set_fontsize(FONT_SIZE)
    for ax in fig.axes:
        ax.tick_params(axis="both", which="both", colors="black", labelsize=FONT_SIZE)
        for spine in ax.spines.values():
            spine.set_color("black")
        # Matplotlib rasterizes long colour bars by default. Preserve paths.
        for collection in ax.collections:
            collection.set_rasterized(False)
    fig.canvas.draw()


def save_figure(fig, stem, outdir=None):
    outdir = Path(outdir or os.environ.get("OUTDIR", DEFAULT_OUT))
    outdir.mkdir(parents=True, exist_ok=True)
    enforce_style(fig)
    images = fig.findobj(match=lambda item: isinstance(item, (AxesImage, FigureImage)))
    mode = "hybrid: raster image panels with vector text/paths" if images else "vector"
    fig.savefig(outdir / (stem + ".pdf"), dpi=600)
    fig.savefig(outdir / (stem + ".svg"), dpi=600)
    fig.savefig(outdir / (stem + ".png"), dpi=220)
    metadata = {
        "font": "Arial", "font_size_pt": FONT_SIZE, "text_colour": "black",
        "width_cm": round(fig.get_figwidth() * 2.54, 3),
        "height_cm": round(fig.get_figheight() * 2.54, 3),
        "content": mode, "raster_image_artists": len(images),
        "files": [stem + suffix for suffix in (".pdf", ".svg", ".png")],
    }
    (outdir / (stem + ".export.json")).write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Saved {outdir / (stem + '.pdf')} ({mode})")
    if images:
        print("  Hybrid preview only: use exported geometry for the vector field edition.")
