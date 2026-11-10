"""Compose the ParaView panels of the R1T tracer into one vector figure (3 x 2 panels + colour bar), panel letters only.
Run: python porous_plug_3d/tracer_R1T/compose_r1t_planes.py   (after render_r1t_planes_paraview.py)
Writes figures/Supplementary_Figure_9.pdf."""
import sys
from pathlib import Path
import pymupdf as fitz

ROOT = Path(__file__).resolve().parents[2]
SRC, OUT = ROOT/'figures'/'r1t_panels', str(ROOT/'figures'/'Supplementary_Figure_9.pdf')
FONT = next(p for p in ['C:/Windows/Fonts/arial.ttf', '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf',
                        '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'] if Path(p).is_file())

FONTI = next(p for p in ['C:/Windows/Fonts/ariali.ttf', '/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf',
                         '/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf'] if Path(p).is_file())
def tlabel(page, x, y, sec):
    """Write 't = N s' with an italic t at (x, y)."""
    page.insert_text((x, y), 't', fontname='ArialItalicT', fontfile=FONTI, fontsize=12)
    w = fitz.Font(fontfile=FONTI).text_length('t', fontsize=12)
    page.insert_text((x + w, y), f' = {sec} s', fontname='ArialRegularT', fontfile=FONT, fontsize=12)
cm = 72/2.54
W = 17*cm; pw = 5.2*cm; ph = pw*760/700; gx = (W-3*pw)/2; top = 0.7*cm; gy = 0.8*cm
cw = 0.8*W; ch = cw*110/500
H = top + 2*ph + gy + 0.1*cm + ch
d = fitz.open(); p = d.new_page(width=W, height=H)
times = ['025.0s', '050.0s', '075.0s', '150.0s', '250.0s', '400.0s']
for k, t in enumerate(times):
    x = (k % 3)*(pw+gx); y = top + (k//3)*(ph+gy)
    src = fitz.open(SRC/f'R1T_plane_t{t}.pdf')
    p.show_pdf_page(fitz.Rect(x, y, x+pw, y+ph), src, 0)
    p.insert_text((x+0.05*cm, y-0.15*cm), f'({"abcdef"[k]})', fontname='ArialRegularT', fontfile=FONT, fontsize=12)
    tlabel(p, x + 0.1*cm, y + 0.45*cm, int(float(t[:-1])))
cx = 0.1*W; cy = top + 2*ph + gy + 0.1*cm
p.show_pdf_page(fitz.Rect(cx, cy, cx+cw, cy+ch), fitz.open(SRC/'R1T_colorbar.pdf'), 0)
d.save(OUT, garbage=4, deflate=True)
print('wrote', OUT, round(W/cm, 1), 'x', round(H/cm, 1), 'cm')
