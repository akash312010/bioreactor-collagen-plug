"""Main-text Fig. 7: (a,b) plug-scale speed and replacement curves (make_plug_transport_figure.py) on top,
(c-e) ParaView tracer snapshots of the 3D run R1T at 25, 75 and 150 s below, with one colour bar. Vector PDF.
Run: python porous_plug_3d/tracer_R1T/compose_figure7.py
(after analysis/porous/make_plug_transport_figure.py and render_r1t_planes_paraview.py)
Writes figures/Figure_Plug_Transport_V3.pdf."""
import sys
from pathlib import Path
import pymupdf as fitz

ROOT = Path(__file__).resolve().parents[2]
TOP, PV, OUT = ROOT/'figures'/'Figure_Plug_Transport_V2.pdf', ROOT/'figures'/'r1t_panels', str(ROOT/'figures'/'Figure_Plug_Transport_V3.pdf')
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
top = fitz.open(TOP); tw, th = top[0].rect.width, top[0].rect.height
W = tw                                   # 17 cm
pw = 4.6*cm; ph = pw*760/700             # snapshot panels keep the ParaView aspect ratio
gx = (W - 3*pw)/4                        # equal gaps, panels centred
y0 = th + 0.75*cm                        # room for the panel letters
cw = 0.8*W; ch = cw*110/500
H = y0 + ph + 0.15*cm + ch + 0.1*cm
d = fitz.open(); p = d.new_page(width=W, height=H)
p.show_pdf_page(fitz.Rect(0, 0, tw, th), top, 0)
for k, t in enumerate(['025.0s', '075.0s', '150.0s']):
    x = gx + k*(pw + gx)
    p.show_pdf_page(fitz.Rect(x, y0, x+pw, y0+ph), fitz.open(PV/f'R1T_plane_t{t}.pdf'), 0)
    p.insert_text((x - 0.45*cm, y0 - 0.2*cm), f'({"cde"[k]})', fontname='ArialRegularT', fontfile=FONT, fontsize=12)
    tlabel(p, x + 0.1*cm, y0 + 0.45*cm, int(float(t[:-1])))
cx = (W - cw)/2; cy = y0 + ph + 0.15*cm
p.show_pdf_page(fitz.Rect(cx, cy, cx+cw, cy+ch), fitz.open(PV/'R1T_colorbar.pdf'), 0)
d.save(OUT, garbage=4, deflate=True)
print('wrote', OUT, round(W/cm, 1), 'x', round(H/cm, 1), 'cm')
