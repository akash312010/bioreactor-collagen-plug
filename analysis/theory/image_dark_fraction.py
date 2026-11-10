"""Dark-area fraction of the source SHG image at the five channel thresholds (porosity sensitivity in the paper).

The threshold intensities of G1-G5 are those used to trace the channels (data/geometry/cad_info.json,
'threshold_value', grey level / 255). For each, the fraction of all image pixels at or below that grey level
is an areal estimate of the dark (non-collagen) space; for an isotropic section it estimates the volume
fraction (Delesse principle). Because dark regions also contain non-collagenous tissue, it is an upper
estimate of the pore space. With k proportional to phi at fixed geometry, the plug pressure demand is
divided by fraction/0.27.
Run:  python analysis/theory/image_dark_fraction.py      -> image_dark_fraction.json
"""
import json
from pathlib import Path
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
img = np.asarray(Image.open(ROOT/'data'/'geometry'/'collagen_SHG_Niestrawska2022_Fig4.png').convert('L')).astype(int)
cad = json.loads((ROOT/'data'/'geometry'/'cad_info.json').read_text())
CH = {28: 'G1', 31: 'G2', 35: 'G3', 40: 'G4', 45: 'G5'}
PHI0, DP_G3_5 = 0.27, 31.56          # assumed porosity; port-fed G3 demand at 5 mL/min from the 3D run R1 (kPa)
out = {}
for e in cad:
    if e['variant'] != 'interior':
        continue
    level = int(round(e['threshold_value']*255))
    frac = float((img <= level).mean())
    out[CH[e['threshold_pct']]] = dict(threshold_pct=e['threshold_pct'], grey_level=level, dark_fraction=frac,
                                       permeability_factor=frac/PHI0, demand_divided_by=frac/PHI0)
out['G3']['portfed_demand_5mLmin_kPa'] = DP_G3_5*PHI0/out['G3']['dark_fraction']
(HERE/'image_dark_fraction.json').write_text(json.dumps(out, indent=1))
for g, v in out.items():
    print(f"{g}: level {v['grey_level']:3d}  fraction {v['dark_fraction']:.3f}  k factor {v['permeability_factor']:.2f}")
print(f"G3 port-fed demand at 5 mL/min with the image fraction: {out['G3']['portfed_demand_5mLmin_kPa']:.1f} kPa")
