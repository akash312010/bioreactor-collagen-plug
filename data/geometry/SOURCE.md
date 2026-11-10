# Source of the pore geometry

`collagen_SHG_Niestrawska2022_Fig4.png` is a second-harmonic-generation image of healthy
aortic adventitial collagen, reproduced from

J. A. Niestrawska, A. Pukaluk, A. R. Babu and G. A. Holzapfel, "Differences in collagen
fiber diameter and waviness between healthy and aneurysmal abdominal aortas,"
Microsc. Microanal. 28, 1649-1663 (2022), Fig. 4, under CC BY 4.0.
Annotations are part of the source image.

The pore channels G1-G5 were traced from the dark interstitial space of this image at five
intensity thresholds (28th to 45th percentile):

| Channel | File | Threshold |
|---|---|---|
| G1 | `collagen_channel_interior_p28.csv` | 28th percentile |
| G2 | `collagen_channel_interior_p31.csv` | 31st percentile |
| G3 (reference) | `collagen_channel_interior_p35.csv` | 35th percentile |
| G4 | `collagen_channel_interior_p40.csv` | 40th percentile |
| G5 | `collagen_channel_interior_p45.csv` | 45th percentile |

The CSV files give the wall and centreline coordinates in millimetres. `cad_info.json` lists, for each
channel, the threshold, length, area, widths and tortuosity of the traced geometry (micrometres).
