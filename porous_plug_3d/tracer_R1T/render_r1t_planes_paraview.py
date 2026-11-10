"""ParaView (pvbatch/pvpython) rendering of the R1T tracer on the axial plane of the 3D device.
Input: R1T_plane_t*.vtu written by plane_to_vtu.py. Output: one vector PDF per time plus a colour-bar PDF.
Turbo colour map, range 0-1, white background, panel letters added when composing.
Run (from the repository root): pvbatch porous_plug_3d/tracer_R1T/render_r1t_planes_paraview.py
Reads porous_plug_3d/tracer_R1T/vtu/, writes figures/r1t_panels/."""
import sys, glob, os
from paraview.simple import *
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'vtu'); OUT = os.path.join(HERE, '..', '..', 'figures', 'r1t_panels')
os.makedirs(OUT, exist_ok=True)
files = sorted(glob.glob(os.path.join(SRC, 'R1T_plane_t*.vtu')))
lut = GetColorTransferFunction('tracer'); lut.AutomaticRescaleRangeMode = 'Never'
lut.ApplyPreset('Turbo', True); lut.RescaleTransferFunction(0.0, 1.0)
GetOpacityTransferFunction('tracer').RescaleTransferFunction(0.0, 1.0)
ycen, half = -43.0, 19.0                     # mm: chamber plus the lowest 7 mm of the inlet conduit
for f in files:
    v = CreateRenderView(); v.ViewSize = [700, 760]; v.OrientationAxesVisibility = 0
    v.Background = [1, 1, 1]
    try: v.UseColorPaletteForBackground = 0
    except AttributeError: pass
    r = XMLUnstructuredGridReader(FileName=[f]); d = Show(r, v)
    ColorBy(d, ('POINTS', 'tracer')); d.LookupTable = lut; d.SetScalarBarVisibility(v, False)
    lut.RescaleTransferFunction(0.0, 1.0)
    o = Show(FeatureEdges(Input=ExtractSurface(Input=r)), v); o.DiffuseColor = [0, 0, 0]; o.AmbientColor = [0, 0, 0]; o.LineWidth = 1.5
    o.ColorArrayName = ['POINTS', '']
    v.CameraPosition = [100, ycen, 0]; v.CameraFocalPoint = [0, ycen, 0]; v.CameraViewUp = [0, 1, 0]
    v.CameraParallelProjection = 1; v.CameraParallelScale = half
    stem = os.path.join(OUT, os.path.basename(f)[:-4])
    ExportView(stem + '.pdf', view=v, Rasterize3Dgeometry=0)
    SaveScreenshot(stem + '.png', v, ImageResolution=[700, 760])
    Delete(v)
# colour bar alone
v = CreateRenderView(); v.ViewSize = [500, 110]; v.OrientationAxesVisibility = 0; v.Background = [1, 1, 1]
try: v.UseColorPaletteForBackground = 0
except AttributeError: pass
r = XMLUnstructuredGridReader(FileName=[files[0]]); d = Show(r, v); ColorBy(d, ('POINTS', 'tracer')); d.LookupTable = lut
d.Opacity = 0.0
sb = GetScalarBar(lut, v); sb.Visibility = 1; sb.Orientation = 'Horizontal'; sb.WindowLocation = 'Any Location'
sb.Position = [0.05, 0.12]; sb.ScalarBarLength = 0.9; sb.Title = 'Tracer concentration'; sb.ComponentTitle = ''
sb.TitleColor = [0, 0, 0]; sb.LabelColor = [0, 0, 0]; sb.TitleFontSize = 16; sb.LabelFontSize = 15
sb.TitleFontFamily = 'Arial'; sb.LabelFontFamily = 'Arial'; sb.AutomaticLabelFormat = 0; sb.LabelFormat = '%.1f'
sb.RangeLabelFormat = '%.0f'; sb.UseCustomLabels = 1; sb.CustomLabels = [0.2, 0.4, 0.6, 0.8]; sb.AddRangeLabels = 1
sb.ScalarBarThickness = 16
ExportView(os.path.join(OUT, 'R1T_colorbar.pdf'), view=v, Rasterize3Dgeometry=0)
SaveScreenshot(os.path.join(OUT, 'R1T_colorbar.png'), v, ImageResolution=[500, 110])
print('wrote', len(files), 'panels to', OUT)
