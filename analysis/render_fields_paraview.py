"""Run with official ParaView 6.2 pvpython --no-mpi --dr --force-offscreen-rendering.
Native polygon cell scalars are rendered flat; only psi contours use point data.
Vector export disables 3D rasterization. No CFD solution is generated here.
"""
from paraview.simple import *
from paraview.simple import _DisableFirstRenderCameraReset
from paraview import servermanager
from pathlib import Path
import json,argparse,xml.etree.ElementTree as ET
import vtk

def context_svg(calc,lut,reader,contours,window,size,path,placed_width):
 """VTK Context2D vector path: modern GL2PS rasterizes 3D props despite its flag.
 Use ParaView's own fetched cells, transfer function, and contour filter in
 vtkPolyDataItem, then vtkSVGExporter. No matplotlib field plotting is used.
 """
 from vtkmodules.util.numpy_support import vtk_to_numpy
 data=servermanager.Fetch(calc)
 surf=vtk.vtkDataSetSurfaceFilter();surf.SetInputData(data);surf.Update()
 tri=vtk.vtkTriangleFilter();tri.SetInputConnection(surf.GetOutputPort());tri.Update()
 poly=tri.GetOutput()
 x0,x1,y0,y1=window
 if contours=='negative':
  for origin,normal in [((x0,0,0),(1,0,0)),((x1,0,0),(-1,0,0)),((0,y0,0),(0,1,0)),((0,y1,0),(0,-1,0))]:
   plane=vtk.vtkPlane();plane.SetOrigin(*origin);plane.SetNormal(*normal)
   clip=vtk.vtkClipPolyData();clip.SetInputData(poly);clip.SetClipFunction(plane);clip.Update()
   new=vtk.vtkPolyData();new.DeepCopy(clip.GetOutput());poly=new
 transform=vtk.vtkTransform();transform.Scale(size[0]/(x1-x0),size[1]/(y1-y0),1);transform.Translate(-x0,-y0,0)
 tf=vtk.vtkTransformPolyDataFilter();tf.SetInputData(poly);tf.SetTransform(transform);tf.Update()
 poly=tf.GetOutput()
 colors=lut.GetClientSideObject().MapScalars(poly.GetCellData().GetArray('display_scalar'),0,0)
 item=vtk.vtkPolyDataItem();item.SetPolyData(poly);item.SetScalarMode(2);item.SetMappedColors(colors)
 ctx=vtk.vtkContextView();ctx.GetRenderWindow().SetOffScreenRendering(1);ctx.GetRenderWindow().SetSize(*size)
 ctx.GetRenderer().SetBackground(1,1,1);ctx.GetScene().AddItem(item)
 if contours:
  cf=vtk.vtkContourFilter();cf.SetInputData(servermanager.Fetch(reader));cf.SetInputArrayToProcess(0,0,0,0,'psi')
  levels=[.1*i for i in range(1,10)] if contours=='positive' else [-.00002,-.0001,-.0003,-.0006,-.001]
  for i,v in enumerate(levels):cf.SetValue(i,v)
  cf.Update();ct=vtk.vtkTransformPolyDataFilter();ct.SetInputConnection(cf.GetOutputPort());ct.SetTransform(transform);ct.Update()
  line=vtk.vtkPolyDataItem();line.SetPolyData(ct.GetOutput());line.SetScalarMode(2)
  black=vtk.vtkUnsignedCharArray();black.SetNumberOfComponents(4);black.SetNumberOfTuples(ct.GetOutput().GetNumberOfCells())
  line_color=255 if contours=='negative' else 0
  for i in range(black.GetNumberOfTuples()):black.SetTuple4(i,line_color,line_color,line_color,255)
  line.SetMappedColors(black);ctx.GetScene().AddItem(line)
 ctx.GetRenderWindow().Render()
 exp=vtk.vtkSVGExporter();exp.SetRenderWindow(ctx.GetRenderWindow());exp.SetFileName(str(path));exp.Write()
 ctx.GetRenderWindow().Finalize()
 # Graphics-only merge: native flat cell colours are unchanged. Nonzero-fill
 # compound paths remove edges inside each colour band. A 0.12 pt matching
 # stroke traps antialias gaps between bands; no field interpolation is used.
 namespace='http://www.w3.org/2000/svg';ET.register_namespace('',namespace)
 root=ET.parse(path).getroot();group=root.find('{'+namespace+'}g');buckets={};ntri=0
 for child in list(group):
  if child.tag=='{'+namespace+'}path' and child.get('fill-opacity')!='0':
   buckets.setdefault(child.get('fill'),[]).append(child.get('d'));group.remove(child);ntri+=1
 scale=placed_width/size[0]
 for index,(color,parts) in enumerate(buckets.items()):
  el=ET.Element('{'+namespace+'}path',{'fill':color,'fill-rule':'nonzero','stroke':color,
   'stroke-width':str(.12/scale),'stroke-linejoin':'round','d':' '.join(parts)})
  group.insert(index+1,el)
 for line_el in group.iter('{'+namespace+'}line'):
  line_el.set('stroke-width',str(.45/scale));line_el.set('stroke-linecap','round')
 # Native external boundary only; the mesh remains visible solely in Fig.1.
 # The device cut stores separate polygon vertices. Merge only coincident
 # geometry for boundary extraction, otherwise every cell edge looks external.
 clean=vtk.vtkCleanPolyData();clean.SetInputData(poly);clean.ToleranceIsAbsoluteOn()
 clean.SetAbsoluteTolerance(1e-8);clean.Update()
 edges=vtk.vtkFeatureEdges();edges.SetInputConnection(clean.GetOutputPort());edges.BoundaryEdgesOn()
 edges.FeatureEdgesOff();edges.ManifoldEdgesOff();edges.NonManifoldEdgesOff();edges.Update()
 edge=edges.GetOutput();outline=[]
 for i in range(edge.GetNumberOfCells()):
  cell=edge.GetCell(i)
  pts=[edge.GetPoint(cell.GetPointId(j)) for j in range(cell.GetNumberOfPoints())]
  outline.append('M'+' L'.join(f'{v[0]:.6f},{size[1]-v[1]:.6f}' for v in pts))
 ET.SubElement(group,'{'+namespace+'}path',{'fill':'none','stroke':'#000000',
  'stroke-width':str(.45/scale),'stroke-linejoin':'round','d':' '.join(outline)})
 ET.ElementTree(root).write(path,encoding='unicode')
 return {'native_triangle_paths':ntri,'compound_colour_paths':len(buckets),'trapping_pt':.12,'contour_and_outline_pt':.45}

ROOT=Path(__file__).resolve().parents[1]
_DisableFirstRenderCameraReset()
OUT=ROOT/'figures/panels';OUT.mkdir(parents=True,exist_ok=True)
DATA=ROOT/'data/fields'
specs=[
 ('G3_speed','G3','speed/0.004164757602156477',[0,4],[-1,91.298103,-13.4423575,12.466524],'positive'),
 ('G3_pressure','G3','p/260.8280671',[0,1],[-1,91.298103,-13.4423575,12.466524],None),
 ('G1_speed','G1','speed/0.0031154692934415207',[0,4],[-1,91.3,-13.4423575,12.466524],'positive'),
 ('G1_pocket','G1','speed/0.0031154692934415207',[0,.08],[37,44.5,-11.5,-3.5],'negative'),
 ('Device_speed','Device_cut','axial_speed',[0,.045],[-19,19,-62,2],None),
]
for name in ['E11_t1','E11_t3','E12_t1','E12_t3','E8_t2','E8_t10']:
 specs.append((name,name,'log10(max(1-uds0,0.001))',[-3,0],[-1,91.3,-13.4423575,12.466524],None))
parser=argparse.ArgumentParser();parser.add_argument('--only');args=parser.parse_args()
manifest=json.loads((OUT/'render_manifest.json').read_text()) if args.only and (OUT/'render_manifest.json').exists() else []
for stem,filename,formula,limits,window,contours in specs:
 if args.only and stem!=args.only:continue
 manifest=[v for v in manifest if v['panel']!=stem]
 print('RENDER',stem,flush=True)
 ResetSession()
 reader=LegacyVTKReader(FileNames=[str(DATA/(filename+'.vtk'))]);reader.UpdatePipeline()
 calc=Calculator(Input=reader);calc.AttributeType='Cell Data';calc.ResultArrayName='display_scalar';calc.Function=formula;calc.UpdatePipeline()
 view=CreateView('RenderView');view.ViewSize=[1200,round(1200*(window[3]-window[2])/(window[1]-window[0]))]
 view.Background=[1,1,1];view.UseColorPaletteForBackground=0;view.OrientationAxesVisibility=0
 view.CameraParallelProjection=1;view.InteractionMode='2D'
 display=Show(calc,view);display.Representation='Surface';display.Ambient=1;display.Diffuse=0
 ColorBy(display,('CELLS','display_scalar'));display.InterpolateScalarsBeforeMapping=0
 lut=GetColorTransferFunction('display_scalar');lut.ColorSpace='RGB'
 assert lut.ApplyPreset('Turbo',True),'Official Turbo preset unavailable'
 lut.Discretize=1;lut.NumberOfTableValues=256
 lut.RescaleTransferFunction(*limits);lut.NanColor=[1,1,1]
 display.SetScalarBarVisibility(view,False)
 if contours:
  con=Contour(Input=reader);con.ContourBy=['POINTS','psi']
  con.Isosurfaces=[.1*i for i in range(1,10)] if contours=='positive' else [-.00002,-.0001,-.0003,-.0006,-.001]
  cd=Show(con,view);ColorBy(cd,None);cd.DiffuseColor=[0,0,0];cd.AmbientColor=[0,0,0];cd.LineWidth=1
 x0,x1,y0,y1=window;cx=(x0+x1)/2;cy=(y0+y1)/2
 view.CameraPosition=[cx,cy,100];view.CameraFocalPoint=[cx,cy,0];view.CameraViewUp=[0,1,0]
 view.CameraParallelScale=(y1-y0)/2
 Render(view)
 ExportView(str(OUT/(stem+'.pdf')),view=view,Rasterize3Dgeometry=0)
 ExportView(str(OUT/(stem+'.svg')),view=view,Rasterize3Dgeometry=0)
 SaveScreenshot(str(OUT/(stem+'.png')),view,ImageResolution=view.ViewSize)
 lut.UpdateVTKObjects()
 values=[]
 for i in range(256):
  rgb=[0.,0.,0.];lut.GetClientSideObject().GetColor(limits[0]+i*(limits[1]-limits[0])/255,rgb);values.append(rgb)
 (OUT/'turbo_colors.json').write_text(json.dumps(values))
 placed_width=125 if stem=='Device_speed' else 113 if stem=='G1_pocket' else 225 if stem.startswith('E') else 337
 vector_stats=context_svg(calc,lut,reader,contours,window,list(view.ViewSize),OUT/(stem+'_vector.svg'),placed_width)
 manifest.append({'panel':stem,'source':'data/fields/'+filename+'.vtk','formula':formula,
  'color_limits':limits,'world_window':window,'view_size':list(view.ViewSize),
  'renderer':'ParaView 6.2.0 official Windows distribution','rasterize_3d_geometry':False,
  'cell_scalar_interpolation':False,'contours':contours,
  'vector_export':'Native ParaView cells and Turbo256 LUT; VTK Context2D export, compound colour paths with matching-colour trapping',
  'colour_preset':'Turbo','colour_bands':256,'vector_path_statistics':vector_stats,
  'gl2ps_export_note':'Modern OpenGL2 GL2PS export is raster despite flag; use _vector.svg and its converted PDF.'})
 (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2))
print('DONE',flush=True)
