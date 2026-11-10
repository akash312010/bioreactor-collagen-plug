"""Convert the Fluent ASCII exports of the R1T axial plane (x = 0, axial_plane/axial_plane_t*.txt) into
triangulated VTK (.vtu) files for ParaView, written to vtu/. Nodes are triangulated in the (z, y) plane;
triangles whose centroid lies outside the device section are removed. Coordinates are written in mm.
Run: python porous_plug_3d/tracer_R1T/plane_to_vtu.py"""
import sys, glob, os, re, numpy as np
from scipy.spatial import Delaunay
HERE = os.path.dirname(os.path.abspath(__file__))
SRC, OUT = os.path.join(HERE, 'axial_plane'), os.path.join(HERE, 'vtu')
os.makedirs(OUT, exist_ok=True)
ri, R, ro = 1.93e-3, 16.664e-3, 2.4995e-3
y_in, y_roof, y_floor, y_out = 0.0, -31.68e-3, -58.553e-3, -60.33e-3
def inside(z, y, tol=1e-7):
    a = np.abs(z)
    return (((y <= y_in+tol) & (y >= y_roof-tol) & (a <= ri+tol)) | ((y <= y_roof+tol) & (y >= y_floor-tol) & (a <= R+tol))
            | ((y <= y_floor+tol) & (y >= y_out-tol) & (a <= ro+tol)))
tri = None
for f in sorted(glob.glob(os.path.join(SRC, 'axial_plane_t*.txt'))):
    d = np.loadtxt(f, skiprows=1); z, y, c = d[:, 3], d[:, 2], d[:, 4]
    if tri is None:
        T = Delaunay(np.c_[z, y]).simplices
        cz, cy = z[T].mean(1), y[T].mean(1)
        T = T[inside(cz, cy)]; P = np.c_[d[:, 1]*0, y, z]; tri = T
    t = float(re.search(r'_t([0-9.]+)s\.txt$', f).group(1))
    name = os.path.join(OUT, f'R1T_plane_t{t:05.1f}s.vtu')
    with open(name, 'w') as o:
        o.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n<UnstructuredGrid>\n')
        o.write(f'<Piece NumberOfPoints="{len(P)}" NumberOfCells="{len(tri)}">\n<PointData Scalars="tracer">\n')
        o.write('<DataArray type="Float32" Name="tracer" format="ascii">\n' + ' '.join(f'{v:.5g}' for v in np.clip(c, 0, 1)) + '\n</DataArray>\n</PointData>\n')
        o.write('<Points>\n<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n'
                + '\n'.join(f'{a:.6e} {b:.6e} {cc:.6e}' for a, b, cc in P*1e3) + '\n</DataArray>\n</Points>\n')   # mm
        o.write('<Cells>\n<DataArray type="Int32" Name="connectivity" format="ascii">\n' + ' '.join(map(str, tri.ravel())) + '\n</DataArray>\n')
        o.write('<DataArray type="Int32" Name="offsets" format="ascii">\n' + ' '.join(map(str, 3*np.arange(1, len(tri)+1))) + '\n</DataArray>\n')
        o.write('<DataArray type="UInt8" Name="types" format="ascii">\n' + ' '.join(['5']*len(tri)) + '\n</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n')
    print(name, len(tri), 'triangles')
