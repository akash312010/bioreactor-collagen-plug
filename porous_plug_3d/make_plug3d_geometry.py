"""Three-body geometry of the bioreactor for the 3D porous-plug check (units: mm).
Bodies share their interface faces (conformal): inlet conduit, chamber (porous plug), outlet conduit.
Axis = y (vertical), inlet at y = 0, flow toward -y. Dimensions are the as-meshed values of the
open-device grids (same as Fig. 1 and the axisymmetric Darcy model)."""
import gmsh, math
ri, R, ro = 1.93, 16.664, 2.4995           # inlet, chamber, outlet radii (mm)
y_in, y_roof, y_floor, y_out = 0.0, -31.68, -58.553, -60.33
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("plug3d")
occ = gmsh.model.occ
# cylinders along -y: addCylinder(x, y, z, dx, dy, dz, r)
cin  = occ.addCylinder(0, y_in,    0, 0, y_roof - y_in,    0, ri)
cch  = occ.addCylinder(0, y_roof,  0, 0, y_floor - y_roof, 0, R)
cout = occ.addCylinder(0, y_floor, 0, 0, y_out - y_floor,  0, ro)
out, _ = occ.fragment([(3, cch)], [(3, cin), (3, cout)])
occ.synchronize()
vols = [t for d, t in gmsh.model.getEntities(3)]
info = {}
for v in vols:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, v)
    info[v] = (ymin, ymax, occ.getMass(3, v))
for v, (a, b, m) in sorted(info.items(), key=lambda kv: -kv[1][1]):
    print(f"volume {v}: y {a:.3f} to {b:.3f} mm, volume {m:.1f} mm^3")
chamber = max(vols, key=lambda v: info[v][2])
print("chamber volume (mL):", info[chamber][2]/1000, " pore volume at phi=0.27 (mL):", 0.27*info[chamber][2]/1000)
gmsh.write("bioreactor_plug3d_3bodies.step")
# quick test mesh to confirm the geometry meshes conformally
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 1.0)
gmsh.model.mesh.generate(3)
print("test tetra mesh nodes:", len(gmsh.model.mesh.getNodes()[0]))
gmsh.finalize()
