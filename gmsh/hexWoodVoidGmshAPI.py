import gmsh
import sys
import math
import meshio
import argparse
import ast

from pathlib import Path

from helpers import parametric_hexagon as ph

def create_hex_polygon(points):
    lines = []
    points_tags = []
    # Create points once and store tags for re-use if needed
    for pt in points:
        points_tags.append(gmsh.model.occ.addPoint(*pt))
    for i in range(len(points)):
        line = gmsh.model.occ.addLine(points_tags[i], points_tags[(i+1) % len(points)])
        lines.append(line)
    wire = gmsh.model.occ.addWire(lines)
    return wire, lines
#----------------------------------------------------------------
def hexagon_edge_midpoints(points):
    """
    Given a list of 6 points defining a hexagon in order,
    return a list of 6 midpoints of its edges.

    Parameters:
        points (list of tuples): List of 6 (x, y) or (x, y, z) tuples.

    Returns:
        list of tuples: List of 6 midpoints.
    """
    if len(points) != 6:
        raise ValueError("Expected 6 points for a hexagon.")
    
    midpoints = []
    for i in range(6):
        p1 = points[i]
        p2 = points[(i + 1) % 6]  # Wrap around to first point
        midpoint = tuple((a + b) / 2 for a, b in zip(p1, p2))
        midpoints.append(midpoint)

    return midpoints
#----------------------------------------------------------------
def isclose3(p1, p2, tol=1e-6):
    return all(math.isclose(a, b, abs_tol=tol) for a, b in zip(p1, p2))
#----------------------------------------------------------------
def find_surface_tag_by_midpoint(target_midpoint, surfaces, tol=1e-6):
    """
   Find the surface tag whose midpoint matches the target coordinates.

    Parameters:
        target_midpoint: tuple (x, y, z)
       surfaces: list of (2, tag) pairs
       tol: tolerance for coordinate match

    Returns:
        tag of the matching surface

    Raises:
        ValueError if no match is found
    """
    for dim, tag in surfaces:
        midpoint = gmsh.model.occ.getCenterOfMass(dim, tag)
        if isclose3(midpoint, target_midpoint, tol):
            return tag
    raise ValueError(f"No surface found matching midpoint {target_midpoint}")
#----------------------------------------------------------------    
def get_front_faces(target_z, surfaces, tol=1e-6):
    """Return all surfaces lying on z = target_z with their areas."""
    faces = []
    for dim, tag in surfaces:
        if dim != 2:
            continue
        x, y, z = gmsh.model.occ.getCenterOfMass(dim, tag)
        if abs(z - target_z) < tol:
            area = gmsh.model.occ.getMass(dim, tag)
            faces.append((tag, area))
    return sorted(faces, key=lambda pair: pair[1])  # sort by area

#----------------------------------------------------------------
#                M A I N
#----------------------------------------------------------------
def main(): 
    parser = argparse.ArgumentParser(description="Hexagon GMSH for FEniCSx")
    
    parser.add_argument("--ID", default="three_hex_rings", help="Problem ID")
    parser.add_argument("--T" , type=float, default=2.0, help="Width")
    parser.add_argument("--D" , type=float,default=2.5, help="Length slanted edge")
    parser.add_argument("--angle"   , type=float, default=60,  help="angle in degree")
    parser.add_argument("--thick"   , type=float, default=0.4, help="Thickness")
    parser.add_argument("--thickOut", type=float, default=0.2, help="Thickness outer ring")
    parser.add_argument("--height"  , type=float, default=1.0, help="Extrusion length")
    parser.add_argument("--MeshSize", type=float, default=0.5, help="Global mesh size")
    parser.add_argument("--IntOrder", type=int  , default=1,   help="Finite Element order (1,2)")
    parser.add_argument("--inner", action="store_true", help="With inner core (default: not included)")
    parser.add_argument("--prisms", action="store_true", help="use prisms (default: tetrahedrons)")

    args = parser.parse_args()
    
    GMSH_ID, T, D, angle_deg    = args.ID, args.T, args.D, args.angle
    thickness, global_mesh_size = args.thick, args.MeshSize
    int_order, height           = args.IntOrder, args.height
    thick_out, inner_core       = args.thickOut, args.inner
    prisms                      = args.prisms

    output_dir = Path("mesh")
    output_dir.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.model.add(GMSH_ID)
    
    hex_outer_points   = ph.create_parametric_hexagon(T, D, angle_deg)
    hex_middle_points  = ph.offset_hexagon(hex_outer_points, thick_out)
    hex_inner_points   = ph.offset_hexagon(hex_outer_points, thickness)

    # Create hex polygons and get their edges (lines)
    hex_outer, lines_outer   = create_hex_polygon(hex_outer_points)
    hex_middle, lines_middle = create_hex_polygon(hex_middle_points)
    hex_inner, lines_inner   = create_hex_polygon(hex_inner_points)

    surf_outer  = gmsh.model.occ.addPlaneSurface([hex_outer])
    surf_middle = gmsh.model.occ.addPlaneSurface([hex_middle])
    surf_inner  = gmsh.model.occ.addPlaneSurface([hex_inner])

    ring_outer = gmsh.model.occ.cut([(2, surf_outer)], [(2, surf_middle)], 
                                      removeObject=True, removeTool=False)[0][0][1]
    ring_middle = gmsh.model.occ.cut([(2, surf_middle)], [(2, surf_inner)],
                                       removeObject=True, removeTool=False)[0][0][1]
    core = surf_inner
    
    gmsh.model.occ.synchronize()

    # Extrude volumes
    vol_outer  = gmsh.model.occ.extrude([(2, ring_outer)], 0, 0,  height)
    vol_middle = gmsh.model.occ.extrude([(2, ring_middle)], 0, 0, height)
    vol_core   = gmsh.model.occ.extrude([(2, core)], 0, 0, height)

    # avoid double nodes (ensure perfect merging)
    core_tag   = [tag for (dim, tag) in vol_core   if dim == 3][0]
    middle_tag = [tag for (dim, tag) in vol_middle if dim == 3][0]
    outer_tag  = [tag for (dim, tag) in vol_outer  if dim == 3][0]

    gmsh.model.occ.fragment([(3, outer_tag), (3, middle_tag), (3, core_tag)], [])
    #gmsh.model.occ.fragment([(3, outer_tag), (3, middle_tag)], [])
    
    if not inner_core: # remove inner core if not requested
           gmsh.model.occ.remove([(3, core_tag)], recursive=True)
    
    gmsh.model.occ.synchronize()

    # --------------------
    # Physical Groups for volumes
    # --------------------
    gmsh.model.addPhysicalGroup(3, [vol_outer[1][1]], tag=1)
    gmsh.model.setPhysicalName(3, 1, "OuterRing")

    gmsh.model.addPhysicalGroup(3, [vol_middle[1][1]], tag=2)
    gmsh.model.setPhysicalName(3, 2, "MiddleRing")

    gmsh.model.addPhysicalGroup(3, [vol_core[1][1]], tag=3)
    gmsh.model.setPhysicalName(3, 3, "InnerCore")

    # --------------------
    # Tag individual edges of innermost hexagon (base)
    # --------------------
    for i, line_tag in enumerate(lines_inner):
        pg_tag = 10 + i  # unique tags starting at 10
        gmsh.model.addPhysicalGroup(1, [line_tag], tag=pg_tag)
        gmsh.model.setPhysicalName(1, pg_tag, f"InnerHexEdge_{i+1}")

    # --------------------
    # Tag individual surfaces of the outer boundary
    # --------------------    

    all_surfaces = gmsh.model.getEntities(2)
    
    outer_mid_points = hexagon_edge_midpoints(hex_outer_points)
    outer_mpts       = [(x, y, height/2.) for (x, y, _) in outer_mid_points]
    
    tag_on  = find_surface_tag_by_midpoint(outer_mpts[0], all_surfaces)
    tag_onw = find_surface_tag_by_midpoint(outer_mpts[1], all_surfaces)
    tag_osw = find_surface_tag_by_midpoint(outer_mpts[2], all_surfaces)
    tag_os  = find_surface_tag_by_midpoint(outer_mpts[3], all_surfaces)
    tag_ose = find_surface_tag_by_midpoint(outer_mpts[4], all_surfaces)
    tag_one = find_surface_tag_by_midpoint(outer_mpts[5], all_surfaces)
        
    # Assign physical groups and names
    pg_on  = gmsh.model.addPhysicalGroup(2, [tag_on])
    pg_onw = gmsh.model.addPhysicalGroup(2, [tag_onw])
    pg_osw = gmsh.model.addPhysicalGroup(2, [tag_osw])
    pg_os  = gmsh.model.addPhysicalGroup(2, [tag_os])
    pg_ose = gmsh.model.addPhysicalGroup(2, [tag_ose])
    pg_one = gmsh.model.addPhysicalGroup(2, [tag_one])

    gmsh.model.setPhysicalName(2, pg_on,  "outer_north_face")
    gmsh.model.setPhysicalName(2, pg_onw, "outer_north_west_face")
    gmsh.model.setPhysicalName(2, pg_osw, "outer_south_west_face")
    gmsh.model.setPhysicalName(2, pg_os,  "outer_south_face")
    gmsh.model.setPhysicalName(2, pg_ose, "outer_south_east_face")
    gmsh.model.setPhysicalName(2, pg_one, "outer_north_east_face")

    front_faces = get_front_faces(0.0, all_surfaces, tol=1e-6)
    back_faces  = get_front_faces(height, all_surfaces, tol=1e-6)

    print("front_faces =", front_faces)
    print("First element:", front_faces[0])

    gmsh.model.addPhysicalGroup(dim=2, tags=[int(tag[0]) for tag in front_faces], tag=1)
    gmsh.model.setPhysicalName(dim=2, tag=1, name="FrontFace")

    gmsh.model.addPhysicalGroup(dim=2, tags=[int(tag[0]) for tag in back_faces], tag=2)
    gmsh.model.setPhysicalName(dim=2, tag=2, name="BackFace")

    for i, (tag, _) in enumerate(front_faces, start=1):
        pg = gmsh.model.addPhysicalGroup(2, [tag])
        gmsh.model.setPhysicalName(2, pg, f"front{i}")
    
    for i, (tag, _) in enumerate(back_faces, start=1):
        pg = gmsh.model.addPhysicalGroup(2, [tag])
        gmsh.model.setPhysicalName(2, pg, f"back{i}")
    
    gmsh.option.setNumber("Mesh.ElementOrder", int_order)

    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", global_mesh_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 2*global_mesh_size)

    if inner_core:
      print("setting background field")
      Field_Radius = (D*np.cos(np.deg2rad(angle_deg))-thickness)*0.75
      cyl_field = gmsh.model.mesh.field.add("Cylinder")
      gmsh.model.mesh.field.setNumber(cyl_field, "Radius", Field_Radius)
      gmsh.model.mesh.field.setNumber(cyl_field, "VIn", Field_Radius)  # coarse mesh in center
      gmsh.model.mesh.field.setNumber(cyl_field, "VOut",global_mesh_size) # finer mesh outside
      gmsh.model.mesh.field.setNumber(cyl_field, "XCenter",  0.0)
      gmsh.model.mesh.field.setNumber(cyl_field, "YCenter",  0.0)
      gmsh.model.mesh.field.setNumber(cyl_field, "ZCenter", -0.1)
      gmsh.model.mesh.field.setNumber(cyl_field, "ZAxis", 2*height)

      # --- Set background field ---
      gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
      gmsh.option.setNumber("Mesh.Algorithm", 6)  # e.g., 6 = DelQuad, 5 = Delaunay

      gmsh.model.mesh.field.setAsBackgroundMesh(cyl_field)
       
    gmsh.model.mesh.generate(3)

    # Save mesh in msh 2.2 format
    #gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    mesh_file = (output_dir / GMSH_ID).with_suffix(".msh")
    gmsh.write(str(mesh_file))


    gmsh.finalize()

    # -------------------
    # Convert to XDMF
    # -------------------
    msh = meshio.read(mesh_file)

    # Write full mesh with physical groups preserved
    xdmf_file = (output_dir / GMSH_ID).with_suffix(".xdmf")
    meshio.write(xdmf_file,msh)

    # -------------------
    # Write Data
    # -------------------
    filename = f"{output_dir}/{GMSH_ID}-corners.dat"

    # The values to write
    
    hex_outer_points_back = [(x, y, height) for (x, y, _) in hex_outer_points] 
    hex_points_all        = hex_outer_points + hex_outer_points_back
    with open(filename, 'w') as f:
        for point in hex_points_all:
            x, y, z = point
            f.write(f"({float(x)}, {float(y)}, {float(z)})\n")
		
if __name__ == "__main__":
    main()
