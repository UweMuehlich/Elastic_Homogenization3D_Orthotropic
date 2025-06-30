import numpy as np

#----------------------------------------------------------------
def hexagon_points(r): # for regular hexagons only
    return [(r * math.cos(math.pi/3 * i), r * math.sin(math.pi/3 * i), 0) for i in range(6)]

#----------------------------------------------------------------
def create_parametric_hexagon(T, D, angle_deg):
    """
    Create a hexagon from top edge length T, side edge length D,
    and slant angle in degrees.

    Parameters:
        T : float
            Length of the horizontal top/bottom edges.
        D : float
            Length of the side (non-horizontal) edges.
        angle_deg : float
            Angle between the slanted edge and the vertical direction (in degrees).

    Returns:
        hex_points : np.ndarray of shape (6, 2)
            Vertices of the hexagon in counter-clockwise order.
    """
    angle_rad = np.radians(angle_deg)
    # Coordinates of the two Points in the first quadrant
    dx =  D * np.sin(angle_rad)
    dy =  D * np.cos(angle_rad) 
    
    points =  [( T/2, dy,0), (-T/2,  dy,0),(-T/2-dx, 0,0),
               (-T/2,-dy,0), ( T/2, -dy,0),( T/2+dx, 0,0)]
  
    return points 

#----------------------------------------------------------------
def offset_hexagon(hex_points, t):
    """
    Compute an inner hexagon offset inward by thickness t.

    Parameters:
        hex_points : array-like, shape (6, 2)
            Vertices of the outer hexagon.
        t : float
            Offset thickness.

    Returns:
        inner_hex : np.ndarray of shape (6, 2)
            Vertices of the inner hexagon.
    """

    def line_intersection(p, u, q, v): # helper function
        A1, A2,B1, B2,C1,C2 =  u[0], u[1], v[0], v[1], q[0] - p[0], q[1]-p[1]
        D = A1 * B2 - A2 * B1

        if np.abs(D) < 1e-10:
            raise ValueError("Lines are parallel.")
        lambda1 = (C1 * B2 - C2 * B1) / D
        lambda2 = (A1 * C2 - A2 * C1) / D

        return  np.array([p[0]+ lambda1*u[0], p[1]+lambda1*u[1]])

    hex_points = np.array(hex_points)

    H     = hex_points[0][1]                       #  |y
    W     = hex_points[5][0] - hex_points[0][0]    #  -------------o v0
    alpha = np.arctan(H/W)                         #  p  -> u     c \
    b     = t / np.sin(alpha)                      #                 \
                                                   #               q  o v5 __ x
    q = np.array([hex_points[5][0] - b, 0,0])
    p = np.array([0.,hex_points[0][1]-t,0])
    u = np.array([1,0,0])
    v = np.array(hex_points[0] - hex_points[5])
    c = line_intersection(p, u, q, v) 

    pp,cc,qq = (p[0], p[1], 0), (c[0], c[1], 0), (q[0], q[1], 0)
    ccm      = (-c[0], -c[1],0)
    points   = [cc, (-c[0], c[1] , 0), (-q[0], 0, 0), ccm, (c[0], -c[1], 0),qq]
  
    return points 

