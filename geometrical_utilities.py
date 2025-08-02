# -*- coding: utf-8 -*-
import numpy as np
import math

"""
Polygon class.Used for computational geometry tasks, such as checking whether
a point lies inside a polygon
"""
class Polygon:
    
    def __init__(self,vertexes):
        """
        Initialize the polygon with a list of vertex points.
        Each vertex is represented as a tuple (x, y).
        """
        self.vertexes = vertexes
        self.sides = len(vertexes)
    
    
    def edges(self): 
        """
        Returns a list of Segment objects representing the edges of the polygon.
        Each edge connects consecutive vertices, with the last vertex connecting back to the first.
        """
        edges=[]
        for i in range(self.sides):
            edges.append(Segment(self.vertexes[i],self.vertexes[(i+1)%self.sides]))
        return edges
    
    def is_point_inside(self, point):
        """
        Checks if a given point is inside the polygon using the ray-casting algorithm.
        Args:
            point (tuple): The point to check, represented as (x, y).
        Returns:
            bool: True if the point is inside or on the edge of the polygon, False otherwise.
        """
        count = 0 # Count of ray intersections with polygon edges
        x, y = point[0], point[1]
        
        for edge in self.edges():
            x1, y1 = edge.origin[0], edge.origin[1]
            x2, y2 = edge.end[0], edge.end[1]
            
            # Check if the point is exactly on an edge
            if min(y1, y2) <= y <= max(y1, y2) and min(x1, x2) <= x <= max(x1, x2):
                if (x2 - x1) * (y - y1) == (x - x1) * (y2 - y1):  # Cross product == 0
                    return True  # On the edge
    
            # Ensure the ray doesn't pass directly through a vertex (adjust y slightly if needed)
            if y1 > y2:
                x1, x2, y1, y2 = x2, x1, y2, y1
            if y == y1 or y == y2:  # Avoid edge case where ray hits vertex
                y += 0.00001
                
            # Check if the horizontal ray intersects the edge    
            if y1 < y < y2:  # Intersection is to the right of the point
                x_intersection = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
                if x_intersection > x:
                    count += 1
        # Return True if the number of intersections is odd (inside the polygon)
        return count % 2 == 1
    
"""
Segment class .
Used for computational geometry tasks, such as finding intersections.
"""      
class Segment:
    def __init__(self,origin, end):
        """
        Initialize a segment with two endpoints.
        Args:
            origin (tuple): Starting point of the segment, represented as (x, y).
            end (tuple): Ending point of the segment, represented as (x, y).
        """
        self.origin=origin
        self.end=end
        
    def intersectsSegment(self,other):
        """
        Determines whether this segment intersects with another segment.
        Args:
           other (Segment): Another Segment object.
        Returns:
           tuple: (bool, float) where the boolean indicates intersection, and the float is the parameter t
                  of this segment at the intersection point.
        """
        
        # equation of the first segment: (x1 + t*(u1-x1), y1 + t*(v1-y1)), t in [0,1]
        (x1,y1)=self.origin
        (u1,v1)=self.end
       
        # equation of the second segment: (x2 + s*(u2-x2), y2 + s*(v2-y2)), s in [0,1]
        (x2,y2)=other.origin
        (u2,v2)=other.end
        
        # Create linear system to solve for parameters t and s
        # x1 + t*(u1-x1) = x2 + s*(u2-x2) <-> (u1-x1)*t + (x2-u2)*s = x2-x1
        # y1 + t*(v1-y1) = y2 + s*(v2-y2) <-> (v1-y1)*t + (y2-v2)*s = y2-y1
        a = np.array([[u1-x1, x2-u2], [v1-y1, y2-v2]])
        b = np.array([x2-x1, y2-y1])
        try:
            x=np.linalg.solve(a,b)
        except:
            return (False, 0) # No solution or parallel lines
        else:
            t, s= x
            if 0 <= t <= 1 and 0 <= s < 1:
                return (True, t) 
            else:
                return (False, 0)
            
    def intersectsPolygon(self,polygon):
        """
         Checks if this segment intersects with any edge of a given polygon.
         Args:
             polygon (Polygon): The polygon to check for intersection.
         Returns:
             tuple: (bool, Segment, float, int) where the boolean indicates intersection, the Segment is
                    the closest intersecting edge, the float is the parameter t of this segment at the closest
                    intersection point, and the int is the number of intersections.
        """
        minimum=2  # Minimum t value for intersections (default to a value higher than 1)
        positive_minimum=2 # Minimum positive t value
        count=0 # Total number of intersections
        segment_int=[] # Closest intersecting edge
        segment_int_positive=[] # Closest intersecting edge with t > 0
        
        for segment in polygon.edges():
            intersection = self.intersectsSegment(segment)
            if intersection[0]: # If an intersection is found
                count=count+1
                t = intersection[1]
                if t < minimum:
                    minimum = t
                    segment_int=segment
                if 0 < t < positive_minimum:
                    positive_minimum = t
                    segment_int_positive=segment
                    
        # Handle no intersections            
        if count==0: 
            return (False,[], minimum, count)
        
        # Handle cases based on the number and type of intersections
        if minimum==0:
            if count==1:
                return (False,[], minimum, count)
            elif count % 2 == 0:
                s=Segment(self.origin, segment_int.end)
                (boolean, seg, m, c) = s.intersectsPolygon(polygon)
                if c>1:
                    return (True, seg, minimum, count)
                else:
                    return (True,segment_int, minimum, count)
            else:
                s=Segment(self.origin, segment_int_positive.end)
                (boolean, seg, m, c) = s.intersectsPolygon(polygon)
                if c>1:
                    return (True,seg, minimum, count)
                else:
                    return (True,segment_int_positive, minimum, count)
        else:
            s=Segment(self.origin, segment_int.end)
            (boolean, seg, m, c) = s.intersectsPolygon(polygon)
            if c>1:
                return (True,seg, minimum, count)
            else:
                return (True,segment_int, minimum, count)

#Counterclockwise rotation that fixes the y-axis  (rotation in the xz-plane)
def rot_fix_y(angle_rad,x,y,z):
    """
    Rotates a point counterclockwise around the y-axis by a given angle in radians.
    Args:
        angle_rad (float): The rotation angle in radians.
        x, y, z (float): The coordinates of the point to rotate.
    Returns:
        tuple: The rotated coordinates (x', y', z').
    """
    return (x*math.cos(angle_rad)-z*math.sin(angle_rad),
            y,
            x*math.sin(angle_rad)+z*math.cos(angle_rad))

# Counterclockwise rotation that fixes the z-axis (rotation in the xy-plane)
def rot_fix_z(angle_rad,x,y,z):
    """
    Rotates a point counterclockwise around the z-axis by a given angle in radians.
    Args:
        angle_rad (float): The rotation angle in radians.
        x, y, z (float): The coordinates of the point to rotate.
    Returns:
        tuple: The rotated coordinates (x', y', z').
    """
    return (x*math.cos(angle_rad)-y*math.sin(angle_rad),
            x*math.sin(angle_rad)+y*math.cos(angle_rad),
            z)

#Transform radians to degrees
def radian2degree(radians):
    """
    Converts an angle from radians to degrees.
    Args:
        radians (float): The angle in radians.
    Returns:
        float: The angle in degrees.
    """
    return radians*180/math.pi

#Transform degrees to radians
def degree2radian(degrees):
    """
    Converts an angle from degrees to radians.
    Args:
        degrees (float): The angle in degrees.
    Returns:
        float: The angle in radians.
    """
    return degrees*math.pi/180

# Stereographic projection of a 3D point onto the x=1 plane. The pole is (-1,0,0)
def stereographic_projection(x,y,z):
    """
    Projects a 3D point onto a 2D plane, namely x=1, using stereographic projection with (-1,0,0) as the pole.
    Args:
        x, y, z (float): The coordinates of the point to project.
    Returns:
        tuple: The projected coordinates
    """
    return (2*y/(x+1), 2*z/(x+1))
    
# Inverse stereographic projection from the x=1 plane to the |(x,y,z)|=1
def inverse_stereographic_projection(u,v):
    """
    Converts 2D radial projection coordinates back to 3D coordinates.
    Args:
        u, v (float): The radial projection coordinates.
    Returns:
        tuple: The 3D coordinates (x, y, z)
    """
    t=4/(u**2+v**2+4)
    return (2*t-1, t*u, t*v)

# Radial projection of a 3D point onto the x=1 plane
def radial_projection(x,y,z):
    """
    Projects a 3D point onto a 2D plane, namely x=1, using radial projection.
    Args:
        x, y, z (float): The coordinates of the point to project.
    Returns:
        tuple: The projected coordinates
    """
    if x>0:
        return (y/x, z/x)
    
# Inverse radial projection from the x=1 plane to the |(x,y,z)|=1, x>0 hemisphere   
def inverse_radial_projection(u,v):
    """
    Converts 2D radial projection coordinates back to 3D coordinates.
    Args:
        u, v (float): The radial projection coordinates.
    Returns:
        tuple: The 3D coordinates (x, y, z), where x > 0.
    """
    r=math.sqrt(1+u**2+v**2)
    return (1/r, u/r, v/r)

# Parallel projection of a 3D point onto the yz-plane    
def parallel_projection(x,y,z):
    """
    Projects a 3D point onto the yz-plane using parallel projection.
    Args:
        x, y, z (float): The coordinates of the point to project.
    Returns:
        tuple: The projected coordinates (y, z).
    """
    if x>0:
        return (y, z)
    
# Inverse parallel projection from the x=1 plane to the |(x,y,z)|=1, x>0 hemisphere   
def inverse_parallel_projection(u,v):
    """
    Converts 2D parallel projection coordinates back to 3D coordinates.
    Args:
        u, v (float): The parallel projection coordinates.
    Returns:
        tuple: The 3D coordinates (x, y, z), where x > 0.
    """
    return (math.sqrt(1-u**2-v**2), u, v)

# Convert polar coordinates to Cartesian coordinates
def polar2cartesian(r, alpha):
    """
    Converts polar coordinates to Cartesian coordinates.
    Args:
        r (float): The radial distance.
        alpha (float): The angle in radians.
    Returns:
        tuple: The Cartesian coordinates (x, y).
    """
    return (r*math.cos(alpha), r*math.sin(alpha))

# Convert Cartesian coordinates to polar coordinates
def cartesian2polar(x,y):
    """
    Converts Cartesian coordinates to polar coordinates.
    Args:
        x, y (float): The Cartesian coordinates.
    Returns:
        tuple: The polar coordinates (r, alpha), where alpha is in radians.
    """
    r = math.sqrt(x**2+y**2) # Radial distance
    alpha = float(np.arctan2(y, x)) # Angle in radians
    return(r, alpha)

# Convert spherical coordinates to Cartesian coordinates
def spherical2cartesian(r, latitude, longitude):
    """
    Converts spherical coordinates to Cartesian coordinates.
    Args:
        r (float): The radial distance.
        latitude (float): The latitude angle in radians (measured from the xy-plane).
        longitude (float): The longitude angle in radians (measured from the x-axis).
    Returns:
        tuple: The Cartesian coordinates (x, y, z).
    """
    return (r*math.cos(latitude)*math.cos(longitude),
            r*math.cos(latitude)*math.sin(longitude),
            r*math.sin(latitude))

# Convert Cartesian coordinates to spherical coordinates
def cartesian2spherical(x,y,z):
    """
    Converts Cartesian coordinates to spherical coordinates.
    Args:
        x, y, z (float): The Cartesian coordinates.
    Returns:
        tuple: The spherical coordinates (r, latitude, longitude), where:
               - r is the radial distance,
               - latitude is the angle in radians measured from the xy-plane,
               - longitude is the angle in radians measured from the x-axis.
    """
    r = math.sqrt(x**2+y**2+z**2) # Radial distance
    if r==0:
        return (0,0,0)
    latitude = math.asin(z/r) # Latitude angle
    if x==0 and y==0:
        return (r,latitude,0) # Undefined longitude, default to 0
    longitude=0
    if y>0:
        longitude = math.acos(x/math.sqrt(x**2+y**2)) # Longitude angle in the positive direction
    else:
        longitude=-math.acos(x/math.sqrt(x**2+y**2)) # Longitude angle in the negative direction
    return (r, latitude, longitude)




