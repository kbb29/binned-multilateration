'''
FINDING THE INTERSECTION COORDINATES (LAT/LON) OF TWO CIRCLES (GIVEN THE COORDINATES OF THE CENTER AND THE RADII)

Many thanks to Ture Pålsson who directed me to the right source, the code below is based on whuber's brilliant logic and
explanation here https://gis.stackexchange.com/questions/48937/calculating-intersection-of-two-circles 

The idea is that;
  1. The points in question are the mutual intersections of three spheres: a sphere centered beneath location x1 (on the 
  earth's surface) of a given radius, a sphere centered beneath location x2 (on the earth's surface) of a given radius, and
  the earth itself, which is a sphere centered at O = (0,0,0) of a given radius.
  2. The intersection of each of the first two spheres with the earth's surface is a circle, which defines two planes.
  The mutual intersections of all three spheres therefore lies on the intersection of those two planes: a line.
  Consequently, the problem is reduced to intersecting a line with a sphere.

Note that "Decimal" is used to have higher precision which is important if the distance between two points are a few
meters.
'''
from decimal import Decimal
from math import cos, sin, sqrt
import math
import numpy as np
from plot import plotBeacons
from multilat_types import Point, Beacon, Limit, Centroid, BoundSquare, get_point_at_distance_and_bearing, get_bearing, great_circle_distance, generate_triangle_points


def limit_intersection(limit1, limit2):
    print(limit1.beacon.point, limit2.beacon.point)
    p1 = limit1.beacon.point.to_tuple()
    p2 = limit2.beacon.point.to_tuple()
    r1_meter = limit1.get_radius()
    r2_meter = limit2.get_radius()


    # p1 = Coordinates of Point 1: latitude, longitude. This serves as the center of circle 1. Ex: (36.110174,  -90.953524)
    # r1_meter = Radius of circle 1 in meters
    # p2 = Coordinates of Point 2: latitude, longitude. This serves as the center of circle 1. Ex: (36.110174,  -90.953524)
    # r2_meter = Radius of circle 2 in meters
    '''
    1. Convert (lat, lon) to (x,y,z) geocentric coordinates.
    As usual, because we may choose units of measurement in which the earth has a unit radius
    '''
    x_p1 = Decimal(cos(math.radians(p1[1]))*cos(math.radians(p1[0])))  # x = cos(lon)*cos(lat)
    y_p1 = Decimal(sin(math.radians(p1[1]))*cos(math.radians(p1[0])))  # y = sin(lon)*cos(lat)
    z_p1 = Decimal(sin(math.radians(p1[0])))                           # z = sin(lat)
    x1 = (x_p1, y_p1, z_p1)

    x_p2 = Decimal(cos(math.radians(p2[1]))*cos(math.radians(p2[0])))  # x = cos(lon)*cos(lat)
    y_p2 = Decimal(sin(math.radians(p2[1]))*cos(math.radians(p2[0])))  # y = sin(lon)*cos(lat)
    z_p2 = Decimal(sin(math.radians(p2[0])))                           # z = sin(lat)
    x2 = (x_p2, y_p2, z_p2)
    '''
    2. Convert the radii r1 and r2 (which are measured along the sphere) to angles along the sphere.
    By definition, one nautical mile (NM) is 1/60 degree of arc (which is pi/180 * 1/60 = 0.0002908888 radians).
    '''
    r1 = Decimal(math.radians((r1_meter/1852) / 60)) # r1_meter/1852 converts meter to Nautical mile.
    r2 = Decimal(math.radians((r2_meter/1852) / 60))
    '''
    3. The geodesic circle of radius r1 around x1 is the intersection of the earth's surface with an Euclidean sphere
    of radius sin(r1) centered at cos(r1)*x1.

    4. The plane determined by the intersection of the sphere of radius sin(r1) around cos(r1)*x1 and the earth's surface
    is perpendicular to x1 and passes through the point cos(r1)x1, whence its equation is x.x1 = cos(r1)
    (the "." represents the usual dot product); likewise for the other plane. There will be a unique point x0 on the
    intersection of those two planes that is a linear combination of x1 and x2. Writing x0 = ax1 + b*x2 the two planar
    equations are;
       cos(r1) = x.x1 = (a*x1 + b*x2).x1 = a + b*(x2.x1)
       cos(r2) = x.x2 = (a*x1 + b*x2).x2 = a*(x1.x2) + b
    Using the fact that x2.x1 = x1.x2, which I shall write as q, the solution (if it exists) is given by
       a = (cos(r1) - cos(r2)*q) / (1 - q^2),
       b = (cos(r2) - cos(r1)*q) / (1 - q^2).
    '''
    q = Decimal(np.dot(x1, x2))

    if q**2 == 1 :
        raise RuntimeError("The centers of the circles can be neither the same point nor antipodal points.")

    a = (Decimal(cos(r1)) - Decimal(cos(r2))*q) / (1 - q**2)
    b = (Decimal(cos(r2)) - Decimal(cos(r1))*q) / (1 - q**2)
    '''
    5. Now all other points on the line of intersection of the two planes differ from x0 by some multiple of a vector
    n which is mutually perpendicular to both planes. The cross product  n = x1~Cross~x2  does the job provided n is 
    nonzero: once again, this means that x1 and x2 are neither coincident nor diametrically opposite. (We need to 
    take care to compute the cross product with high precision, because it involves subtractions with a lot of
    cancellation when x1 and x2 are close to each other.)
    '''
    n = np.cross(x1, x2)
    '''
    6. Therefore, we seek up to two points of the form x0 + t*n which lie on the earth's surface: that is, their length
    equals 1. Equivalently, their squared length is 1:  
    1 = squared length = (x0 + t*n).(x0 + t*n) = x0.x0 + 2t*x0.n + t^2*n.n = x0.x0 + t^2*n.n
    '''
    x0_1 = [a*f for f in x1]
    x0_2 = [b*f for f in x2]
    x0 = [sum(f) for f in zip(x0_1, x0_2)]
    '''
        The term with x0.n disappears because x0 (being a linear combination of x1 and x2) is perpendicular to n.
        The two solutions easily are   t = sqrt((1 - x0.x0)/n.n)    and its negative. Once again high precision
        is called for, because when x1 and x2 are close, x0.x0 is very close to 1, leading to some loss of
        floating point precision.
    '''
    if (np.dot(x0, x0) > 1) or (np.dot(n,n) == 0): # This is to secure that (1 - np.dot(x0, x0)) / np.dot(n,n) > 0
        if (np.dot(n,n) == 0):
            raise RuntimeError("The centers of the circles can be neither the same point nor antipodal points.")
        else:
            raise RuntimeError("The circles do not intersect")

    print(x0, n)
    print(np.dot(x0, x0))
    print(np.dot(n, n))
    print((1 - np.dot(x0, x0)) / np.dot(n,n))
    print(sqrt((1 - np.dot(x0, x0)) / np.dot(n,n)))
    t = Decimal(sqrt((1 - np.dot(x0, x0)) / np.dot(n,n)))
    t1 = t
    t2 = -t

    i1 = x0 + t1*n
    i2 = x0 + t2*n
    '''
    7. Finally, we may convert these solutions back to (lat, lon) by converting geocentric (x,y,z) to geographic
    coordinates. For the longitude, use the generalized arctangent returning values in the range -180 to 180
    degrees (in computing applications, this function takes both x and y as arguments rather than just the
    ratio y/x; it is sometimes called "ATan2").
    '''

    i1_lat = math.degrees( math.asin(i1[2]))
    i1_lon = math.degrees( math.atan2(i1[1], i1[0] ) )
    ip1 = Point(i1_lat, i1_lon)

    i2_lat = math.degrees( math.asin(i2[2]))
    i2_lon = math.degrees( math.atan2(i2[1], i2[0] ) )
    ip2 = Point(i2_lat, i2_lon)
    return [ip1, ip2]

'''
Example: the output of below is  [(36.989311051533505, -88.15142628069133), (38.2383796094578, -92.39048549120287)]
         intersection_points = intersection((37.673442, -90.234036), 107.5*1852, (36.109997, -90.953669), 145*1852)
         print(intersection_points)
'''

if __name__ == '__main__':
    b1 = Beacon(Point(37.673442, -90.234036), (199090.0, 199090 + 90000))
    b2 = Beacon(Point(36.109997, -90.953669), (268540, 268540 + 90000))
    l1 = Limit(b1, True)
    l2 = Limit(b2, True)
    #intersection_points = intersection((37.673442, -90.234036), 107.5*1852, (36.109997, -90.953669), 145*1852)
    intersection_points = limit_intersection(l1, l2)
    print(intersection_points)
    plotBeacons([b1, b2], preds=intersection_points)

    b1 = Beacon(Point(37.673442, -90.234036), (1000.0, 2000))
    b2 = Beacon(get_point_at_distance_and_bearing(b1.point, 2000.0, 45.0), (2000.0, 5000))
    l1 = Limit(b1, True)
    l2 = Limit(b2, True)
    #intersection_points = intersection((37.673442, -90.234036), 107.5*1852, (36.109997, -90.953669), 145*1852)
    intersection_points = limit_intersection(l1, l2)
    plotBeacons([b1, b2], preds=intersection_points)
    
