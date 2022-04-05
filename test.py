import math
from scipy.optimize import minimize

class Point():
    def __init__(self, lat, lon):
        self.lat = lat
        self.lon = lon

    def __repr__(self) -> str:
        return f'[{self.lat}, {self.lon}]'

class Beacon:
    def __init__(self, point, limits):
        self.point = point
        self.limits = limits

def get_point_at_distance_and_bearing(point, distance, bearing):
    R = 6378100 #Radius of the Earth
    brng = math.radians(bearing) #Bearing is 90 degrees converted to radians.
    d = distance #Distance in km

    #lat2  52.20444 - the lat result I'm hoping for
    #lon2  0.36056 - the long result I'm hoping for.

    lat1 = math.radians(point.lat) #Current lat point converted to radians
    lon1 = math.radians(point.lon) #Current long point converted to radians

    lat2 = math.asin( math.sin(lat1)*math.cos(d/R) +
        math.cos(lat1)*math.sin(d/R)*math.cos(brng))

    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),
                math.cos(d/R)-math.sin(lat1)*math.sin(lat2))

    lat2 = math.degrees(lat2)
    lon2 = math.degrees(lon2)


    return Point(lat2, lon2)

def generate_triangle_points(center, radius):
    return [get_point_at_distance_and_bearing(center, radius, bearing) for bearing in [0,120,240]]


def great_circle_distance(ptA, ptB):
    # Degrees to radians
    phi1    = math.radians(ptA.lat)
    lambda1 = math.radians(ptA.lon)

    phi2    = math.radians(ptB.lat)
    lambda2 = math.radians(ptB.lon)

    delta_lambda = math.fabs(lambda2 - lambda1)

    central_angle = \
        math.atan2 \
        (
            # Numerator
            math.sqrt
            (
                # First
                math.pow
                (
                    math.cos(phi2) * math.sin(delta_lambda)
                    , 2.0
                )
                +
                # Second
                math.pow
                (
                    math.cos(phi1) * math.sin(phi2) -
                    math.sin(phi1) * math.cos(phi2) * math.cos(delta_lambda)
                    , 2.0
                )
            ),
            # Denominator
            (
                math.sin (phi1) * math.sin(phi2) +
                math.cos (phi1) * math.cos(phi2) * math.cos(delta_lambda)
            )
        )

    R = 6371009 # m
    return R * central_angle

# Mean Square Error
# locations: [ (lat1, long1), ... ]
# distances: [ distance1, ... ]
def overlap_error(x, beacons):
    mse = 0.0
    pt = Point(x[0], x[1])
    for beacon in beacons:
        distance_calculated = great_circle_distance(pt, beacon.point)
        err = 0
        if distance_calculated < beacon.limits[0]:
            err = beacon.limits[0] - distance_calculated
        elif distance_calculated > beacon.limits[1]:
            err = distance_calculated - beacon.limits[1]
        mse += err #math.pow(err, 2.0)
    return mse / len(beacons)

def grad_func(x, beacons):
    grad = [0.0,0.0]
    for beacon in beacons:
        dist = great_circle_distance(Point(x[0], x[1]),beacon.point)
        if dist < beacon.limits[0]:
            #then gradient is magnitude 1 in the direction from the point to the beacon center
            lat_component = beacon.point.lat - x[0]
            lon_component = (beacon.point.lon - x[1]) * math.cos(math.radians(x[0]))
        elif dist > beacon.limits[1]:
            #then gradient is magnitude 1 in the direciton from the beacon center to the point
            lat_component = x[0] - beacon.point.lat
            lon_component = (x[1] - beacon.point.lon) * math.cos(math.radians(x[0]))
        else:
            lat_component = 0.0
            lon_component = 0.0
        #if we are inside the bounds, then this contributes 0 to the grad
        
        #now normalize lat_component and lon_component so that the squares sum to 1
        if lat_component != 0.0 or lon_component != 0.0:
            denom = math.sqrt(lat_component**2 + lon_component**2)
            lat_component /= denom
            lon_component /= denom

            grad[0] += lat_component
            grad[1] += lon_component
    
    grad[0] /= len(beacons)
    grad[1] /= len(beacons)
    return grad


def obfuscate_distance(dist, buckets):
    for bucket in buckets:
        if dist >= bucket[0] and dist < bucket[1]:
            return bucket

def simulate_service_distances(point, beacon_points, buckets):
    return [Beacon(bp, obfuscate_distance(great_circle_distance(bp, point), buckets)) for bp in beacon_points ]

buckets = [[0,500], [500,1000], [1000, 2000], [2000, 5000], [5000, 10000], [10000, 20000]]

center = Point(45,45)
beacon_points = generate_triangle_points(center, 1300)

print(beacon_points)
beacons = simulate_service_distances(center, beacon_points, buckets)

furthest_beacon_point = max(beacons, key=lambda b: b.limits[1]).point

#start = Point(52.20472, 0.14056)
#distance = 15
#bearing = 90 
#res = get_point_at_distance_and_bearing(start, distance, bearing)
#disp = great_circle_distance(start, res)
#err = great_circle_distance(res, Point(52.20444, 0.36056))
#print(disp, err)

# initial_location: (lat, long)
# locations: [ (lat1, long1), ... ]
# distances: [ distance1,     ... ] 
result = minimize(
    overlap_error,                         # The error function
    [furthest_beacon_point.lat, furthest_beacon_point.lon],            # The initial guess
    args=(beacons), # Additional parameters for mse
    method='L-BFGS-B',           # The optimisation algorithm
    jac=grad_func,
    options={
        'ftol':1e-5,         # Tolerance
        'maxiter': 1e+7      # Maximum iterations
    })
location = result.x
print(result)

