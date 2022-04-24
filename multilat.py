import math
from scipy.optimize import minimize
from typing import Tuple, List
from plot import plotBeacons
from geographiclib.geodesic import Geodesic

def get_bearing(pt1, pt2):
    bearing = Geodesic.WGS84.Inverse(pt1.lat, pt1.lon, pt2.lat, pt2.lon)['azi1']
    if bearing < 0:
        bearing += 360
    return bearing

class Point():
    def __init__(self, lat, lon):
        self.lat = lat
        self.lon = lon
    
    def to_tuple(self):
        return (self.lat, self.lon)

    def __repr__(self) -> str:
        return f'[{self.lat}, {self.lon}]'

class Beacon:
    def __init__(self, point, limits):
        self.point = point
        self.limits = limits

    def is_within_limits(self, point):
        dist = great_circle_distance(self.point, point)
        return dist >= self.limits[0] and dist <= self.limits[1]

class Centroid:
    def __init__(self, point, radius):
        self.point = point
        self.radius = radius
        self.n = 0

    def add_point(self, point):
        self.point.lat += point.lat
        self.point.lon += point.lon
        self.n += 1

    def divide(self):
        self.point.lat /= self.n
        self.point.lon /= self.n
    
    def set_radius(self, points):
        self.radius = max([great_circle_distance(p, self.point) for p in points])

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


def error_component_regular(beacon, pt, zero_gradient_within_limits=True):
        distance_calculated = great_circle_distance(pt, beacon.point)
        err = 0.0
        if distance_calculated < beacon.limits[0]:
            err = beacon.limits[0] - distance_calculated
        elif distance_calculated > beacon.limits[1]:
            err = distance_calculated - beacon.limits[1]
        #we are within the bounds. apply a gentle gradient to nudge the optimizer towards the middle of the bounds
        elif not zero_gradient_within_limits:
            #set the gradient to 0.2 towards the middle
            #error at bounds = 0
            #error at the center of the bounds = 0.2 * (upper - lower) / 2
            #error everywhere = dist_to_nearest_limit * 0.2
            #dist_to_nearest_limit = min(abs(distance_calculated - beacon.limits[0]), abs(beacon.limits[1] - distance_calculated))
            limits_center = (beacon.limits[1] + beacon.limits[0]) / 2
            dist_from_limits_center = abs(distance_calculated - limits_center)
            #err = dist_to_nearest_limit * -0.2
            d = beacon.limits[1] - limits_center
            a = 1 / (2*d)
            c = -1 * d / 2
            err = a * dist_from_limits_center**2 + c

        return err

def error_component_maximize(beacon, pt):
        distance_calculated = great_circle_distance(pt, beacon.point)
        err = 0.0
        if distance_calculated < beacon.limits[1]:
            err = 0.5 * (beacon.limits[1] - distance_calculated)
        elif distance_calculated > beacon.limits[1]:
            err = distance_calculated - beacon.limits[1]
        return err

def error_component_minimize(beacon, pt):
        distance_calculated = great_circle_distance(pt, beacon.point)
        err = 0.0
        if distance_calculated < beacon.limits[0]:
            err = beacon.limits[0] - distance_calculated
        elif distance_calculated > beacon.limits[0]:
            err = 0.5 * (distance_calculated - beacon.limits[0])
        return err

# Mean Square Error
# locations: [ (lat1, long1), ... ]
# distances: [ distance1, ... ]
def overlap_error(x, beacons, beacons_maximize, beacons_minimize, zero_gradient_within_limits) -> float:
    mse = 0.0
    pt = Point(*x)
    errs = [error_component_regular(beacon, pt, zero_gradient_within_limits=zero_gradient_within_limits) for beacon in beacons]
    errs += [error_component_maximize(beacon, pt) for beacon in beacons_maximize]
    errs += [error_component_minimize(beacon, pt) for beacon in beacons_minimize]
    num_beacons = len(beacons) + len(beacons_maximize) + len(beacons_minimize)
    return sum(errs) / num_beacons

#for when we want to optimize the distance only
def overlap_error_by_distance_from_beacon(x, beacons, target_idx, bearing, minimize=False):
    pt = get_point_at_distance_and_bearing(beacons[target_idx].point, x[0], bearing)
    total = 0.0
    for i,beacon in enumerate(beacons):
        if i == target_idx: 
            if minimize:
                total += error_component_minimize(beacon, pt)
            else:
                total += error_component_maximize(beacon, pt)
        else:
            total += error_component_regular(beacon, pt, gradient_within_limits=False)
    return total

def grad_component_regular(beacon, pt, zero_gradient_within_limits=True) -> Tuple[float,float]:
        dist = great_circle_distance(pt,beacon.point)
        if dist < beacon.limits[0]:
            #then gradient is magnitude 1 in the direction from the point to the beacon center
            lat_component = beacon.point.lat - pt.lat
            lon_component = (beacon.point.lon - pt.lon) * math.cos(math.radians(pt.lat))
            return grad_component_norm(lat_component, lon_component, norm_to=1)
        elif dist > beacon.limits[1]:
            #then gradient is magnitude 1 in the direciton from the beacon center to the point
            lat_component = pt.lat - beacon.point.lat
            lon_component = (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))
            return grad_component_norm(lat_component, lon_component, norm_to=1)
        elif zero_gradient_within_limits:
            return 0.0, 0.0
        else:
            #we are within the bounds. apply a gentle gradient to nudge the optimizer towards the middle of the bounds
            #set the gradient to 0.2 towards the middle
            #here we are closer to the beacon than the middle, so grad is mag 0.2 in the direction from the point to the beacon center
            limits_center = (beacon.limits[1] + beacon.limits[0]) / 2
            dist_from_limits_center = dist - limits_center
            #err = dist_to_nearest_limit * -0.2
            d = beacon.limits[1] - limits_center
            a = 1 / d
            mag = a * dist_from_limits_center

            #now determine the direction of the gradient based on whether it is beyond the center of the limits or not
            if dist < limits_center:
                #if the point is closer to the beacon than the limits_center, then the grad is increasing in direction from point to the beacon
                lat_component = beacon.point.lat - pt.lat
                lon_component = (beacon.point.lon - pt.lon) * math.cos(math.radians(pt.lat))
            else:
                #if the point is beyond the limits_center, then the grad is increasing in direction from beacon to the point
                lat_component = pt.lat - beacon.point.lat
                lon_component = (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))
            
            return grad_component_norm(lat_component, lon_component, norm_to=mag)
            

def grad_component_maximize(beacon, pt) -> Tuple[float,float]:
        dist = great_circle_distance(pt,beacon.point)
        if dist < beacon.limits[1]:
            #then gradient is magnitude 0.5 in the direction from the point to the beacon center
            lat_component = (beacon.point.lat - pt.lat)
            lon_component = (beacon.point.lon - pt.lon) * math.cos(math.radians(pt.lat))
            return grad_component_norm(lat_component, lon_component, norm_to=0.5)
        elif dist > beacon.limits[1]:
            #then gradient is magnitude 1 in the direction from the beacon center to the point
            lat_component = pt.lat - beacon.point.lat
            lon_component = (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))
            return grad_component_norm(lat_component, lon_component, norm_to=1)
        
        #if we are inside the bounds, then this contributes 0 to the grad
        return 0.0, 0.0

def grad_component_minimize(beacon, pt) -> Tuple[float,float]:
        dist = great_circle_distance(pt,beacon.point)
        if dist < beacon.limits[0]:
            #then gradient is magnitude 1 in the direction from the point to the beacon center
            lat_component = beacon.point.lat - pt.lat
            lon_component = (beacon.point.lon - pt.lon) * math.cos(math.radians(pt.lat))
            return grad_component_norm(lat_component, lon_component, norm_to=1)
        elif dist > beacon.limits[0]:
            #then gradient is magnitude 0.5 in the direction from the beacon center to the point
            lat_component = (pt.lat - beacon.point.lat)
            lon_component = (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))
            return grad_component_norm(lat_component, lon_component, norm_to=0.5)
        
        return 0.0, 0.0
        #if we are inside the bounds, then this contributes 0 to the grad

def grad_component_norm(lat_component, lon_component, norm_to=1):
    if lat_component != 0.0 or lon_component != 0.0:
        denom = math.sqrt((lat_component**2 + lon_component**2)/norm_to**2)
        return lat_component/denom, lon_component/denom
    return 0.0, 0.0


def grad_func(x, beacons, maximize_beacons, minimize_beacons, zero_gradient_within_limits) -> List[float]:
    grad = [0.0,0.0]
    pt = Point(*x)
    component_tuples = [grad_component_regular(beacon, pt, zero_gradient_within_limits=zero_gradient_within_limits) for beacon in beacons]
    component_tuples += [grad_component_maximize(beacon, pt) for beacon in maximize_beacons]
    component_tuples += [grad_component_minimize(beacon, pt) for beacon in minimize_beacons]

    num_total_beacons = len(beacons) + len(minimize_beacons) + len(maximize_beacons)
    grad[0] = sum([latc for latc,_ in component_tuples]) / num_total_beacons
    grad[1] = sum([lonc for _,lonc in component_tuples]) / num_total_beacons

    return grad


def obfuscate_distance(dist, buckets):
    for bucket in buckets:
        if dist >= bucket[0] and dist < bucket[1]:
            return bucket

def simulate_service_distances(point, beacon_points, buckets):
    return [Beacon(bp, obfuscate_distance(great_circle_distance(bp, point), buckets)) for bp in beacon_points ]

#print(beacon_points)


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

def find_outer_bounds_for_beacon(beacons, reference_point, target_idx):
    bounds = []
    optimize_beacon = beacons[target_idx]
    reference_bearing = get_bearing(optimize_beacon.point, reference_point)
    standard_beacons = beacons[:target_idx] + beacons[target_idx+1:]
    #min bound
    start_point = get_point_at_distance_and_bearing(optimize_beacon.point, optimize_beacon.limits[0] * 0.5, reference_bearing)
    result = minimize(
        overlap_error,                         # The error function
        [start_point.lat, start_point.lon],            # The initial guess
        args=(standard_beacons, [], [optimize_beacon], True), # Additional parameters for mse
        method='L-BFGS-B',           # The optimisation algorithm
        #jac=grad_func,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })

    print("minimize minimize\n")
    print(result)
    bounds.append(Point(*result.x))

    #max bound
    start_point = get_point_at_distance_and_bearing(optimize_beacon.point, optimize_beacon.limits[1] * 1.5, reference_bearing)
    result = minimize(
        overlap_error,                         # The error function
        [start_point.lat, start_point.lon],            # The initial guess
        args=(standard_beacons, [optimize_beacon], [], True), # Additional parameters for mse
        method='L-BFGS-B',           # The optimisation algorithm
        #jac=grad_func,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })

    print("minimize maximize\n")
    print(result)
    bounds.append(Point(*result.x))
    return bounds

def find_outer_bounds_for_beacon_alt(beacons, reference_point, target_idx):
    bounds = []
    bearing = get_bearing(beacons[target_idx].point, reference_point)
    result = minimize(
        overlap_error_by_distance_from_beacon,                         # The error function
        [beacons[target_idx].limits[0]],            # The initial guess
        args=(beacons, target_idx, bearing, True), # Additional parameters for mse
        method='L-BFGS-B',           # The optimisation algorithm
        #jac=grad_func,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })
    print("minimize minimize\n")
    print(result)
    bounds.append(get_point_at_distance_and_bearing(beacons[target_idx].point, result.x, bearing))
    result = minimize(
        overlap_error_by_distance_from_beacon,                         # The error function
        [beacons[target_idx].limits[0]],            # The initial guess
        args=(beacons, target_idx, bearing, False), # Additional parameters for mse
        method='L-BFGS-B',           # The optimisation algorithm
        #jac=grad_func,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })

    print("maximize minimize\n")
    print(result)
    bounds.append(get_point_at_distance_and_bearing(beacons[target_idx].point, result.x, bearing))
    return bounds


def calculate_furthest_point(x, points):
    c = Centroid(Point(*x), 0)
    c.set_radius(points)
    return c.radius


def do_multilat(beacons, tag=''):
    furthest_beacon_point = max(beacons, key=lambda b: b.limits[1]).point
    #get a point in the intersection

    result = minimize(
        overlap_error,
        [furthest_beacon_point.lat, furthest_beacon_point.lon],
        args=(beacons, [], [], False),
        method='L-BFGS-B',           # The optimisation algorithm
        jac=grad_func,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })

    print("initial minimize\n", result)
    intersection_point = Point(*result.x)

    if not all([b.is_within_limits(intersection_point) for b in beacons]):
        raise RuntimeError(f'the first minimization did not settle within the intersection {result}')

    #bounds  = [find_outer_bound_for_beacon(beacons, intersection_point, optimize_index=i, maximize=True ) for i in range(len(beacons))]
    #bounds += [find_outer_bound_for_beacon(beacons, intersection_point, optimize_index=i, maximize=False) for i in range(len(beacons))]
    bounds = [pt for bounds in [find_outer_bounds_for_beacon(beacons, intersection_point, i) for i in range(len(beacons))] for pt in bounds]
    #print(bounds)
    plotBeacons(beacons, actual=intersection_point, preds=bounds, options={'tag': f'debug-{tag}'})
    centroid = Centroid(Point(0.0, 0.0), 0.0)
    [centroid.add_point(p) for p in bounds]
    centroid.divide()

    result = minimize(calculate_furthest_point, [centroid.point.lat ,centroid.point.lon], bounds, method='L-BFGS-B', options={'ftol': 1e-5, 'maxiter': 1e6})
    return Centroid(Point(*result.x), result.fun), bounds

if __name__ == '__main__':

    buckets = [[0,500], [500,1000], [1000, 2000], [2000, 5000], [5000, 10000], [10000, 20000]]
    center = Point(45,45)
    beacon_points = generate_triangle_points(center, 1300)
    for n,actual in enumerate([center, get_point_at_distance_and_bearing(center, 1300, 300)]):
        beacons = simulate_service_distances(actual, beacon_points, buckets)
        centroid, bounds = do_multilat(beacons, tag=str(n))

        plotBeacons(beacons, actual=actual, preds=bounds, centroid=centroid, options={'tag': str(n)})
