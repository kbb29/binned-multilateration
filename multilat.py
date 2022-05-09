import math
from matplotlib.pyplot import magnitude_spectrum, plot
from scipy.optimize import minimize
from typing import Tuple, List
from plot import plotBeacons
from multilat_types import Point, Beacon, Limit, Centroid, BoundSquare, IntersectionPoint, get_point_at_distance_and_bearing, get_bearing, great_circle_distance, generate_triangle_points, normalize_bearing
from limit_intersection import limit_intersection
from itertools import combinations

class IntersectionError(RuntimeError):
    def __init__(self, message, result):
        self.result = result
        self.message = message
        super().__init__(self.message)


def normalize_gradient_components(lat_component, lon_component, norm_to=1):
    if lat_component != 0.0 or lon_component != 0.0:
        denom = math.sqrt((lat_component**2 + lon_component**2)/norm_to**2)
        return lat_component/denom, lon_component/denom
    return 0.0, 0.0

def error_component_regular_quadratic(beacon, pt):
        distance_calculated = great_circle_distance(pt, beacon.point)
        err = 0.0
        limits_center = (beacon.limits[1] + beacon.limits[0]) / 2
        dist_from_limits_center = abs(distance_calculated - limits_center)
        d = beacon.limits[1] - limits_center
        if distance_calculated < beacon.limits[0] or distance_calculated > beacon.limits[1]:
            a = 1 / (2*d)
            c = d / -2
            err = a * dist_from_limits_center**2 + c
        #we are within the bounds. apply a gentle gradient to nudge the optimizer towards the middle of the bounds
        else:
            #set the gradient to 0.2 towards the middle
            #error at bounds = 0
            #error at the center of the bounds = 0.2 * (upper - lower) / 2
            #error everywhere = dist_to_nearest_limit * 0.2
            a = 1 / (4*d**3)
            c = d / -4
            err = a * dist_from_limits_center**4 + c
          
        return err


def error_component_regular_linear(beacon, pt):
        distance_calculated = great_circle_distance(pt, beacon.point)
        err = 0.0
        if distance_calculated < beacon.limits[0]:
            err = beacon.limits[0] - distance_calculated
        elif distance_calculated > beacon.limits[1]:
            err = distance_calculated - beacon.limits[1]
        #we are within the bounds. apply a gentle gradient to nudge the optimizer towards the middle of the bounds
        else:
            #set the gradient to 0.2 towards the middle
            #error at bounds = 0
            #error at the center of the bounds = 0.2 * (upper - lower) / 2
            #error everywhere = dist_to_nearest_limit * 0.2
            #dist_to_nearest_limit = min(abs(distance_calculated - beacon.limits[0]), abs(beacon.limits[1] - distance_calculated))
            limits_center = (beacon.limits[1] + beacon.limits[0]) / 2
            dist_from_limits_center = abs(distance_calculated - limits_center)
            #err = dist_to_nearest_limit * -0.2
            d = beacon.limits[1] - limits_center
            a = 1 / (4*d**3)
            c = d / -4
            err = a * dist_from_limits_center**4 + c
            #a = 1 / (2*d)
            #c = d / -2
            #err = a * dist_from_limits_center**2 + c

        return err

def error_component_find_bounds_optimize(beacon, pt, maximize):
    dist = great_circle_distance(pt, beacon.point)
    err = 0.0
    if maximize:
        if dist > beacon.limits[1]:
            err = (dist - beacon.limits[1])**2
        else:
            err = (beacon.limits[1] - dist)**1.2
    else:
        if dist < beacon.limits[0]:
            err = (beacon.limits[0] - dist)**2
        else:
            err = (dist - beacon.limits[0])**1.2
    return err


def error_component_find_bounds_other(beacon, pt):
    dist = great_circle_distance(pt, beacon.point)
    err = 0.0
    if dist > beacon.limits[1]:
        err = (dist - beacon.limits[1])**2
    elif dist < beacon.limits[0]:
        err = (beacon.limits[0] - dist)**2
    #we are within the bounds. 
    else:
        err = 0.0
    return err



# Mean Square Error
# locations: [ (lat1, long1), ... ]
# distances: [ distance1, ... ]
def overlap_error(x, beacons) -> float:
    mse = 0.0
    pt = Point(*x)
    errs = [error_component_regular_quadratic(beacon, pt) for beacon in beacons]
    num_beacons = len(beacons) 
    return sum(errs) / num_beacons

def overlap_error_find_bounds(x, beacons, reference_point, avg_reference_dist, maximize) -> float:
    pt = Point(*x)
    errs = [error_component_find_bounds_other(beacon, pt) for beacon in beacons]
    dist_from_reference = great_circle_distance(pt, reference_point)
    if maximize:
        errs.append(avg_reference_dist - dist_from_reference )
    else:
        errs.append(dist_from_reference - avg_reference_dist)
    return sum(errs) / (len(beacons) + 1)

def grad_component_regular(beacon, pt) -> Tuple[float,float]:
        dist = great_circle_distance(pt,beacon.point)
        if dist < beacon.limits[0]:
            #then gradient is magnitude 1 in the direction from the point to the beacon center
            direction = -1
            magnitude = 1
        elif dist > beacon.limits[1]:
            #then gradient is magnitude 1 in the direciton from the beacon center to the point
            direction = 1
            magnitude = 1
        else:
            #we are within the bounds. apply a gentle gradient to nudge the optimizer towards the middle of the bounds
            #set the gradient to 0.2 towards the middle
            #here we are closer to the beacon than the middle, so grad is mag 0.2 in the direction from the point to the beacon center
            limits_center = (beacon.limits[1] + beacon.limits[0]) / 2
            dist_from_limits_center = dist - limits_center
            #err = dist_to_nearest_limit * -0.2
            d = beacon.limits[1] - limits_center
            a = 1 / d
            magnitude = a * dist_from_limits_center

            #now determine the direction of the gradient based on whether it is beyond the center of the limits or not
            if dist < limits_center:
                #if the point is closer to the beacon than the limits_center, then the grad is increasing in direction from point to the beacon
                lat_component = beacon.point.lat - pt.lat
                lon_component = (beacon.point.lon - pt.lon) * math.cos(math.radians(pt.lat))
                direction = -1
            else:
                #if the point is beyond the limits_center, then the grad is increasing in direction from beacon to the point
                direction = 1

        lat_component = direction * (pt.lat - beacon.point.lat)
        lon_component = direction * (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))
            
        return normalize_gradient_components(lat_component, lon_component, norm_to=magnitude)
  
def grad_component_find_bounds_optimize(beacon, pt, maximize) -> Tuple[float, float]:
    dist = great_circle_distance(pt, beacon.point)
    if maximize:
        if dist > beacon.limits[1]:
            magnitude = 2*(dist - beacon.limits[1])
            direction = 1
        else:
            magnitude = 1.2*(beacon.limits[1] - dist)**0.2
            direction = -1
    else:
        if dist < beacon.limits[0]:
            magnitude = 2*(beacon.limits[0] - dist)
            direction = -1
        else:
            magnitude = 1.2*(dist - beacon.limits[0])**0.2
            direction = 1
        
    lat_component = direction * (pt.lat - beacon.point.lat)
    lon_component = direction * (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))

    return normalize_gradient_components(lat_component, lon_component, norm_to=magnitude)

def grad_component_find_bounds_other(beacon, pt) -> Tuple[float, float]:
    dist = great_circle_distance(pt, beacon.point)
    if dist > beacon.limits[1]:
        magnitude = 2*(dist - beacon.limits[1])
        direction = 1
    elif dist < beacon.limits[0]:
        magnitude = 2*(beacon.limits[0] - dist)
        direction = -1
    else:
        return 0.0, 0.0

    lat_component = direction * (pt.lat - beacon.point.lat)
    lon_component = direction * (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))

    return normalize_gradient_components(lat_component, lon_component, norm_to=magnitude)

def get_grad_components(fm, to, magnitude=1):
    lat_component = to.lat - fm.lat
    lon_component = (to.lon - fm.lon) * math.cos(math.radians(to.lat))
    return normalize_gradient_components(lat_component, lon_component, norm_to=magnitude)


def grad_func(x, beacons) -> List[float]:
    grad = [0.0,0.0]
    pt = Point(*x)
    component_tuples = [grad_component_regular(beacon, pt) for beacon in beacons]

    grad[0] = sum([latc for latc,_ in component_tuples]) / len(beacons)
    grad[1] = sum([lonc for _,lonc in component_tuples]) / len(beacons)

    return grad

def grad_func_find_bounds(x, beacons, reference_point, maximize) -> List[float]:
    grad = [0.0, 0.0]
    pt = Point(*x)
    component_tuples = [grad_component_find_bounds_other(beacon, pt) for beacon in beacons]
    if maximize:
        component_tuples.append( get_grad_components(pt, reference_point, magnitude=1) )
    else:
        component_tuples.append( get_grad_components(reference_point, pt, magnitude=1))

    grad[0] = sum([latc for latc,_ in component_tuples]) / len(beacons)
    grad[1] = sum([lonc for _,lonc in component_tuples]) / len(beacons)

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

def find_outer_bounds_for_reference_point(beacons, intersection_point, reference_point):
    bounds = []
    avg_reference_dist = great_circle_distance(intersection_point, reference_point)
    #reference_bearing = get_bearing(optimize_beacon.point, reference_point)
    #min bound
    #start_point = get_point_at_distance_and_bearing(optimize_beacon.point, optimize_beacon.limits[0] * 0.5, reference_bearing)
    start_point = intersection_point
    result = minimize(
        overlap_error_find_bounds,         # The error function
        start_point.to_tuple(),            # The initial guess
        args=(beacons, reference_point, avg_reference_dist,  False), # Additional parameters for mse
        method='L-BFGS-B',           # The optimisation algorithm
        #jac=grad_func_find_bounds,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })

    #print("minimize minimize\n")
    #print(result)
    #bounds.append(Point(*result.x))

    #max bound
    #start_point = get_point_at_distance_and_bearing(optimize_beacon.point, optimize_beacon.limits[1] * 1.5, reference_bearing)
    #start_point = reference_point
    result = minimize(
        overlap_error_find_bounds,         # The error function
        start_point.to_tuple(),            # The initial guess
        args=(beacons, reference_point, avg_reference_dist, True), # Additional parameters for mse
        method='L-BFGS-B',           # The optimisation algorithm
        #jac=grad_func_find_bounds,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })

    print("minimize maximize\n")
    print(result)
    bounds.append(Point(*result.x))
    return bounds

    

def find_closest_limit_to_intersection(beacons, intersection_point):
    closest_limit_with_distance = min([b.get_nearest_limit_with_distance_from_point(intersection_point) for b in beacons], key=lambda l: l[1])
    closest_limit = closest_limit_with_distance[0]
    return closest_limit

def find_next_boundary_point_anticlockwise(beacons, current_boundary: Limit, current_boundary_point: Point, tag=None) -> IntersectionPoint:
    bearing_of_current_boundary_point = get_bearing(current_boundary.beacon.point, current_boundary_point)
    other_beacons = [b for b in beacons if b != current_boundary.beacon]
    other_beacon_limits = [l for b in other_beacons for l in b.get_limits()]
    intersection_points = [ip for l in other_beacon_limits for ip in limit_intersection(current_boundary, l) if ip.point != current_boundary_point]
    intersection_points_with_bearing = [(ipt, normalize_bearing(current_boundary.get_bearing_of_point(ipt.point) - bearing_of_current_boundary_point) ) for ipt in intersection_points]
    #plotBeacons(beacons, preds=[ip.point for ip in intersection_points], actual=current_boundary_point, options={'tag': str(tag)})
    
    if current_boundary.limit_index == 0:
        closest_ip = min(intersection_points_with_bearing, key=lambda ptb: ptb[1])[0]
    else:
        closest_ip = max(intersection_points_with_bearing, key=lambda ptb: ptb[1])[0]
    return closest_ip

def map_bounds_of_intersection(beacons, intersection_point, tag=''):
    closest_limit = find_closest_limit_to_intersection(beacons, intersection_point)
    bearing_of_closest_limit_point = get_bearing(closest_limit.beacon.point, intersection_point)
    closest_limit_point = get_point_at_distance_and_bearing(closest_limit.beacon.point, closest_limit.get_radius(), bearing_of_closest_limit_point)
    
    first_ip = find_next_boundary_point_anticlockwise(beacons, closest_limit, closest_limit_point)
    current_ip = first_ip
    current_limit = current_ip.limit1 if current_ip.limit2 == closest_limit else current_ip.limit2
        
    bounds = []
    while True:
        bounds.append(current_ip.point)
        if len(bounds) > len(beacons)**2 * 4: # this is a safety net
            raise RuntimeError('boundary mapping did not return to start point')
        next_ip = find_next_boundary_point_anticlockwise(beacons, current_limit, current_ip.point, tag=str(len(bounds)))
        if next_ip == first_ip: # we are back at the start, break
            break
        #if len(bounds) > 3: break
        next_limit = next_ip.limit1 if next_ip.limit2 == current_limit else next_ip.limit2

        current_ip = next_ip
        current_limit = next_limit
        
    #plotBeacons(beacons, actual=intersection_point, preds=[closest_limit_point] + bounds, options={'tag': str(tag)})
    return bounds

def compute_intersection_points_for_beacons(beacon_pair):
    (beacon1, beacon2) = beacon_pair
    ips  = limit_intersection( beacon1.get_lower_limit(), beacon2.get_upper_limit() )
    ips += limit_intersection( beacon1.get_lower_limit(), beacon2.get_lower_limit() )
    ips += limit_intersection( beacon1.get_upper_limit(), beacon2.get_upper_limit() )
    ips += limit_intersection( beacon1.get_upper_limit(), beacon2.get_lower_limit() )
    return ips

def flatten_generator(gen):
    for thing in gen:
        for subthing in thing:
            yield subthing

def compute_limit_intersections(beacons: List[Beacon]): # -> List[IntersectionPoint]:
    all_beacon_combinations = combinations(beacons, 2)
    all_intersection_points = map(compute_intersection_points_for_beacons, all_beacon_combinations)
    return flatten_generator(all_intersection_points)

def is_limit_intersection_on_target(beacons: List[Beacon], limit_intersection: IntersectionPoint) -> bool:
    beacons_to_check = filter(lambda b: b != limit_intersection.limit1.beacon and b != limit_intersection.limit2.beacon, beacons)
    return all(b.is_within_limits(limit_intersection.point) for b in beacons_to_check)
    
def filter_limit_intersections(beacons: List[Beacon], limit_intersections: List[IntersectionPoint]) -> List[IntersectionPoint]:
    return filter(lambda ip: is_limit_intersection_on_target(beacons, ip), limit_intersections)


def calculate_furthest_point(x, points):
    c = Centroid(Point(*x), 0)
    c.set_radius(points)
    return c.radius


def do_multilat(beacons, tag='', start_point=None):
    closest_beacon = min(beacons, key=lambda b: b.limits[1])
    if start_point is None:
        furthest_beacon_point = max(beacons, key=lambda b: b.limits[1]).point
        start_point = furthest_beacon_point
    #get a point in the intersection

    result = minimize(
        overlap_error,
        start_point.to_tuple(),
        args=(beacons),
        method='L-BFGS-B',           # The optimisation algorithm
        #jac=grad_func,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })

    #print("initial minimize\n", result)
    intersection_point = Point(*result.x)

    if not all([b.is_within_limits(intersection_point) for b in beacons]):
        plotBeacons(beacons, actual=intersection_point, options={'tag': str(tag)})
        raise IntersectionError(f'the first minimization did not settle within the intersection', result)

    #now try to map the boundaries of the intersection
    bounds = map_bounds_of_intersection(beacons, intersection_point, tag=tag)

    plotBeacons(beacons, actual=intersection_point, preds=bounds, options={'tag': f'debug-{tag}'})

    #now calculate the smallest circle which encloses the boundary points we mapped
    #this is an optimization based approach.
    #there are more efficient mathematical algorithms for doing this.
    result = minimize(calculate_furthest_point, intersection_point.to_tuple(), bounds, method='L-BFGS-B', options={'ftol': 1e-5, 'maxiter': 1e6})
    return Centroid(Point(*result.x), result.fun), bounds

def do_multilat2(beacons: List[Beacon], tag=''):
    all_intersections = compute_limit_intersections(beacons)
    boundary_intersections = filter_limit_intersections(beacons, all_intersections)
    bounds = [b.point for b in boundary_intersections]
    bounds_centroid = Centroid.new_empty()
    bounds_centroid.add_points(bounds)
    bounds_centroid.divide()

    #now calculate the smallest circle which encloses the boundary points we mapped
    #this is an optimization based approach.
    #there are more efficient mathematical algorithms for doing this.
    result = minimize(calculate_furthest_point, bounds_centroid.point.to_tuple(), bounds, method='L-BFGS-B', options={'ftol': 1e-5, 'maxiter': 1e6})
    return Centroid(Point(*result.x), result.fun), bounds

if __name__ == '__main__':

    buckets = [[0,500], [500,1000], [1000, 2000], [2000, 5000], [5000, 10000], [10000, 20000]]
    center = Point(45,45)
    beacon_points = generate_triangle_points(center, 1100)
    for n,actual in enumerate([center, get_point_at_distance_and_bearing(center, 1300, 300)]): 
        beacons = simulate_service_distances(actual, beacon_points, buckets)
        try:
            centroid, bounds = do_multilat2(beacons, tag=str(n))
            
            plotBeacons(beacons, actual=actual, preds=bounds, centroid=centroid, options={'tag': str(n)})
        except IntersectionError as e:
            print(e)
            pass
    beacon_points = generate_triangle_points(center, 1300) # + [get_point_at_distance_and_bearing(center, 15900, 180)]
    for n,actual in enumerate([Point(45.02245801238342, 45.03176042407004) ]):
        beacons = simulate_service_distances(actual, beacon_points, buckets)
        
        centroid, bounds = do_multilat2(beacons, tag=str(n+4))
        
        plotBeacons(beacons, actual=actual, preds=bounds, centroid=centroid, options={'tag': str(n+4)})