import math
from scipy.optimize import minimize
import pathlib, chevron, os
from typing import Tuple, List

template = (pathlib.Path(__file__).parent / 'map.mustache').open('r').read()



def generateDonut(lat, lon, outerradius, innerradius, options={}) :
    return f'''L.donut([{lat}, {lon}], {{
        color: '{options.get('color', 'red')}',
        fillColor: '{options.get('color', 'red')}',
        fillOpacity: 0.4,
        radius: {outerradius},
        innerRadius: {innerradius}
    }})'''

def generateCircle(lat, lon, radius, options={}):
    return f'''L.circle([{lat}, {lon}], {{
        color: '{options.get('color', "red")}',
        fillColor: '{options.get('color', "red")}',
        fillOpacity: {options.get('opacity', "0.0")},
        radius: {radius}
    }})'''

def generateMarker(lat, lon, title):
    return f"L.marker([{lat}, {lon}], {{title: '{title}'}})"

def generateAddToFeaturesAndMap(featureCode):
    return f"features.push({featureCode}.addTo(map));\n"

def plotBeacons(beacons, preds=[], centroid=None, actual=None, options={}):
    featureJS = ""
    for beacon in beacons:
        featureJS += generateAddToFeaturesAndMap(generateDonut(beacon.point.lat, beacon.point.lon, beacon.limits[1], beacon.limits[0], {'color': "red"}))
        featureJS += generateAddToFeaturesAndMap(generateMarker(beacon.point.lat, beacon.point.lon, "kevin"))

    if actual:
        featureJS += generateAddToFeaturesAndMap(generateMarker(actual.lat, actual.lon, 'actual'))
    
    for pred in preds:
        featureJS += generateAddToFeaturesAndMap(generateMarker(pred.lat, pred.lon, 'pred'))

    if centroid:
        featureJS += generateAddToFeaturesAndMap(generateCircle(centroid.point.lat, centroid.point.lon, centroid.radius, {'color': 'blue', 'opacity': "0.1"}))

    output = chevron.render(template, {'features': featureJS});
    
    tag = options.get('tag', 'test')

    fn = pathlib.Path('/tmp') / f'map-{tag}.html'
    print('fn', fn)
    with fn.open('w') as fh:
        fh.write(output)

    print('fn2', fn)
    os.system(f'xdg-open {fn}')

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


def error_component_regular(beacon, pt):
        distance_calculated = great_circle_distance(pt, beacon.point)
        err = 0.0
        if distance_calculated < beacon.limits[0]:
            err = beacon.limits[0] - distance_calculated
        elif distance_calculated > beacon.limits[1]:
            err = distance_calculated - beacon.limits[1]
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
def overlap_error(x, beacons, beacons_maximize, beacons_minimize) -> float:
    mse = 0.0
    pt = Point(*x)
    errs = [error_component_regular(beacon, pt) for beacon in beacons]
    errs += [error_component_maximize(beacon, pt) for beacon in beacons_maximize]
    errs += [error_component_minimize(beacon, pt) for beacon in beacons_minimize]
    mse = sum(errs)
    return mse / len(beacons)

def grad_component_regular(beacon, pt) -> Tuple[float,float]:
        dist = great_circle_distance(pt,beacon.point)
        if dist < beacon.limits[0]:
            #then gradient is magnitude 1 in the direction from the point to the beacon center
            lat_component = beacon.point.lat - pt.lat
            lon_component = (beacon.point.lon - pt.lon) * math.cos(math.radians(pt.lat))
        elif dist > beacon.limits[1]:
            #then gradient is magnitude 1 in the direciton from the beacon center to the point
            lat_component = pt.lat - beacon.point.lat
            lon_component = (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))
        else:
            lat_component = 0.0
            lon_component = 0.0
        #if we are inside the bounds, then this contributes 0 to the grad
        return lat_component, lon_component

def grad_component_maximize(beacon, pt) -> Tuple[float,float]:
        dist = great_circle_distance(pt,beacon.point)
        if dist < beacon.limits[1]:
            #then gradient is magnitude 0.5 in the direction from the point to the beacon center
            lat_component = 0.5 * (beacon.point.lat - pt.lat)
            lon_component = 0.5 * (beacon.point.lon - pt.lon) * math.cos(math.radians(pt.lat))
        elif dist > beacon.limits[1]:
            #then gradient is magnitude 1 in the direction from the beacon center to the point
            lat_component = pt.lat - beacon.point.lat
            lon_component = (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))
        else:
            lat_component = 0.0
            lon_component = 0.0
        #if we are inside the bounds, then this contributes 0 to the grad
        return lat_component, lon_component

def grad_component_minimize(beacon, pt) -> Tuple[float,float]:
        dist = great_circle_distance(pt,beacon.point)
        if dist < beacon.limits[0]:
            #then gradient is magnitude 1 in the direction from the point to the beacon center
            lat_component = beacon.point.lat - pt.lat
            lon_component = (beacon.point.lon - pt.lon) * math.cos(math.radians(pt.lat))
        elif dist > beacon.limits[0]:
            #then gradient is magnitude 0.5 in the direction from the beacon center to the point
            lat_component = 0.5 * (pt.lat - beacon.point.lat)
            lon_component = 0.5 * (pt.lon - beacon.point.lon) * math.cos(math.radians(pt.lat))
        else:
            lat_component = 0.0
            lon_component = 0.0
        #if we are inside the bounds, then this contributes 0 to the grad
        return lat_component, lon_component

def grad_func(x, beacons, maximize_beacons, minimize_beacons) -> List[float]:
    grad = [0.0,0.0]
    pt = Point(*x)
    component_tuples = [grad_component_regular(beacon, pt) for beacon in beacons]
    component_tuples += [grad_component_maximize(beacon, pt) for beacon in maximize_beacons]
    component_tuples += [grad_component_minimize(beacon, pt) for beacon in minimize_beacons]

    for lat_component, lon_component in component_tuples:    
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

def find_outer_bound_for_beacon(beacons, start_point, optimize_index=None, maximize=True):
    optimize_beacon = beacons[optimize_index]
    maximize_beacons = [optimize_beacon] if maximize else []
    minimize_beacons = [optimize_beacon] if not maximize else []
    standard_beacons = beacons[:optimize_index] + beacons[optimize_index+1:]
    result = minimize(
        overlap_error,                         # The error function
        [start_point.lat, start_point.lon],            # The initial guess
        args=(standard_beacons, maximize_beacons, minimize_beacons), # Additional parameters for mse
        method='L-BFGS-B',           # The optimisation algorithm
        #jac=grad_func,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })

    print(result)
    pt = Point(*result.x)
    return pt

def calculate_furthest_point(x, points):
    c = Centroid(Point(*x), 0)
    c.set_radius(points)
    return c.radius


def do_multilat(beacons):
    furthest_beacon_point = max(beacons, key=lambda b: b.limits[1]).point
    #get a point in the intersection

    result = minimize(
        overlap_error,
        [furthest_beacon_point.lat, furthest_beacon_point.lon],
        args=(beacons, [], []),
        method='L-BFGS-B',           # The optimisation algorithm
        jac=grad_func,
        options={
            'ftol':1e-5,         # Tolerance
            'maxiter': 1e+7      # Maximum iterations
        })

    print(result)
    intersection_point = Point(*result.x)

    bounds  = [find_outer_bound_for_beacon(beacons, intersection_point, optimize_index=i, maximize=True ) for i in range(len(beacons))]
    bounds += [find_outer_bound_for_beacon(beacons, intersection_point, optimize_index=i, maximize=False) for i in range(len(beacons))]
    centroid = Centroid(Point(0.0, 0.0), 0.0)
    [centroid.add_point(p) for p in bounds]
    centroid.divide()

    result = minimize(calculate_furthest_point, [centroid.point.lat ,centroid.point.lon], bounds, method='L-BFGS-B', options={'ftol': 1e-5, 'maxiter': 1e6})
    return Centroid(Point(*result.x), result.fun), bounds

for n,actual in enumerate([center, get_point_at_distance_and_bearing(center, 1300, 300)]):
    beacons = simulate_service_distances(actual, beacon_points, buckets)
    centroid, bounds = do_multilat(beacons)

    plotBeacons(beacons, actual=actual, preds=bounds, centroid=centroid, options={'tag': str(n)})
