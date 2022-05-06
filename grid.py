
from multilat import IntersectionError, Point, Beacon, BoundSquare, grad_func, generate_triangle_points, get_bearing, get_point_at_distance_and_bearing, do_multilat, great_circle_distance, simulate_service_distances, plotBeacons
from plot import plotBeacons
import math

def gridgen(center, dimension, num_points):
    bs = BoundSquare(center, dimension)
    n = math.floor(math.sqrt(num_points))
    latinc = (bs.north - bs.south) / (n - 1)
    loninc = (bs.east - bs.west) / (n - 1)
    for lati in range(n):
        lat = latinc * lati + bs.south
        for loni in range(n):
            lon = loninc * loni + bs.west
            yield Point(lat, lon)

if __name__ == '__main__':

    buckets = [[0,500], [500,1000], [1000, 2000], [2000, 5000], [5000, 10000], [10000, 20000]]
    center = Point(45,45)
    beacon_points = generate_triangle_points(center, 1300)
    num_points = 10**2
    centroids = []
    num_not_settled = 0
    num_settled_2nd_attempt = 0
    num_not_in_centroid = 0
    for i,pt in enumerate(gridgen(center, 5000, num_points)):
        beacons = simulate_service_distances(pt, beacon_points, buckets)
        furthest_beacon_point = max(beacons, key=lambda b: b.limits[1]).point
        try:
            centroid, bounds = do_multilat(beacons, start_point=furthest_beacon_point)
        except IntersectionError as e:
            next_attempt = get_point_at_distance_and_bearing(center, great_circle_distance(center, furthest_beacon_point), get_bearing(furthest_beacon_point, center))
            try:
                centroid, bounds = do_multilat(beacons, start_point=next_attempt)
                num_settled_2nd_attempt += 1
            except IntersectionError as e:
                print(f'gradient at unsettled {grad_func(e.result.x, beacons)}')
                print(e.result)
                plotBeacons(beacons, preds=[Point(*e.result.x)], actual=pt, options={'tag': str(i)})
                num_not_settled += 1
            else:
                centroids.append(centroid)
        else:
            centroids.append(centroid)

        if not centroid.includes(pt):
            print(f'{i} ({pt.to_tuple()}): point not in centroid!')

            num_not_in_centroid += 1
            plotBeacons(beacons, actual=pt, preds=bounds, centroid=centroid)

    radii = [c.radius for c in centroids]
    avg_radius = sum(radii) / len(radii)
    max_radius = max(radii)

    print(len(radii))
    print('not settled', num_not_settled)
    print('settled 2nd attempt', num_settled_2nd_attempt)
    print('not in centroid',  num_not_in_centroid)
    print(avg_radius, max_radius, min(radii))

