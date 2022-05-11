
from binned_multilateration import Point, Beacon, BoundSquare, generate_triangle_points, get_point_at_distance_and_bearing, do_multilat, great_circle_distance, simulate_service_distances, plotBeacons
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
    num_not_in_centroid = 0
    for i,pt in enumerate(gridgen(center, 5000, num_points)):
        beacons = simulate_service_distances(pt, beacon_points, buckets)
        furthest_beacon_point = max(beacons, key=lambda b: b.limits[1]).point
        centroid, bounds = do_multilat(beacons)
        centroids.append(centroid)

        if not centroid.includes(pt):
            print(f'{i} ({pt.to_tuple()}): point not in centroid!')

            num_not_in_centroid += 1
            plotBeacons(beacons, actual=pt, preds=bounds, centroid=centroid)

    radii = [c.radius for c in centroids]
    avg_radius = sum(radii) / len(radii)
    max_radius = max(radii)

    print(len(radii))
    print('not in centroid',  num_not_in_centroid)
    print(avg_radius, max_radius, min(radii))

