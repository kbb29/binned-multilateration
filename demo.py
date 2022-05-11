from binned_multilateration import Point, simulate_service_distances, do_multilat, get_point_at_distance_and_bearing, generate_triangle_points, plotBeacons


if __name__ == '__main__':

    buckets = [[0,500], [500,1000], [1000, 2000], [2000, 5000], [5000, 10000], [10000, 20000]]
    center = Point(45,45)
    beacon_points = generate_triangle_points(center, 1100)
    for n,actual in enumerate([center, get_point_at_distance_and_bearing(center, 1300, 300)]): 
        beacons = simulate_service_distances(actual, beacon_points, buckets)
        centroid, bounds = do_multilat(beacons, tag=str(n))
            
        plotBeacons(beacons, actual=actual, preds=bounds, centroid=centroid, options={'tag': str(n)})
        
    beacon_points = generate_triangle_points(center, 1300) # + [get_point_at_distance_and_bearing(center, 15900, 180)]
    for n,actual in enumerate([Point(45.02245801238342, 45.03176042407004) ]):
        beacons = simulate_service_distances(actual, beacon_points, buckets)
        
        centroid, bounds = do_multilat(beacons, tag=str(n+4))
        
        plotBeacons(beacons, actual=actual, preds=bounds, centroid=centroid, options={'tag': str(n+4)})
