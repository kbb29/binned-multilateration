import unittest, math

from sklearn.ensemble import GradientBoostingRegressor
from multilat import Point, Beacon, generate_triangle_points, grad_func, great_circle_distance, overlap_error, overlap_error_find_bounds, error_component_find_bounds_optimize, error_component_find_bounds_other, error_component_regular, grad_func_find_bounds, grad_component_find_bounds_optimize, grad_component_find_bounds_other, grad_component_regular, great_circle_distance, get_point_at_distance_and_bearing, normalize_gradient_components

class TestGradMethods(unittest.TestCase):

    #overlap_error(x, beacons, beacons_maximize, beacons_minimize)

    def generate_test_beacons(self):
        self.center = Point(45.0, 45.0)
        self.radius = 750
        self.points = generate_triangle_points(self.center, self.radius)
        self.beacons = [Beacon(p, [500, 1000]) for p in self.points]
        self.beacon_separation = 1297.5938693081803

    def test_overlap_func(self):
        self.generate_test_beacons()
        
        for beacon in self.beacons:
            minimum = (beacon.limits[1] - beacon.limits[0]) / -4
            self.assertAlmostEqual(error_component_regular(beacon, self.center), minimum, delta=minimum/-200)
        
        self.assertAlmostEqual(overlap_error(self.center.to_tuple(), self.beacons), -125, delta=150/200)

    def test_overlap_at_point0(self):
        self.generate_test_beacons()

        self.assertEqual(great_circle_distance(self.beacons[0].point, self.points[0]), 0.0)
        beacon_separation = self.beacon_separation
        self.assertAlmostEqual(great_circle_distance(self.beacons[1].point, self.points[0]), beacon_separation)
        self.assertAlmostEqual(great_circle_distance(self.beacons[2].point, self.points[0]), beacon_separation)
        self.assertAlmostEqual(great_circle_distance(self.points[2], self.points[1]), beacon_separation)
        self.assertEqual(error_component_regular(self.beacons[0], self.points[0]), 500.0)
        self.assertAlmostEqual(error_component_regular(self.beacons[1], self.points[0]), beacon_separation - 1000)
        self.assertAlmostEqual(error_component_regular(self.beacons[2], self.points[0]), beacon_separation - 1000)
    
        self.assertAlmostEqual(overlap_error(self.points[0].to_tuple(), self.beacons), (500.0 + (beacon_separation - 1000)*2) / len(self.beacons) )

    def do_error_component_find_bounds_maximize_test(self, beacon, distance, bearing):
        point = get_point_at_distance_and_bearing(beacon.point, distance, bearing)
        dist = great_circle_distance(point, beacon.point)
        if distance > beacon.limits[1]:
            self.assertAlmostEqual(error_component_find_bounds_optimize(beacon, point, True), (dist - beacon.limits[1])**2)
        else:
            self.assertAlmostEqual(error_component_find_bounds_optimize(beacon, point, True), (beacon.limits[1] - dist)**1.2)

    def do_error_component_find_bounds_minimize_test(self, beacon, distance, bearing):
        point = get_point_at_distance_and_bearing(beacon.point, distance, bearing)
        dist = great_circle_distance(point, beacon.point)
        if dist > beacon.limits[0]:
            self.assertAlmostEqual(error_component_find_bounds_optimize(beacon, point, False), (dist - beacon.limits[0])**1.2)
        else:
            self.assertAlmostEqual(error_component_find_bounds_optimize(beacon, point, False), (beacon.limits[0] - dist)**2)

    def test_error_component_maximize(self):
        self.generate_test_beacons()

        b = self.beacons[0]
        self.do_error_component_find_bounds_maximize_test(b, 750, 0)
        self.do_error_component_find_bounds_maximize_test(b, 750, 90)
        self.do_error_component_find_bounds_maximize_test(b, 750, 305)
        self.do_error_component_find_bounds_maximize_test(b, 1250, 0)
        self.do_error_component_find_bounds_maximize_test(b, 1250, 90)
        self.do_error_component_find_bounds_maximize_test(b, 1250, 305)

    def test_error_component_minimize(self):
        self.generate_test_beacons()

        b = self.beacons[0]
        self.do_error_component_find_bounds_minimize_test(b, 750, 0)
        self.do_error_component_find_bounds_minimize_test(b, 750, 90)
        self.do_error_component_find_bounds_minimize_test(b, 750, 305)
        self.do_error_component_find_bounds_minimize_test(b, 1250, 0)
        self.do_error_component_find_bounds_minimize_test(b, 1250, 90)
        self.do_error_component_find_bounds_minimize_test(b, 1250, 305)
        

    def test_overlap_maximize(self):
        self.generate_test_beacons()
        self.assertAlmostEqual(overlap_error_find_bounds(self.center.to_tuple(), self.beacons, 0, True), 250**1.2/3, delta=1.1)
        b = self.beacons[1]
        expected_error = ((self.beacon_separation - self.beacons[0].limits[1])**2 + self.beacons[1].limits[0]**2 + (self.beacon_separation - self.beacons[2].limits[1])**2) / 3
        self.assertAlmostEqual(overlap_error_find_bounds(self.beacons[1].point.to_tuple(), self.beacons, 0, True), expected_error, places=4)
        expected_error = ((self.beacon_separation - self.beacons[0].limits[1])**2 + self.beacons[1].limits[1]**1.2 + (self.beacon_separation - self.beacons[2].limits[1])**2) / 3
        self.assertAlmostEqual(overlap_error_find_bounds(self.beacons[1].point.to_tuple(), self.beacons, 1, True), expected_error, places=4)

    def test_overlap_minimize(self):
        self.generate_test_beacons()
        self.assertAlmostEqual(overlap_error_find_bounds(self.center.to_tuple(), self.beacons, 0, False), 250**1.2/3, delta=1.1)
        b = self.beacons[1]
        expected_error = ((self.beacon_separation - self.beacons[0].limits[0])**1.2 + self.beacons[1].limits[0]**2 + (self.beacon_separation - self.beacons[2].limits[1])**2) / 3
        self.assertAlmostEqual(overlap_error_find_bounds(self.beacons[1].point.to_tuple(), self.beacons, 0, False), expected_error, places=4)
        expected_error = ((self.beacon_separation - self.beacons[0].limits[1])**2 + self.beacons[1].limits[0]**2 + (self.beacon_separation - self.beacons[2].limits[1])**2) / 3
        self.assertAlmostEqual(overlap_error_find_bounds(self.beacons[1].point.to_tuple(), self.beacons, 1, False), expected_error, places=4)

    def assertTupleAlmostEqual(self, a, b, places=None, delta=None):
        self.assertEqual(len(a), len(b))

        opts = {}
        if delta is not None:
            opts['delta'] = delta
            if places is not None:
                raise TypeError('specify places or delta, not both')
        elif places is not None:
            opts['places'] = places

        
        for i in range(len(a)):
            self.assertAlmostEqual(a[i], b[i], **opts)



    def do_gradient_test_at_point(self, test_point):
        #generate points 10m north south east and west of this point
        test_point_n = get_point_at_distance_and_bearing(test_point, 1, 0) 
        test_point_s = get_point_at_distance_and_bearing(test_point, 1, 180)
        test_point_e = get_point_at_distance_and_bearing(test_point, 1, 90)
        test_point_w = get_point_at_distance_and_bearing(test_point, 1, 270)
        
        val_n = overlap_error(test_point_n.to_tuple(), self.beacons)
        val_s = overlap_error(test_point_s.to_tuple(), self.beacons)
        val_e = overlap_error(test_point_e.to_tuple(), self.beacons)
        val_w = overlap_error(test_point_w.to_tuple(), self.beacons)

        #print([b.is_within_limits(test_point) for b in self.beacons])
        grad_approx = [(val_n-val_s)/2, (val_e-val_w)/2]
        grad = grad_func(test_point.to_tuple(), self.beacons)   

        #print(val_n, val_s, val_e, val_w)
        #print(grad, grad_approx)
        self.assertTupleAlmostEqual(grad, grad_approx, delta=0.01)     

    def test_grad(self):
        self.generate_test_beacons()

        self.do_gradient_test_at_point(self.center)

        test_point = get_point_at_distance_and_bearing(self.center, 500, 0.0) #this point is 250m south of the top point of the triangle
        self.do_gradient_test_at_point(test_point)
       
        test_point = get_point_at_distance_and_bearing(self.center, 200, 0.0) #this point is 550m south of the top point of the triangle
        self.do_gradient_test_at_point(test_point)

    def do_gradient_test_find_bounds_at_point(self, test_point, idx_to_optimize, maximize):
        #generate points 10m north south east and west of this point
        test_point_n = get_point_at_distance_and_bearing(test_point, 1, 0) 
        test_point_s = get_point_at_distance_and_bearing(test_point, 1, 180)
        test_point_e = get_point_at_distance_and_bearing(test_point, 1, 90)
        test_point_w = get_point_at_distance_and_bearing(test_point, 1, 270)
        
        val_n = overlap_error_find_bounds(test_point_n.to_tuple(), self.beacons, idx_to_optimize, maximize)
        val_s = overlap_error_find_bounds(test_point_s.to_tuple(), self.beacons, idx_to_optimize, maximize)
        val_e = overlap_error_find_bounds(test_point_e.to_tuple(), self.beacons, idx_to_optimize, maximize)
        val_w = overlap_error_find_bounds(test_point_w.to_tuple(), self.beacons, idx_to_optimize, maximize)
        
        #print([b.is_within_limits(test_point) for b in self.beacons])
        grad_approx = [(val_n-val_s)/2, (val_e-val_w)/2]
        grad = grad_func_find_bounds(test_point.to_tuple(), self.beacons, idx_to_optimize, maximize)   

        #print(val_n, val_s, val_e, val_w)
        #print(grad, grad_approx)
        self.assertTupleAlmostEqual(grad, grad_approx, places=0)     

    def test_grad_find_bounds_maximize(self):
        self.generate_test_beacons()
        self.do_gradient_test_find_bounds_at_point(self.center, 0, True)
        self.do_gradient_test_find_bounds_at_point(self.center, 1, True)
        self.do_gradient_test_find_bounds_at_point(self.center, 2, True)
        
        test_point = get_point_at_distance_and_bearing(self.center, 500, 0.0) #this point is 250m south of the top point of the triangle
        self.do_gradient_test_find_bounds_at_point(test_point, 0, True)
        self.do_gradient_test_find_bounds_at_point(test_point, 1, True)
        self.do_gradient_test_find_bounds_at_point(test_point, 2, True)
        
        test_point = get_point_at_distance_and_bearing(self.center, 200, 0.0) #this point is 550m south of the top point of the triangle
        self.do_gradient_test_find_bounds_at_point(test_point, 0, True)
        self.do_gradient_test_find_bounds_at_point(test_point, 1, True)
        self.do_gradient_test_find_bounds_at_point(test_point, 2, True)

    def test_grad_find_bounds_minimize(self):
        self.generate_test_beacons()
        
        self.do_gradient_test_find_bounds_at_point(self.center, 0, False)
        self.do_gradient_test_find_bounds_at_point(self.center, 1, False)
        self.do_gradient_test_find_bounds_at_point(self.center, 2, False)
        
        test_point = get_point_at_distance_and_bearing(self.center, 500, 0.0) #this point is 250m south of the top point of the triangle
        self.do_gradient_test_find_bounds_at_point(test_point, 0, False)
        self.do_gradient_test_find_bounds_at_point(test_point, 1, False)
        self.do_gradient_test_find_bounds_at_point(test_point, 2, False)
       
        test_point = get_point_at_distance_and_bearing(self.center, 200, 0.0) #this point is 550m south of the top point of the triangle
        self.do_gradient_test_find_bounds_at_point(test_point, 0, False)
        self.do_gradient_test_find_bounds_at_point(test_point, 1, False)
        self.do_gradient_test_find_bounds_at_point(test_point, 2, False)



    def test_grad_norm(self):
        self.assertTupleAlmostEqual(normalize_gradient_components(4,5,norm_to=1), (0.6246950475544243, 0.7808688094430304))
        self.assertTupleAlmostEqual(normalize_gradient_components(4,5,norm_to=0.2), (0.12493900951088485, 0.15617376188860607))
        self.assertTupleAlmostEqual(normalize_gradient_components(4,5,norm_to=0.5), (0.31234752377721214, 0.3904344047215152))
        
        self.assertTupleAlmostEqual(normalize_gradient_components(4,0,norm_to=0.5), (0.5, 0))
        self.assertTupleAlmostEqual(normalize_gradient_components(0,3,norm_to=0.2), (0, 0.2))

if __name__ == '__main__':
    unittest.main()