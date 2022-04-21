import unittest, math
from multilat import Point, Beacon, generate_triangle_points, grad_func, great_circle_distance, overlap_error, error_component_regular, error_component_maximize, error_component_minimize, great_circle_distance, get_point_at_distance_and_bearing

class TestGradMethods(unittest.TestCase):

    #overlap_error(x, beacons, beacons_maximize, beacons_minimize)

    def generate_test_beacons(self):
        self.center = Point(45.0, 45.0)
        self.radius = 750
        self.points = generate_triangle_points(self.center, self.radius)
        self.beacons = [Beacon(p, [500, 1000]) for p in self.points]

    def test_overlap_func(self):
        self.generate_test_beacons()
        
        for beacon in self.beacons:
            self.assertEqual(error_component_regular(beacon, self.center), 0.0)
        
        self.assertEqual(overlap_error([45.0, 45.0], self.beacons, [], []), 0.0)

    def test_overlap_at_point0(self):
        self.generate_test_beacons()

        self.assertEqual(great_circle_distance(self.beacons[0].point, self.points[0]), 0.0)
        beacon_separation = 1297.5938693081803
        self.assertAlmostEqual(great_circle_distance(self.beacons[1].point, self.points[0]), beacon_separation)
        self.assertAlmostEqual(great_circle_distance(self.beacons[2].point, self.points[0]), beacon_separation)
        self.assertAlmostEqual(great_circle_distance(self.points[2], self.points[1]), beacon_separation)
        self.assertEqual(error_component_regular(self.beacons[0], self.points[0]), 500.0)
        self.assertAlmostEqual(error_component_regular(self.beacons[1], self.points[0]), beacon_separation - 1000)
        self.assertAlmostEqual(error_component_regular(self.beacons[2], self.points[0]), beacon_separation - 1000)
    
        self.assertAlmostEqual(overlap_error([self.points[0].lat, self.points[0].lon], self.beacons, [], []), (500.0 + (beacon_separation - 1000)*2) )

    def test_error_component_maximize(self):
        self.generate_test_beacons()

        b = self.beacons[0]
        outer_limit = 1000
        inner_limit = 500
        self.assertAlmostEqual(round(error_component_maximize(b, get_point_at_distance_and_bearing(b.point, 750, 0))), 0.5* (outer_limit - 750))
        self.assertAlmostEqual(round(error_component_maximize(b, get_point_at_distance_and_bearing(b.point, 750, 90))), 0.5* (outer_limit - 750))
        self.assertAlmostEqual(round(error_component_maximize(b, get_point_at_distance_and_bearing(b.point, 750, 305))), 0.5* (outer_limit - 750))
        
        self.assertAlmostEqual(round(error_component_maximize(b, get_point_at_distance_and_bearing(b.point, 1250, 0))), (1249 - outer_limit))
        self.assertAlmostEqual(round(error_component_maximize(b, get_point_at_distance_and_bearing(b.point, 1250, 90))), (1249 - outer_limit))
        self.assertAlmostEqual(round(error_component_maximize(b, get_point_at_distance_and_bearing(b.point, 1250, 305))), (1249 - outer_limit))


    def test_error_component_minimize(self):
        self.generate_test_beacons()

        b = self.beacons[0]
        outer_limit = 1000
        inner_limit = 500
        self.assertAlmostEqual(round(error_component_minimize(b, get_point_at_distance_and_bearing(b.point, 750, 0))), 0.5 * (750 - inner_limit))
        self.assertAlmostEqual(round(error_component_minimize(b, get_point_at_distance_and_bearing(b.point, 750, 90))), 0.5 * (750 - inner_limit))
        self.assertAlmostEqual(round(error_component_minimize(b, get_point_at_distance_and_bearing(b.point, 750, 305))), 0.5 * (750 - inner_limit))
        

        self.assertAlmostEqual(round(error_component_minimize(b, get_point_at_distance_and_bearing(b.point, 250, 0))), (inner_limit - 250))
        self.assertAlmostEqual(round(error_component_minimize(b, get_point_at_distance_and_bearing(b.point, 250, 90))), (inner_limit - 250))
        self.assertAlmostEqual(round(error_component_minimize(b, get_point_at_distance_and_bearing(b.point, 250, 305))), (inner_limit - 250))
        

    def test_overlap_maximize(self):
        self.generate_test_beacons()
        self.assertEqual(int(overlap_error([45.0, 45.0], self.beacons[1:], [self.beacons[0]], [])), 125)

    def test_grad(self):
        self.generate_test_beacons()

        test_point = get_point_at_distance_and_bearing(self.center, 500, 0.0) #this point is 250m south of the top point of the triangle
        #generate points 10m north south east and west of this point
        test_point_n = get_point_at_distance_and_bearing(test_point, 10, 0) 
        test_point_s = get_point_at_distance_and_bearing(test_point, 10, 180)
        test_point_e = get_point_at_distance_and_bearing(test_point, 10, 90)
        test_point_w = get_point_at_distance_and_bearing(test_point, 10, 270)
        
        val_n = overlap_error(test_point_n.to_tuple(), self.beacons, [], [])
        val_s = overlap_error(test_point_s.to_tuple(), self.beacons, [], [])
        val_e = overlap_error(test_point_e.to_tuple(), self.beacons, [], [])
        val_w = overlap_error(test_point_w.to_tuple(), self.beacons, [], [])

        grad_approx = [(val_n-val_s)/20, (val_e-val_w)/20]
        grad = grad_func(test_point.to_tuple(), self.beacons, [], [])   

        #print(grad, grad_approx)
        self.assertAlmostEqual(grad[0], grad_approx[0], 2)
        self.assertAlmostEqual(grad[1], grad_approx[1], 2)     
if __name__ == '__main__':
    unittest.main()