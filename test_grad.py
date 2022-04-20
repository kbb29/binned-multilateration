import unittest, math
from multilat import Point, Beacon, generate_triangle_points, great_circle_distance, overlap_error, error_component_regular, great_circle_distance, get_point_at_distance_and_bearing

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
    

        self.assertAlmostEqual(overlap_error([self.points[0].lat, self.points[0].lon], self.beacons, [], []), (500.0 + (beacon_separation - 1000)*2) / 3 )

    def test_grad(self):
        self.assertEqual(True, False)

        
if __name__ == '__main__':
    unittest.main()