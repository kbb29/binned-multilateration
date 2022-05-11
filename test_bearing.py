import unittest, math
from binned_multilateration import Point, Beacon, great_circle_distance, get_point_at_distance_and_bearing, get_bearing, generate_triangle_points

class TestBearingMethods(unittest.TestCase):

    #overlap_error(x, beacons, beacons_maximize, beacons_minimize)

    def generate_test_beacons(self):
        self.center = Point(45.0, 45.0)
        self.radius = 750
        self.points = generate_triangle_points(self.center, self.radius)
        self.beacons = [Beacon(p, [500, 1000]) for p in self.points]

    def test_get_bearing(self):
        self.assertAlmostEqual(0, get_bearing(Point(45,45), Point(46,45)), delta=1)
        self.assertAlmostEqual(90, get_bearing(Point(45,45), Point(45,46)), delta=1)
        self.assertAlmostEqual(180, get_bearing(Point(45,45), Point(44,45)), delta=1)
        self.assertAlmostEqual(270, get_bearing(Point(45,45), Point(45,44)), delta=1)

    def do_bearing_test(self, origin, distance, bearing):
        dest = get_point_at_distance_and_bearing(origin, distance, bearing)
        self.assertAlmostEqual(get_bearing(origin, dest), bearing, delta = bearing/300)
        self.assertAlmostEqual(great_circle_distance(origin, dest), distance, delta=distance/500 )

    def test_get_point(self):
        self.do_bearing_test(Point(45,45), 1000, 0)
        self.do_bearing_test(Point(45,45), 1000, 90)
        self.do_bearing_test(Point(45,45), 1000, 270)
        self.do_bearing_test(Point(45,45), 1000, 180)
        self.do_bearing_test(Point(45,45), 1000, 36.45)

        
if __name__ == '__main__':
    unittest.main()