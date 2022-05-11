import unittest, math
from binned_multilateration import limit_intersection, Beacon, Limit, Point, IntersectionPoint, get_point_at_distance_and_bearing

class TestBearingMethods(unittest.TestCase):
    def test_given(self):
        b1 = Beacon(Point(37.673442, -90.234036), (199090.0, 199090 + 90000))
        b2 = Beacon(Point(36.109997, -90.953669), (268540, 268540 + 90000))
        l1 = Limit(b1, True)
        l2 = Limit(b2, True)
        intersection_points = limit_intersection(l1, l2)
        self.assertEqual(intersection_points[0], IntersectionPoint(Point(36.98931105153329, -88.15142628069124), l1, l2))
        self.assertEqual(intersection_points[1], IntersectionPoint(Point(38.238379609457695, -92.39048549120311), l1, l2))

    def test_smaller(self):
        b1 = Beacon(Point(37.673442, -90.234036), (1000.0, 2000))
        b2 = Beacon(get_point_at_distance_and_bearing(b1.point, 2000.0, 45.0), (2000.0, 5000))
        l1 = Limit(b1, True)
        l2 = Limit(b2, True)
        intersection_points = limit_intersection(l1, l2)
        self.assertEqual(intersection_points[0], IntersectionPoint(Point(37.68117930503295, -90.2398426645948), l1, l2))
        self.assertEqual(intersection_points[1], IntersectionPoint(Point(37.668846060772346, -90.22426101493188), l1, l2))
        
        
if __name__ == '__main__':
    unittest.main()
