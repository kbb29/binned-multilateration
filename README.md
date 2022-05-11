# binned-multilateration

![](images/trilat.png)

This is a node module for triangulating/trilaterating/multilaterating an object when you know how far it is from two or more beacons of known location *and* the distances are obfuscated into buckets eg. 0m-500m, 500m-1000m etc.

It could also be called annulus-trilateration, annuli-trilateration or torus-trilateration, but I've gone for binned trilateration because torus is incorrect and annulus makes certain people snigger...

This is a python version of [bucket-trilateration](https://github.com/kbb29/bucket-trilateration).  However, the JS and Python algorithms are very different.  This one (python) has complexity greater than n**2, so will not perform well with large numbers of beacons.  But this one is more likely to get the correct answer.


## usage

in this directory:

```pip3 install .```


then you can do something like this (or run demo.py):

```
from binned_multilateration import simulate_service_distances, get_point_at_distance_and_bearing, do_multilat, generate_triangle_points, Point, plotBeacons

#simulates a service by returning the binned distances for an actual point on the earth
def dummy_service(beacon_points):
  # these are the bins into which the service approximates the distance
  buckets = [[0,500], [500,1000], [1000, 2000], [2000, 5000], [5000, 10000], [10000, 20000]]
  actual = get_point_at_distance_and_bearing(Point(45, 45), 1300, 300)
  beacons = simulate_service_distances(actual, beacon_points, buckets)
  return beacons


beacon_points = generate_triangle_points(Point(45,45), 1100)

#in reality the beacon distance values would be returned by some service which we do not control
#for the purpose of this example, we will use the dummy_service function above
beacons = dummy_service(beacon_points)
centroid, bounds = do_multilat(beacons)

#opens a viz in browser
plotBeacons(beacons, preds=bounds, centroid=centroid)


```


## algorithm details

The algorithm is based around the concept of beacons.  A beacons is a fixed point on earth.  You know that the distance from the beacon to some target object is between two limits eg. between 500m and 1000m away.  You have this kind of distance estimate for multiple beacons and you want to know the area in which the target object could be located.

1. first we calculate the positions of all the points where two beacon limits intersect.
1. Then we calculate which of these limit intersections lie on the boundary of the target area.
1. Then we find the smallest circle which encloses all the target boundary points and return this.


