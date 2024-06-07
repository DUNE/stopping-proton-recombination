import numpy as np
from scipy import interpolate


def f(x):
    x_points = [-90.0,-70.0,-50.0,-30.0,-10.0]
    y_points = [0.517,0.532,0.548,0.564,0.587]
    tck = interpolate.splrep(x_points, y_points)
    return interpolate.splev(x, tck)


# for x, these are the mean, min and max values
mu = -43.3
min_x, max_x = -71.8, -16.7
minus_one, plus_one = -49.1, -37.5


print(f"Average E field for tracks: {f(mu)} kV/cm, range of values [{f(min_x)},{f(max_x)}]")
print(f"Average E field for tracks: {f(mu)} kV/cm, one sigma range of values [{f(minus_one)},{f(plus_one)}]")
