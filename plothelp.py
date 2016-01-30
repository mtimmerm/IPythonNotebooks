import matplotlib
import matplotlib.pyplot as plt
import numpy
import math

def hermite_interp(dar,var,i):
    n = math.floor(i)
    if n>=len(dar)-1:
        n=len(dar)-2
    t = i-n

    # Hermite basis functions
    h00 = (2.0 * t**3) - (3.0 * t**2) + 1.0
    h10 = t**3.0  - (2.0 * t**2) + t
    h01 = (-2.0* t**3) + (3.0 * t**2)
    h11 = t**3 - t**2
    # Compute the interpolated value of "y"
    return h00*dar[n] + h10*var[n] + h01*dar[n+1] + h11*var[n+1]

# first derivative of hermite spline
def hermite_interp1(dar,var,i):
    n = math.floor(i)
    if n>=len(dar)-1:
        n=len(dar)-2
    t = i-n

    h00 = (6.0 * t**2) - (6.0 * t)
    h10 = 3.0*t**2  - (4.0 * t) + 1
    h01 = (-6.0* t**2) + (6.0 * t)
    h11 = 3.0*t**2 - 2.0*t
    return h00*dar[n] + h10*var[n] + h01*dar[n+1] + h11*var[n+1]

#second derivative of hermite spline
def hermite_interp2(dar,var,i):
    n = math.floor(i)
    if n>=len(dar)-1:
        n=len(dar)-2
    t = i-n

    # Hermite basis functions
    h00 = (12.0 * t) - 6.0
    h10 = 6.0*t  - 4.0
    h01 = (-12.0* t) + 6.0
    h11 = 6.0*t - 2.0
    # Compute the interpolated value of "y"
    return h00*dar[n] + h10*var[n] + h01*dar[n+1] + h11*var[n+1]
