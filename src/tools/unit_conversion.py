import numpy as np
import math


def ms2kt(cas_ms):

    cas_kt = 1.94384 * cas_ms

    return cas_kt


def kt2ms(cas_kt):

    cas_ms = cas_kt / 1.94384

    return cas_ms


def ft2m(ft):

    m = 0.3048 * ft

    return m


def m2ft(m):

    ft = m / 0.348

    return ft


def rad2deg(rad):

    deg = rad / 0.0174533

    return deg


def deg2rad(deg):

    rad = 0.0174533 * deg

    return rad


def deg2perc(deg):

    precentage = np.tan(np.radians(deg)) * 100

    return precentage


def perc2rad(perc):

    rad = (2 * math.pi) * perc / 100

    return rad


def rad2perc(rad):

    perc = 100 * rad / (2 * math.pi)

    return perc