#!/usr/bin/env python3

import sys

def GetArgs(ArgList):
    from optparse import OptionParser

    # define the parser
    Usage = "usage: %prog [options]"
    parser = OptionParser(Usage)
    parser.add_option("-n", "--numsteps", dest="numsteps", 
                      type="int", default=10,
                      help="Number of steps")

    # parse the argument list
    (options, args) = parser.parse_args()

    # check options and arguments
    if (options.numsteps < 0):
        parser.error("Option -n requires positive integer: {0:10d}".format(options.numsteps))

    # send back option values
    return ( options, args )


def f(x):
    return x**4 - 2*x + 1


def trap(N, a, b, h):
    s = 0.5 * (f(a) + f(b))

    for k in range(1, N):
        s += f(a+k*h)

    return h * s


def simp(N, a, b, h):
    s = f(a) + f(b)

    s_odd = 0.0
    for k in range(1, N, 2):
        s_odd += f(a + k*h)
    s_odd *= 4.0

    s_even = 0.0
    for k in range(2, N, 2):
        s_even += f(a + k*h)
    s_even *= 2.0

    return (s + s_odd + s_even) * h / 3.0


# MAIN

(options, args) = GetArgs(sys.argv)

N = options.numsteps
a = 0.0
b = 2.0
h = (b - a)/N

result_trap = trap(N, a, b, h)
result_simp = simp(N, a, b, h)

print("Number of steps = {0:10d}".format(N))
print("")

print("Result from trapezoid method = {0:20.10f}".format(result_trap))
print("Result from Simpson's method = {0:20.10f}".format(result_simp))
