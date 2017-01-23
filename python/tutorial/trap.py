#!/usr/bin/env python3

def f(x):
    return x**4 - 2*x + 1

N = 100
a = 0.0
b = 2.0
h = (b - a)/N

s = 0.5 * (f(a) + f(b))

for k in range(1, N):
    s += f(a+k*h)

result = h * s

print(result)
