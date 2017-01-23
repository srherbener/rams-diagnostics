#!/usr/bin/env python3


def f(x):
    return x**4 - 2*x + 1

N = 100
a = 0.0
b = 2.0
h = (b - a)/N

s = f(a) + f(b)

s_odd = 0.0
for k in range(1, N, 2):
    s_odd += f(a + k*h)
s_odd *= 4.0

s_even = 0.0
for k in range(2, N, 2):
    s_even += f(a + k*h)
s_even *= 2.0

result = (s + s_odd + s_even) * h / 3.0

print(result)
