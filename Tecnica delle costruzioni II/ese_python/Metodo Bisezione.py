import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import pylab

def h(x):
  return x**2

class Bisezione():
  def __init__(self, f, y, a, b):
    an = a
    bn = b
    fcn = f((an+bn)/2) - y

    tau = 0.001

    while abs(fcn) > tau:
      cn = (an+bn)/2
      fan = f(an) - y
      fbn = f(bn) - y
      fcn = f(cn) - y
      if fan*fcn < 0:
        an = an
        bn = cn
      if fbn*fcn < 0:
        an = cn
        bn = bn
    self.cn = cn

a = 0
b = 1e4
y = 16
bisezione = Bisezione(h, y, a, b)
print bisezione.cn