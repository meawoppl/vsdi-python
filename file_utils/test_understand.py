from scipy import *
from pylab import *


xs = linspace(-10, 10, 1000)
sig1 = sin(pi * xs) + 0.1 * normal(size=1000)
sig2 = sin(pi * xs + pi) + 0.1 * normal(size=1000)

figure()
plot(xs, sig1)
plot(xs, sig2)

figure()
psd(sig1, Fs=10)
psd(sig2, Fs=10)

show()






