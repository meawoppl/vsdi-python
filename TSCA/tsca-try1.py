from scipy import *
from scipy import signal
from pylab import *

def make_autocorrelation_matrix(signal1):
    mat = signal1[:,newaxis] * signal1[newaxis,:]
    mat /= mat.max()
    return mat

size_of_image = (30, 30)
time_slices = arange(0,1000)

x, y = mgrid[-1.0:1.0:size_of_image[0]*1j, -1.0:1.0:size_of_image[1]*1j]

# The Geometric Shapes
in_circle = 1.0 * (sqrt(x**2 + y**2) < 0.5)
x_grating = cos(x * 2 * pi)

# Their individual time Courses
circ_strength = cos(time_slices * 0.21342 / (2 * pi))
wave_strength = cos(time_slices * 0.13420 / (2 * pi))

# X x Y x T
time_course = in_circle[:,:,newaxis] * circ_strength + x_grating[:,:,newaxis] * wave_strength

# Correlation matrices
Cx = make_autocorrelation_matrix(circ_strength)
Cy = make_autocorrelation_matrix(wave_strength)

gamma1, gamma2 = 1, 1

intermediate_matrix = matrix([[sum(Cx*Cx), sum(Cy*Cx)],
                              [sum(Cx*Cy), sum(Cy*Cy)]])

vec = array([[gamma1*trace(Cx)],
             [gamma2*trace(Cx)]])

alphae = intermediate_matrix.I * vec
alpha1, alpha2 = alphae.flat

Q = matrix( (Cx * alpha1) + (Cy * alpha2) )
Z = matrix( time_course.reshape((size_of_image[0]*size_of_image[1], -1)) )

vals, vecs = eig(Z*Q*Z.T)

vals_as = vals.argsort()

comp1 = real(vecs[:,0]).reshape(size_of_image)
comp2 = real(vecs[:,1]).reshape(size_of_image)

imshow(comp1)
show()

savetxt("c1.csv", comp1)
savetxt("c2.csv", comp2)
