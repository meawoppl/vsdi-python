from pylab import *
from scipy import *
from scipy import integrate

# t  = linspace(0,10*pi,1000)
# tc = sin(t)
# tc = 1.0 * (tc > 0)

# space_time_course = sig[:,:,newaxis] * tc


def two_dee_normed_gaussian(sigma_in_px, size):
    assert len(size) == 2,  "Read the function label jerk"
    assert size[0] == size[1],  "Only square ATM"
    assert (size[0] % 2) == 1, "Probably should be odd"
    
    # Xs and ys
    xs, ys = mgrid[-size[0]//2:size[0]//2, -size[1]//2:size[1]//2]
    gauss = (sigma_in_px * sqrt(2*pi))**-1 * exp(-0.5 * (xs**2 + ys**2) / sigma_in_px**2) 

    # Normalize so that it has a mean of 1 under the curve
    gauss /= gauss.sum()

    return gauss

# Make something semi-realistic
times  = linspace(0,1.2, 132)
on_off = (sin(times*2*pi) < 0.0) * 1.0

def dV_dt(V, t, I=on_off, It=times, g0=1, C=1, G_sigma=10, B_sigma=20):
    ''' 
    This is the equation describing the population gain control model.
    This gets integrated by odeint to determine (n-dimensional) response 
    characters to a stimulus in time:
      I is a 1d array describing the behavior of the signal in time It
      g0 is a constant related to the equalibrium condition
      G_sigma is the standard devation of the cerebral PSF.
      B_sigma is the standarrd deviation of the normalizing pool.
    '''

    # ODEint expects a 1d vector, so (TODO: Ugleeee)
    edge_width = int(round(sqrt(V.size)))
    V = V.reshape((edge_width, edge_width))

    # Figure out if the stimulus is on or off, or transitioning
    Ic = interp(t, I, It)

    # Make the smoothed stimulus 
    G = two_dee_normed_gaussian(G_sigma, V.shape)

    # Make the normalization pool gaussian
    B = two_dee_normed_gaussian(B_sigma, V.shape)

    # Compute : C (dv/dt) = I*G + (g0(1+I*B)) V
    dv_dt = ( (Ic*G) + (g0 * (1 + (Ic*B))) * V ) / C

    # Squish it back to the vector
    return dv_dt.flatten()


times_of_interest = linspace(0.0, 1.0, 100)
starting = zeros((101, 101)).flatten()

t = integrate.odeint(dV_dt, starting, times_of_interest)
