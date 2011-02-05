from pylab import *
from scipy import *
from scipy import integrate

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
times  = linspace(0,1.0, 1000)
on_off = (sin(times*10*pi) < 0.0) * 1.0

def dV_dt(V, t, I=on_off, It=times, g0=1, C=1, k=1e5, A_sigma=10, B_sigma=20):
    ''' 
    This is the equation describing the population gain control model.
    This gets integrated by odeint to determine (n-dimensional) response 
    characters to a stimulus in time:
      I is a 1d array describing the contrast of the signal in time 
      It is the time course that matches the above (units seconds)
      g0 is a constant related to the equalibrium condition
      G_sigma is the standard deviation of the cerebral PSF.
      B_sigma is the standard deviation of the normalizing pool.
    '''

    # ODEint expects a 1d vector, so (TODO: Ugleeee)
    edge_width = int(round(sqrt(V.size)))
    V = V.reshape((edge_width, edge_width))

    # Figure out if the stimulus is on or off, or transitioning
    Ic = interp(t, It, I)

    # Make the smoothed stimulus 
    A = two_dee_normed_gaussian(A_sigma, V.shape) * Ic

    # Make the normalization pool gaussian
    B = two_dee_normed_gaussian(B_sigma, V.shape) * Ic * k

    # Relate this to RC circuit as the following
    #    dV     V
    # C ---- + ---  = 0
    #    dt     R
    #                      OR . . .
    #    dV       /  V  \
    #   ---- = - | ----  |
    #    dt       \ RC  /

    # Then we thing about conductance instead of resistance (1/R -> N)

    #    dV       / NV  \
    #   ---- = - | ----  |
    #    dt       \  C  /

    # But our conductor has a baseline and is a function local pooled 
    # Properties as: 
    # dv_dt =  10 * ((Ic + V_resting) - V)

    dv_dt = A - (g0 * (1+B) * V)

    # Squish it back to the vector
    return dv_dt.flatten()

orig_shape = (101, 101)
times_of_interest = linspace(0.0, 1.0, 1000)
starting = zeros(orig_shape).flatten()

results = integrate.odeint(dV_dt, starting, times_of_interest, hmax=0.01)
results = results.reshape(tuple([-1]) + orig_shape)

close("all")

figure()
subplot(2,1,1)
# Time Course of the stimulus

plot(times, on_off)

ylabel("Stimulus On-Off")
axis((0.0, 1.0, 0.0, 1.2))

subplot(2,1,2)
# TC of the response
plot(times_of_interest, results[:,50,50]*1000)

xlabel("Time in Seconds")
ylabel("Center Point Response")

savefig("hells-yea.pdf")
show()
