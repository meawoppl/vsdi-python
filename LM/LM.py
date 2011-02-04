from scipy import *
from scipy import optimize, ndimage
from warnings import warn


def calc_fq(times, pulse_delay, on_time, off_time, pulse_number):
    '''This develops a group of baisis functions somewhat like our
       oscilating stimuli.  All inputs in seconds'''
    # Basically, this is the standard experiment type we run.
    #                      Experiment 
    # /=======================^===========================================\
    # <--Delay--><--on--><--off--><--on--><--off--><--on--><--off--> . . . 

    answers = zeros_like(times)

    # Calculate the "on" windows in seconds
    openings = []               # Start of the n'th window
    closings = []               # End   of the n'th window

    # Populate the above arrays
    for pn in range(pulse_number):
        offset = pulse_delay + ((on_time + off_time) * pn)
        openings.append(offset)
        closings.append(offset + on_time)

    # Based on the start and stop times, alter the signal
    for opening_time, closing_time in zip(openings, closings):
        during_this_pulse = (times > opening_time) & (times < closing_time)
        answers[during_this_pulse] = 1

    return ndimage.gaussian_filter1d(answers, 1)


def calc_ra(times, a1, a2, a3, a4):
    # Sums of times, basically, times at which the funciton branches
    # Calculating this upfront minimizes later confusion
    a12   = a1   + a2
    a123  = a12  + a3
    a1234 = a123 + a4

    #group1 = times <= a1
    group2 = (a1 < times) & (times <= a12)
    group3 = (a12 < times) & (times <= a123)
    group4 = (a123 < times) & (times <= a1234)
    print "g2:%i, g3:%i, g4:%i" % ((group2*1).sum(), (group3*1).sum(), (group4*1).sum())
    #group5 = (a1234 < times)

    # Somewhere to collect the answers
    answers = zeros_like(times)

    # Commented out the group1 and group5 (they are already zero)
    # But left them in to mimic the functional behavior

    # answers[group1] = 0
    answers[group2] = (0.5*(1-cos(pi*((times[group2] - a1[group2])/a2[group2]))))
    answers[group3] = 1
    answers[group4] = (0.5*(1+cos(pi*((times[group4]-a1[group4]-a2[group4]-a3[group4])/a4[group4]))))
    # answers[group5] = 0    
    return answers

def est_exp_decay(time_vector, signals, hb_in_hz=2):
    # Time x Signal-Number
    # TODO, docsting and assert docs
    assert signals.ndim == 2
    assert signals.shape[0] == time_vector.size

    # First_frame normalize the signal
    first_frame = signals[0,:]
    sig_norm = (signals - first_frame) / first_frame

    def sum_sq_err(parms):
        a, b, c, tau, phi = parms
        fit_fun = a * (exp(-time_vector/tau) - 1) + b * sin(2 * time_vector - phi) + c
        err = (sig_norm - fit_fun[:,newaxis])**2
        print a, b, c, tau, phi, err.sum()
        return err.sum()

    c_guess = array([0.01,0,sig_norm.mean(),1,0])

    return optimize.fmin(sum_sq_err, c_guess)


class LM_Decomposer:
    def __init__(self, trial_time_array):
        # Stupid check myself, before . . . 
        assert trial_time_array.min() == 0, "Must have a zero time point"
        assert trial_time_array.ptp() >  0, "Time must progress forward!"

        # Save the passed time array, and calculate its span
        self.times = trial_time_array
        self.time_span = trial_time_array.ptp()
        
        # This is where we will gather arrays (length time) that
        # Represent the various different compoennts we will 
        # Decompose the system into
        self.signals = []
        self.noises  = []
        self.noise_t = []

        # Get the constant signal out of the way
        self.constant = ones_like(trial_time_array)

        # Make sum to one, so coeff on this component is mean signal
        self.constant /= self.constant.sum() 

    def add_signal(self, sig):
        # Stupid check each signal before we add it.  
        # Must be the same length as the time array
        assert sig.size == self.times.size, "Must be same length as time array"
        self.signals.append(sig)

    def add_noise(self, sig, typ):
        # Like signals, but label the types to do stuff later
        assert sig.size == self.times.size, "Must be same length as time array"
        self.noises.append(sig)
        self.noise_t.append(typ)

    def add_sinusoid(self, freq_in_hz):
        # Both sine and cosine are added 
        # The linear combintion deals with the phase shift)
        self.add_noise(sin(2*pi*freq_in_hz * self.times), typ="sin")
        self.add_noise(cos(2*pi*freq_in_hz * self.times), typ="cos")

    def add_our_signal(self, seomfhti):
        '''This adds a signal more like the signals that we see in our experiments'''
        pass

    def add_signal_group(self, a_mins, a_maxes, a_steps):
        ''' This adds a signal group like in the paper.  See figure 2 c&d.'''
        assert c_[a_mins, a_maxes, a_steps].shape == (4,3), "Not the required shape.  See doc."
        assert a_mins.size == 4, "Not the required shape.  See doc."
        assert "int" in str(a_steps.dtype), "Number of steps must be an integer."
        # Were gonna explode the start-stop-step data into a tuple of slices
        # From there, we can make mgrid generate all possabilites
        # And compute all timecourse in one function call

        print "Computing all %i possibilities." % (a_steps.prod() * self.times.size)

        # Transform the start stop and step sizes into "slices"
        slices = []
        for start, stop, step in zip(a_mins, a_maxes, a_steps):
            slices += [slice(start, stop, step*1j)] 
                  
        # Add a slice for time onto the front
        slices = [slice(self.times.min(), self.times.max(), self.times.size * 1j)] + slices       

        # Now each of these is a (t x A1 x A2 x A3 x A4) array of values
        # Spanning the slices defined above
        ts, a1s, a2s, a3s, a4s = mgrid[slices]

        # Flatten to time x number of signals (A1 * A2 * A3 * A4)
        self.sofar = array( calc_ra(ts, a1s, a2s, a3s, a4s) ).reshape((self.times.size, -1))

        # Decompose that into the baisis vectors (v is what we care about)
        u, s, v = svd(lmd.sofar.T)

        # Grab the top 10 and add them as signals
        cutoff_n = 10
        for x in range(cutoff_n):
            self.add_signal( v[x,:] )


    def set_exp_decay(self, time_const):
        # Semi-stupid check
        if hasattr(self, "exp_signal"):
            warn("Set decay exponent called more than once.  Deleting previous definition")
        self.exp_signal = "Set"

        # Add the exp falloff as a noise 
        self.add_noise(exp(-times/time_const) - 1., "exp")


if __name__ == "__main__":
    # Now lets average the exp. falloff of our dataset.  
    import tables
    from pylab import *


    exp_times = linspace(0, 1200, 132)
    lmd = LM_Decomposer(exp_times)

    # # This is just like from the paper to confirm that the
    # # SVD vector extraction is working as expected
    # a_mins  = array([20,  30,  450, 60], dtype=float64)
    # a_maxs  = array([120, 200, 590, 100], dtype=float64)
    # a_steps = array([11,  18, 15, 4], dtype=int32)
    # lmd.add_signal_group(a_mins, a_maxs, a_steps)

    
    h5 = tables.openFile("/media/Local Disk/Frodo20100716-run0.h5")
    datae = h5.root.block[200:300,200:300,:,:]

    # Collapse the first two space axes
    datae = datae.mean(axis=0).mean(axis=0)

    conditions = h5.root.condition[:]

    uniq_conditions = set(conditions)
    print uniq_conditions

    for u in uniq_conditions:
        current_condition = (conditions == u)
        trace = datae[:,current_condition].mean(axis=1)
        normed_trace = (trace - trace[0])/trace[0]
        plot(normed_trace, label=str(u))



    # a, b, c, tau, f, phi = parms
    # fit_fun = a * (exp(-time_vector/tau) - 1) + b * sin(f * time_vector - phi) + c



    legend()
    show()

    

    
