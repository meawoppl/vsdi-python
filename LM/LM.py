from scipy import *
from scipy import optimize
from warnings import warn

def calc_ra(times, a1, a2, a3, a4):
    a12 = a1 + a2
    a123 = a12 + a3
    a1234 = a123 + a4

    #group1 = times <= a1
    group2 = (a1 < times) & (times <= a12)
    group3 = (a12 < times) & (times <= a123)
    group4 = (a123 < times) & (times <= a1234)
    print "g2:%i, g3:%i, g4:%i" % ((group2*1).sum(), (group3*1).sum(), (group4*1).sum())
    #group5 = (a1234 < times)

    answers = zeros_like(times)

    # answers[group1] = 0
    answers[group2] = (0.5*(1-cos(pi*((times[group2] - a1[group2])/a2[group2]))))
    answers[group3] = 1
    answers[group4] = (0.5*(1+cos(pi*((times[group4]-a1[group4]-a2[group4]-a3[group4])/a4[group4]))))
    # answers[group5] = 0    
    return answers


def est_exp_decay(time_vector, signals):
    # Time x Signal-Number
    # TODO, docsting and assert docs
    assert signals.ndim == 2
    assert signals.shape[0] == time_vector.size

    # First_frame normalize the signal
    first_frame = signals[0,:]
    sig_norm = (signals - first_frame) / first_frame

    def sum_sq_err(parms):
        a, b, c, tau, f, phi = parms
        fit_fun = a * (exp(-time_vector/tau) - 1) + b * sin(f * time_vector - phi) + c
        err = (sig_norm - time_vector[:Newaxis])**2
        return err.sum()

    c_guess = ones(6)

    print optimize.fmin(sum_sq_err, c_guess)


class LM_Decomposer:
    def __init__(self, trial_time_array):
        # Stupid check myself, before . . . 
        assert trial_time_array.min() == 0, "Must have a zero time point"
        assert trial_time_array.ptp() >  0, "Time must progress forward!"

        self.times = trial_time_array
        self.time_span = trial_time_array.ptp()
        
        self.sinusoids = []
        self.cosusoids = []
        
        # Get the constant signal out of the way
        self.constant = ones_like(trial_time_array)
        # Make sum to one, so coeff on this component is mean signal
        self.constant /= self.constant.sum() 

    def add_sinusoid(self, freq_in_hz):
        self.sinusoids.append(sin(2*pi*freq_in_hz * self.times))
        self.cosusoids.append(cos(2*pi*freq_in_hz * self.times))

    def add_signal_group(self, a_mins, a_maxes, a_steps):
        assert c_[a_mins, a_maxes, a_steps].shape == (4,3), "Not the required shape.  See doc."
        assert a_mins.size == 4, "Not the required shape.  See doc."
        assert "int" in str(a_steps.dtype), "Number of steps must be an integer."
        # Were gonna explode the start-stop-step data into a tuple of slices
        # From there, we can make mgrid generate all possabilites
        # And compute all timecourse in one function call

        print "Computing all %i possabilities." % (a_steps.prod() * self.times.size)

        slices = []
        for start, stop, step in zip(a_mins, a_maxes, a_steps):
            slices += [slice(start, stop, step*1j)] 
                  
        slices = [slice(self.times.min(), self.times.max(), self.times.size * 1j)] + slices       
        ts, a1s, a2s, a3s, a4s = mgrid[slices]

        self.sofar = array( calc_ra(ts, a1s, a2s, a3s, a4s) ).reshape((self.times.size, -1))

    def set_exp_decay(self, time_const):
        if hasattr(self, "exp_trend"):
            warn("Set decay exponent called more than once.  Deleting previous definition")
        self.exp_signal = exp(-times/time_const) - 1.




if __name__ == "__main__":
    # Now lets average the exp. falloff of our dataset.  
    import tables
    from pylab import *

    a_mins  = array([20,  30,  450, 60], dtype=float64)
    a_maxs  = array([120, 200, 590, 100], dtype=float64)
    a_steps = array([11,  18, 15, 4], dtype=int32)

    exp_times = linspace(0, 1200, 132)
    lmd = LM_Decomposer(exp_times)

    lmd.add_signal_group(a_mins, a_maxs, a_steps)

    
    h5 = tables.openFile("/media/Local Disk/Frodo20100716-run0.h5")

    datae = h5.root.block[200:300,200:300,:,:]

    # Collapse the first two space axes
    datae = datae.mean(axis=0).mean(axis=0)

    conditions = h5.root.condition[:]

    uniq_conditions = set(conditions)

    for u in uniq_conditions:
        current_condition = (conditions == u)
        plot(datae[:,current_condition].mean(axis=1), label=str(u))

    show()

    
