def blank_fit_function(time_vector, a, b, c, tau, phi, hb_in_hz):
    return a * (exp(-time_vector/tau) - 1) + b * sin(2 * pi * hb_in_hz * time_vector - phi) + c

def est_exp_decay(time_vector, mean_signal):
    # TODO, docsting and assert docs
    assert mean_signal.shape == time_vector.shape
    # First_frame normalize the signal
    first_frame = mean_signal[0]
    # This little function returns the sum of squared errors
    def sum_sq_err(parms):
        # Asplode the array (vector required for the optimizer)
        a, b, c, tau, phi, hb_in_hz = parms
        
        # Calculare the expected values
        expected_values = blank_fit_function(time_vector, a, b, c, tau, phi, hb_in_hz)

        # Sum of squared errors
        err = (mean_signal - expected_values)**2

        # TODO: Debugging statement
        print a, b, c, tau, phi, hb_in_hz, err.sum()
        return err.sum()

    # Starting guess . . . 
    c_guess = array([0.01,0, mean_signal.mean(), 1., 0., 2.])
    return optimize.fmin(sum_sq_err, c_guess)


@DoOnce
def compute_mean_centerblock_conditions(hdf5_filename):
    '''This functions takes an hdf5 file as an argument, and 
    returns a dictionary of (center averaged) time courses 
    keyed to the experiment condition'''

    # Open a hdf5 object of the give filename
    h5 = tables.openFile(hdf5_filename)

    # Pull out the middle of the block for all times and all trials
    datae = h5.root.block[200:300,200:300,:,:]

    # Collapse the first two space axes remainder is (time x trials)
    datae = datae.mean(axis=0).mean(axis=0)
    
    # Grab the conditions out of the array, and find the unique ones
    conditions = h5.root.condition[:]
    uniq_conditions = set(conditions)

    # Define a dict to collect the datas in 
    cn_to_meantrace = {}

    # For each condition, collapse on trial
    for u in uniq_conditions:
        # Which are this condition number
        current_condition = (conditions == u)

        # Make an array that are only that condition and average all trials
        trace = datae[:,current_condition].mean(axis=1)

        # Make it dF/F and record
        cn_to_meantrace[u] = (trace - trace[0])/trace[0]
    
    return cn_to_meantrace

if __name__ == "__main__":
    # Now lets average the exp. falloff of our dataset.  
    import tables
    from pylab import *
    
    # Laptop
    cn_to_tc = compute_mean_centerblock_conditions("/home/meawoppl/Frodo20100716-run0.h5", do_debug=True)

    # # Office
    # cn_to_tc = compute_mean_centerblock_conditions("/media/Local Disk/Frodo20100716-run0.h5")

    exp_times = linspace(0, 1.2, 132)
    lmd = LM_Decomposer(exp_times)

    # # This is just like from the paper to confirm that the
    # # SVD vector extraction is working as expected
    # a_mins  = array([20,  30,  450, 60], dtype=float64)
    # a_maxs  = array([120, 200, 590, 100], dtype=float64)
    # a_steps = array([11,  18, 15, 4], dtype=int32)
    # lmd.add_signal_group(a_mins, a_maxs, a_steps)

    figure()
    for cn, ntr in cn_to_tc.iteritems():
        plot(exp_times, ntr, label=str(cn), alpha=0.8)

        if cn == 0:
            parms = est_exp_decay(exp_times, ntr)
            a, b, c, tau, phi, hb_in_hz = parms
            # hb_in_hz = 2.

            fit_vals = blank_fit_function(exp_times, a, b, c, tau, phi, hb_in_hz)
            plot(exp_times, fit_vals, "black", linestyle="--", alpha=0.9,  label="Blank Fit")

    legend()
    equation_text = r"FIT = $a \, e^{(-t/\tau) - 1} + b \, \sin(2 \pi f \, t + \phi) +c$"

    var_text1 = r"$a=%.03E \,\,\,\, \tau=%.03f \,\,\,\,  b=%.03E$" % (a, tau, b)
    var_text2 = r"$f=%.03f \,\,\,\, \phi=%.03E \,\,\,\,  c=%.03E$" % (hb_in_hz, phi, c)
    text(0.05, -0.0045, equation_text)
    text(0.05,-0.0047, var_text1)
    text(0.05,-0.0049, var_text2)

    title("Center Block Mean per Condiiton (Frodo20100716-run0)")

    xlabel("Time in Seconds")
    ylabel("$\Delta f/f$")

    savefig("blank-fitting-params.pdf")

    show()
