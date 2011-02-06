from scipy import *
from scipy import optimize
from warnings import warn
from do_once import DoOnce

def blank_fit_function(time_vector, a, b, c, tau, phi, hb_in_hz):
    return a * (exp(-time_vector/tau) - 1) + b * sin(2 * pi * hb_in_hz * time_vector - phi) + c

def est_exp_decay(time_vector, mean_signal, debug=False):
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

        # # Handy Fit Debugging statement
        if debug: print a, b, c, tau, phi, hb_in_hz, err.sum()
        return err.sum()

    # Starting guess . . . 
    c_guess = array([0.01,0, mean_signal.mean(), 1., 0., 2.])

    # Run the Simplex Optimizer
    a, b, c, tau, phi, hb_in_hz = optimize.fmin(sum_sq_err, c_guess)

    # Stupid check
    assert (hb_in_hz < 4.0), "Heartbeat fit returned as > 4 Hz (%f)" % hb_in_hz
    assert (hb_in_hz > 1.0), "Heartbeat fit returned as < 1 Hz (%f)" % hb_in_hz
    
    # Pass back found values
    return a, b, c, tau, phi, hb_in_hz 


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


def graph_hb_trend(h5_file, img_file):
    # Open the hdf5 file and get the centerblock trends
    cn_to_tc = compute_mean_centerblock_conditions(h5_file, do_debug=True)

    # TODO!  Unhardcode this!
    exp_times = linspace(0, 1.2, 132)

    figure()
    for cn, ntr in cn_to_tc.iteritems():
        # If it is the blank condition
        if cn == 0: 
            # Fit the mean blank with exp sinusoidy thing (get the paramaters)
            a, b, c, tau, phi, hb_in_hz = est_exp_decay(exp_times, ntr)

            # Compute the fit shape from those paramaters
            fit_vals = blank_fit_function(exp_times, a, b, c, tau, phi, hb_in_hz)
            
            # Plot the fit, in bark black with a dashed link
            plot(exp_times, fit_vals, "black", linestyle="--", alpha=0.9,  label="Blank Fit")

            # Make the blank conditions stand out
            plot_alpha = 1.0
        else:
            # But still plot the non-blank
            plot_alpha = 0.3

        # This plots the data for both the above with different alpha values
        plot(exp_times, ntr, label="Cond. "+str(cn), alpha=plot_alpha)

    legend()

    # Print out the fit equation on the plot
    equation_text = r"FIT = $a \, e^{(-t/\tau) - 1} + b \, \sin(2 \pi f \, t + \phi) +c$"

    # These describe the varible values used
    var_text1 = r"$a=%.03E \,\,\,\, \tau=%.03f \,\,\,\,  b=%.03E$" % (a, tau, b)
    var_text2 = r"$f=%.03f \,\,\,\, \phi=%.03E \,\,\,\,  c=%.03E$" % (hb_in_hz, phi, c)

    # Lower left cornerish
    text(0.05, -0.0045, equation_text)
    text(0.05,-0.0047, var_text1)
    text(0.05,-0.0049, var_text2)

    # Title and axis labels
    title("Center Block Mean per Condiiton (Frodo20100716-run0)")
    xlabel("Time in Seconds")
    ylabel("$\Delta f/f$")

    savefig(img_file)


if __name__ == "__main__":
    # Now lets average the exp. falloff of our dataset.  
    import tables, sys

    # Stupid check the CLI information
    assert len(sys.argv) == 3, "\nUsage: python %s [hdf5-file] [output-image]" % sys.argv[0]
    from pylab import *

    graph_hb_trend(sys.argv[1], sys.argv[2])
