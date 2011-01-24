from tables import *
from ECG_reader import ECG_Decode

# Open the file and print the signal numbers
ecg_data = ECG_Decode("../Frodo20100803.h5")
ecg_data.print_event_count()
#ecg_data.plot_event_timeline()

# Event 12 starts, we grab 400ms before, and event 11 ends, 600 ms after also!
ecg_data.divy_experients(12, 400, 11, 600)

# Signal 40 (Break fixation) or less than 1000 ms exp length is no good!
ecg_data.mark_trials(40, 1000)  

# Mark the peaks in regions above 3000
ecg_data.mark_peaks_of_trace(absolute = 3000)

start_times = ecg_data.exp_starts


good_exps = start_times[ecg_data.good_trials]

print len(good_exps)

print len(ecg_data.h5.listNodes("/blk"))

intens = []
t_intens = []
for n, b in enumerate(ecg_data.h5.listNodes("/blk")):
    print n, b.name
    i = b[:,:,0].mean()
    t_intens.append(b[:].mean(axis=0).mean(axis=0))
    intens.append( i )
    print i



intens = array(intens)
figure("Mean Intensity")
plot(good_exps, intens)




