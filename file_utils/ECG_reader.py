from scipy import ndimage
from pylab import *
from tables import openFile

# import matplotlib, sys

# matplotlib.rcParams['font.family'] = 'serif'
# matplotlib.rcParams['font.size'] = 12
# matplotlib.rcParams['axes.labelsize'] = 14
# matplotlib.rcParams["text.usetex"] = True

# Event definitions from .log file
# TRIAL_START=10
# TRIAL_END=11
# DB_TRIGGER=12
# OI_TRIG_ON=30
# OI_TRIG_OFF=31
# OI_GO=32
# BREAK_FIX=40
# SACCADE_START=41
# SACCADE_END=42
# REWARD=60
# NO_REWARD=61
# FP_ON=70
# FP_OFF=71
# STIM_ON=72
# STIM_OFF=73
# T1_ON=74
# T1_OFF=75
# T2_ON=76
# T2_OFF=77

num_to_event = {
    0:"NO_EVENT",
    10:"TRIAL_START",
    11:"TRIAL_END",
    12:"DB_TRIGGER",
    30:"OI_TRIG_ON",
    31:"OI_TRIG_OFF",
    32:"OI_GO",
    40:"BREAK_FIX",
    41:"SACCADE_START",
    42:"SACCADE_END",
    60:"REWARD",
    61:"NO_REWARD",
    70:"FP_ON",
    71:"FP_OFF",
    72:"STIM_ON",
    73:"STIM_OFF",
    74:"T1_ON",
    75:"T1_OFF",
    76:"T2_ON",
    77:"T2_OFF"}

class openECG():
    def __init__(self, h5_exp_file):
        self.h5 = openFile(h5_exp_file)
        self.events_array = self.h5.root.htb.Events[:]

        # Find all the non-zero's in the events array
        # These are the actual events
        self.non_zero_event_times = self.events_array.nonzero()[0]
        self.non_zero_event_numbers = self.events_array[ self.non_zero_event_times ]
        
    def print_event_count(self):
        # Print Unique event numbers (helps to understand what is going on)
        unique_event_numbers = list( set(self.non_zero_event_numbers) )
        unique_event_numbers.sort()
        for u in unique_event_numbers:
            print "Event #%i - (%s) was recorded %i times" % (u, num_to_event[u], (1*(self.non_zero_event_numbers==u)).sum() )

    def plot_event_timeline(self):
        # Using latex labeling on this takes _forever_
        matplotlib.rcParams["text.usetex"] = False
        figure()
        title("Events and Event Descriptions")

        # Plot the actual values of the signals 
        plot(self.non_zero_event_times, self.non_zero_event_numbers)
        
        # Label each point with a name corresponding to the signal
        for i, n in zip(self.non_zero_event_times, self.non_zero_event_numbers):
            event_label = str(n) + " - " + num_to_event[n]
            text(i, n, event_label)

        show()

    def get_event_times(self, event_number):
        return self.non_zero_event_times[ self.non_zero_event_numbers==event_number ]

    def divy_experients(self, start_event, time_before, stop_event, time_after):
        # This value is noted so it can be used to compute 
        # spike locations with respect to imaging start times
        self.experiment_offset = time_before

        # We also store the global ekg mean
        self.global_ekg_max = -inf

        # This divides up the EKG traces based on event log
        # It makes no assessment of whether it was a valid experiment!
        self.exp_starts = self.get_event_times( start_event ) # 12
        self.exp_stopps = self.get_event_times( stop_event  ) # 11

        # Stupid check!
        assert len(self.exp_starts) == len(self.exp_stopps), "There must be the same number of start and stop signals!"
        
        # Step over the trials, and rip out the ECG data
        self.traces = []
        for n, (st, et) in enumerate(zip(self.exp_starts, self.exp_stopps)):
            # Compute the time to start and end cutting
            start_index = st - time_before
            end_index   = et + time_after

            # Make sure were not wrapping . . . 
            if start_index < 0: 
                print "Index corrected . . .", start_index
                start_index = 0

            # Cut the signal out of the block
            ekg_trace = self.h5.root.htb.EKG[ start_index : end_index ]
            self.traces.append(ekg_trace)

            # Find its max, and compare to the global max.  Store if larger
            trace_max = ekg_trace.max()
            if trace_max > self.global_ekg_max:
                self.global_ekg_max = trace_max


        # Note the number or trials
        self.ekg_count = len(self.traces)

    def mark_trials(self, abort_event_number, abort_length):
        no_good_trials = []
        for n, (st, et) in enumerate(zip(self.exp_starts, self.exp_stopps)):
            # Check to see if there was a fix break inside the trial
            #                   v   after start        v       V     before end     V 
            in_time_bracket = (self.non_zero_event_times > st) & (self.non_zero_event_times < et )
            abort_event     = (self.non_zero_event_numbers == abort_event_number )

            # Mark trials with an abort event inside the trigger windows
            if any(in_time_bracket & abort_event):
                no_good_trials += [n]
                print "Trial number", n, "is no good!"
                continue

            # Mark trials that are too short . . . necessary?
            # I don't understand why these experiments get recorded, but they do!
            if (et - st < abort_length):
                no_good_trials += [n]
                print "Trial number", n, "is no good too short!  Only",  et - st, "ms"
                continue

        # Preen out the trials that at no good . . .
        # Make some structures to facilitate slicing
        self.bad_trials = zeros(self.ekg_count)
        self.bad_trials[array(no_good_trials)] = 1
        self.good_trials = logical_not( self.bad_trials )

    def mark_peaks_of_trace(self, absolute=None, relative=None):
        '''Absolute means that the peak-finding is based on the raw value.  
           Relative is a fraction of the maximum. '''
        # Calculate the cutoff
        # The logic is inflated for sanity check 
        if (absolute is not None) and (relative is None):
            ekg_cutoff = absolute
        elif (relative is not None) and (absolute is None):
            ekg_cutoff = self.global_ekg_max * relative
        else:
            # Probably not necessary, as ekg_cutoff would be undefined
            raise RuntimeError("Must sepcify ekg cutoff in one way!")
        
        # Stuff to store internally. 
        self.spike_times = []   # The times in EKG ms offset
        self.spike_number = []  # The number of spikes
        self.spike_heights = [] # The amplitudes

        # For each EKG trace
        for n, tr in enumerate(self.traces):
            # Find the area above the cutoff
            screen = (tr > ekg_cutoff)

            # Label each and record the number of labels
            labels, spike_count = ndimage.label(screen)
            
            # Find the extrema; locations and amplitudes
            indexs_of_spikes = range(1, spike_count + 1)

            # This conditional should not be necessary.  
            # TODO: scipy bug-tracker is down . . .  report this one
            if indexs_of_spikes != []:
                ext = ndimage.extrema(tr, labels=labels, index=indexs_of_spikes)
            else:
                ext = ([],[],[],[])
            
            # Break out the results of the extrema, ignore min-related values
            min_vals, max_vals, min_locs, max_locs = ext

            # Store the values on a per-experiment basis
            # The arrays appended _could_ be empty, but will be an object . . . 
            # so the offsets will still be correct.
            self.spike_times.append( array(max_locs).flatten() )
            self.spike_number.append( spike_count )
            self.spike_heights.append( array(max_vals).flatten() ) 

    # A bit obvious, just convenience functions
    def is_good_trial(self, trial_no):
        return self.good_trials[trial_no]
    def is_bad_trial(self, trial_no):
        return self.bad_trials[trial_no]

    def iterate_trials(self):
        # Iterator that has the form
        # image_data, condition, trace, peaks_times, peak_vals
        
        blk_node_list = [node.name for node in self.h5.listNodes("/blk")]
        blk_node_list.sort()

        # Sync up the trial and condition datae
        # TODO, handle this in the file loader/packer
        cond_sorter = self.h5.root.trial_num[:].argsort()
        conditions =  self.h5.root.condition[:][cond_sorter]

        blk_node_number = 0
        for tn, (trace, spi_t, spi_h) in enumerate( zip(self.traces, 
                                                        self.spike_times, 
                                                        self.spike_heights)):
            # Skip bad trials . . . 
            # They dont line up with the blk files, do don't increment
            if self.is_bad_trial(tn):
                continue

            # Grab the right image block and experimental condition
            blk_name = blk_node_list[blk_node_number]
            cond = conditions[blk_node_number]

            # Get the node of interest to this block
            image_data_node = getattr( self.h5.root.blk, blk_name)
            # return the node instead of the data in case we don't need this block
            # image_data = image_data_node[:]
            
            # Yield the results and increment the block file counter!  
            yield image_data_node, cond, trace, spi_t, spi_t
            blk_node_number += 1
        

# # Open the file and print the signal numbers
# ecg_data = ECG_Decode("../Frodo20100803.h5")
# ecg_data.print_event_count()
# ecg_data.plot_event_timeline()

# # Event 12 starts, we grab 400ms before, and event 11 ends, 600 ms after also!
# ecg_data.divy_experients(12, 400, 11, 600)

# # Signal 40 (Break fixation) or less than 1000 ms exp length is no good!
# ecg_data.mark_trials(40, 1000)  

# # Mark the peaks in regions above 3000
# ecg_data.mark_peaks_of_trace(absolute = 3000)




# good_spikes = [ecg_data.spike_times[n] for n in range(len(ecg_data.spike_times)) if ecg_data.is_good_trial(n) ]
# good_spikes = array(good_spikes)

# 1/0


# figure()
# title("Offset from start of ECG Recording")
# hist(good_spikes[:,0], bins=25)
# xlabel("Milliseconds")
# ylabel("Count")
# savefig("offset_prune.pdf")


# figure()
# for trace, spi_t, spi_h in zip(ecg_data.traces, 
#                                ecg_data.spike_times, 
#                                ecg_data.spike_heights):
#     plot(trace)
#     plot(spi_t, spi_h, "bo", alpha=0.6)
    
# exp_start_ms = ecg_data.experiment_offset
# exp_lenght_ms = int(134./110 * 1000)
# exp_end_ms = exp_start_ms + exp_lenght_ms

# # Label and color some interesting stuff
# text(100, 7500, "Pre-Imaging")
# axvspan(0, exp_start_ms, facecolor="r", alpha = 0.2)

# axvline([exp_start_ms], alpha=0.5, linestyle="--", color="black")

# text(1000, 7500, "Imaging")
# axvspan(exp_start_ms, exp_end_ms, facecolor="b", alpha = 0.2)

# axvline([exp_end_ms], alpha=0.5, linestyle="--", color="black")

# text(1850, 7500, "Post-Imaging")
# axvspan(exp_end_ms, exp_end_ms+600, facecolor="r", alpha = 0.2)

# xlabel("Time in Milleseconds")
# ylabel("Raw ECG Data")

# title("Reference Slide")

# axis((0.0, exp_end_ms+600, -4000.0, 8000.0))

# # savefig("hb-assembly-ref.pdf")

# figure()

# # Assemble mean heartbeat space blank
# interp_rez = 132
# heart_beat_video = zeros((504,504,interp_rez))
# blank_count = 0


# for img, cond, trace, p_time, p_vals in ecg_data.iterate_trials():
#     # Skip non blank trials for now . . .
#     if cond != 0 or p_time[0] > 330:
#         continue
#     print img
#     print p_time

#     # Grab the whole image, and deterend based on mean intensity
#     image_data = img[:]
#     frame_mean = image_data.mean(axis=0).mean(axis=0)
#     frames = arange(frame_mean.size)
#     slope, y_int = polyfit(frames, frame_mean, 1)
#     image_data = image_data - (frames * slope + y_int)

#     heart_beat_frames = around( ( ( array( p_time ) - 400) / 1000. ) * 110)
#     print heart_beat_frames

#     # For each pix.
#     hbs_interp = linspace(0.0, 1.0, interp_rez)

#     print heart_beat_frames[1], heart_beat_frames[2], heart_beat_frames[2] - heart_beat_frames[1]

#     # This nukes previous def.
#     # heart_beat_video = image_data[:, :, heart_beat_frames[1]:heart_beat_frames[2]]
#     heart_beat_video += image_data[:, :, :]

#     # for x in range(504):
#     #     print x
#     #     sys.stdout.flush()
#     #     for y in range(504):
#     #         # Grab that time course inside the heart beats
#     #         tc = image_data[x, y, heart_beat_frames[1]:heart_beat_frames[2]]

#     #         # Make a heartbeat space array to match
#     #         hbs = linspace(0.0, 1.0, tc.size)

#     #         # Iterpolate out 120 frames, and add to the video block
#     #         heart_beat_video[x, y, :] = interp( hbs_interp, hbs, tc)
  
#     blank_count += 1
#     break




# # Mean b/t blanks
# heart_beat_video /= blank_count

# heart_beat_video = heart_beat_video - heart_beat_video[:,:,0][:,:,newaxis]

# # 1/0

# figure(figsize=(4,4), dpi=128)
# # mx, mi = heart_beat_video[350:450,150:250,:].max(), heart_beat_video[350:450,150:250,:].min()
# mx, mi = heart_beat_video.max(), heart_beat_video.min()
# for frame in range(heart_beat_video.shape[2]):
#     frac_frame = 1.* frame / heart_beat_video.shape[2]
#     print "Rendering frame", frame, frac_frame
#     sys.stdout.flush()

#     clf()
#     imshow(heart_beat_video[:,:,frame], vmin=mi, vmax=mx)
#     # imshow(heart_beat_video[350:450,150:250,frame], vmin=mi, vmax=mx)
#     gray()
#     title("Frame: %i" % frame)
#     xlabel("Width in pix")
#     ylabel("Height in pix")
#     axis("image")
#     colorbar()

#     savefig("frames/1trialonly-%03i.png" % frame)

#     # plot( frame_mean )
#     # plot( frames * slope + y_int )

#     # print , p_time[2]
#     # sys.stdout.flush()

# # Project the blank back into experiment time between heartbeats
# # Crop that to the imaging time, and use in TSCA analysis

# 1/0

# # Heartbeat space figure
# trials = [name for name in dir(h5.root.blk) if not name.startswith("_")]

# # Stupid check!
# print len(trials), len(traces)
# assert len(trials) == len(traces), "We need to have as many experiments at ECG traces!"

# # Were gonna make a 504x504x100 block that is the mean heartbeat artifact . . .  maybe

# frames_per_second = 110
# time_per_slice = 1./frames_per_second
# interp_accumulator = []

# for spi, tri, tra in zip(spikes_in_exp, trials, traces):
#     # For now, only blank trials
#     if not tri.endswith("00"):
#         continue

#     hb1_ms, hb2_ms = spi

#     hb1_frame = int( round( (hb1_ms / 1000.) * frames_per_second ))
#     hb2_frame = int( round( (hb2_ms / 1000.) * frames_per_second ))

#     trial_data = getattr(h5.root.blk, tri)[:,:,hb1_frame-1:hb2_frame+1]
#     trial_data = trial_data.mean(axis=0).mean(axis=0)

#     fr_count = trial_data.shape[0]
#     fake_time = linspace(0,1,fr_count)

#     ldt_trace = detrend_linear(trial_data)
#     interp_accumulator.append( interp( linspace(0.0,1.0,1000), fake_time, ldt_trace ) )

#     figure(10)
#     plot(fake_time, detrend_linear(trial_data))

#     figure(11)
#     plot(fake_time, detrend_mean(trial_data))

#     figure(12)
#     plot(fake_time, detrend_none(trial_data))

#     figure(13)
#     plot(detrend_linear(trial_data))

#     figure(14)
#     plot(detrend_mean(trial_data))

#     figure(15)
#     plot(detrend_none(trial_data))

    


# figure(10)
# title("Linearly Detrended Mean over All Area")
# xlabel("Time in Heartbeat Space")
# ylabel("Intensity")
# savefig("lin-detrend-heartbeat-space.pdf")

# figure(11)
# title("Mean Detrended Mean over All Area")
# xlabel("Time in Heartbeat Space")
# ylabel("Intensity")
# savefig("mean-detrend-heartbeat-space.pdf")

# figure(12)
# title("Not Detrended Mean over All Area")
# xlabel("Time in Heartbeat Space")
# ylabel("Intensity")
# savefig("heartbeat-space.pdf")



# figure(13)
# title("Linearly Detrended Mean over All Area")
# xlabel("Time in Frames")
# ylabel("Intensity")
# savefig("lin-detrend-frame-space.pdf")

# figure(14)
# title("Mean Detrended Mean over All Area")
# xlabel("Time in Frames")
# ylabel("Intensity")
# savefig("mean-detrend-frame-space.pdf")

# figure(15)
# title("Not Detrended Mean over All Area")
# xlabel("Time in Frames")
# ylabel("Intensity")
# savefig("frame-space.pdf")

# show()
