from numpy import array
import os, re

# This is the regex that grabs the condition number, experiment number, and trial number
fn_re = re.compile("C([0-9]+)_?(?:[0-9]+_)E([0-9]+)B([0-9]+)")

class openDIR:
    def __init__(self, directory):
        '''Iterater over the directory and return objects that 
        help deal with experiments, trials, and logfiles.'''
        self.exp_paths  = []
        self.exp_number = []
        self.log_paths  = []
        self.total_trial_count = 0

        # Iterate over the directory structure
        for n, (dpath, dnames, fnames) in enumerate( os.walk(directory)):
            # Count the number of blk files in the directory
            number_of_blk_files = len([f for f in fnames if f.endswith('.BLK')])
            self.total_trial_count += number_of_blk_files

            # Any file w/o blk files is cruft, skip it
            if number_of_blk_files == 0:
                continue

            # Log the expiment numbers, paths, and logfile paths
            self.exp_number.append( number_of_blk_files )
            self.exp_paths.append( dpath )
            self.log_paths.append( self.get_LOG(dpath) )

    def get_LOG(self, path):
        '''Get the .log file associated with a given experimental path.'''
        # list the run folder, par is down to just .jog files, and join with the full path
        d = [ os.path.join(path, lg) for lg in os.listdir(path) if lg.endswith(".log") ]

        # Stupiud check.  One logfile only
        assert len(d) == 1, "There can be only one! (logfile): \n %s" % "\n\t".join(d)
        return d[0]

    def get_HTB(self, path):
        '''Get the .log file associated with a given experimental path.'''
        # list the run folder, par is down to just .jog files, and join with the full path
        d = [ os.path.join(path, lg) for lg in os.listdir(path) if lg.endswith(".htb") ]

        # Stupiud check.  One logfile only
        assert len(d) == 1, "There can be only one! (htbfile): \n %s" % "\n\t".join(d)
        return d[0]


    def BLK_iter(self, path):
        '''Iterate over the .BLK files associated with an experimental path.'''
        # Grab the files in the directory with ".BLK" suffix
        block_file_names = [os.path.join(path, blk) for blk in os.listdir(path) if blk.endswith(".BLK")]
        
        # Make lists to report datae with
        condition_numbers  = []
        experiment_numbers = []
        trial_numbers      = []

        # Loop over all the block files
        for filename in block_file_names:
            # Regex out the numbers we care about, and int convert them
            cn, en, tn = [int(val) for val in fn_re.findall(filename)[0]]

            # Stick them in lists
            condition_numbers.append( cn )
            experiment_numbers.append( en )
            trial_numbers.append( tn )

        # Now, lets sort this shit on the condition number
        condition_numbers = array(condition_numbers)
        experiment_numbers = array(experiment_numbers)
        trial_numbers = array(trial_numbers)
        block_file_names = array(block_file_names)

        tn_sorter = trial_numbers.argsort()
        condition_numbers = condition_numbers.take(tn_sorter)
        experiment_numbers = experiment_numbers.take(tn_sorter)
        trial_numbers = trial_numbers.take(tn_sorter)
        block_file_names = block_file_names.take(tn_sorter)

        # That was fun.
        return zip(condition_numbers, experiment_numbers, trial_numbers, block_file_names)

    def EXP_iter(self):
        '''Iterate over the experiment files assoaciated in the dirctory analyzed.'''
        return iter( self.exp_paths )

if __name__ == "__main__":
    # Run a quick test of things. . . 
    openDIR("/media/jupiter/data/archive/working/frodo")
