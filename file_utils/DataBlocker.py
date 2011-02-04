from BLK_reader import openBLK
from HTB_reader import openHTB
from LOG_reader import openLOG
from DIR_reader import openDIR

import tables, os
from tables.nodes import filenode

from numpy import *

def dir_to_h5(directory, h5_obj, h5_path="/"):
    # Use the openDIR object to deal with the ugliness
    dir_listing = openDIR( directory )

    # You know what they say about assertions
    assert (len(dir_listing.exp_paths) == 1), "This function only does 1 run at a time:" + str(dir_listing.exp_paths)
    print dir_listing.total_trial_count, "blockfiles detected, compressing . . ."

    # Create a leaf in the h5 for HTB and LOG file related data
    htb_leaf = h5_obj.createGroup(h5_path, "htb")
    log_leaf = h5_obj.createGroup(h5_path, "log")

    # Create filenodes
    htb_fn = filenode.newNode(h5_obj, where=htb_leaf, name='htbfile')
    log_fn = filenode.newNode(h5_obj, where=log_leaf, name='logfile')

    # Get file paths
    log_file = dir_listing.get_LOG(dir_listing.exp_paths[0])
    htb_file = dir_listing.get_HTB(dir_listing.exp_paths[0])

    # Store the raw files (they are pretty small, so why not)
    htb_fn.write(open(htb_file).read())
    log_fn.write(open(log_file).read())

    # Close to force flush
    htb_fn.close()
    log_fn.close()

    # Open the HTB file, and break out and write each array
    htb_file = openHTB(htb_file)
    for a_name, htb_arr in zip( htb_file.db_names, htb_file.db_array ):
        # Figure out the type in terms of tables symantics
        atm = tables.Atom.from_dtype(htb_arr.dtype)
    
        # Takes spaces out of the name, so we can use natural naming
        nat_name = a_name.replace(" ","")

        # Create the array and throw the data in it
        leaf_arr = h5_obj.createCArray(htb_leaf, nat_name, 
                                      atm, htb_arr.shape, 
                                      a_name)
        leaf_arr[:] = htb_arr

    # Block file bits
    # Conv. lists
    condition = []
    trial_num = []

    # Pack in the block files
    array_created = False
    for n, (cn, bn, tn, f) in enumerate(dir_listing.BLK_iter(dir_listing.exp_paths[0])):
        print "Starting", n, "of", dir_listing.total_trial_count
        print "\tLoading .BLK (%s) cn:%i tn:%i" % (f, cn, tn)

        # Open the block file
        bf = openBLK(f)
        bf.load_data()

        # Create a c-array
        if not array_created:
            new_shape = bf.data_array.shape + tuple([0])
            ear = h5_obj.createEArray("/", "block", tables.Int32Atom(), new_shape, "Data from '%s'" % directory)
            array_created = True


        print "\tWriting to h5."
        # Make the stings for the name and description of each data bloc
        ear.append( bf.data_array[:,:,:,newaxis] )

        # Mem clear necessary?  TODO
        ear.flush()
        h5_obj.flush()
        del bf

        # For later conv.
        condition.append(cn)
        trial_num.append(tn)
        
        print "\tDone."

    # Create the condition and trial number arrays.  
    # (i.e avoid redoing that ugly file name munging each time)
    h5_obj.createArray(h5_path, "condition", array(condition))
    h5_obj.createArray(h5_path, "trial_num", array(trial_num))


    # Yay done!
    print "Done!"
    # Other file closing is dealt with through object destruction later


if __name__ == "__main__":
    import sys
    usage = '''python DataBlocker.py [dir-to-block] [h5-out]'''

    assert len(sys.argv) == 3, usage

    # Setup the compression, and open an h5 file as such
    filters = tables.Filters(complevel=3, complib='zlib', shuffle=True)
    h5 = tables.openFile(sys.argv[2], "w", filters = filters)
    
    # Push the dir into it
    dir_to_h5(sys.argv[1], h5)

    # This should happen anyway, but why not
    h5.flush()
    h5.close()
