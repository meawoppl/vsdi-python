from DIR_reader import openDIR
from DataBlocker import dir_to_h5

import os, tables, shutil

dr = openDIR("/media/jupiter/data/archive/working/frodo/")

temp_file = "/tmp/tmp-data2.h5"

b0rked = []

for path in dr.exp_paths:
    if "run" not in path: 
        print "skipping", path
        continue
    print "Doing:", path
    # Include the folder name it came from and the run number in the filename
    output_stuff = path.split("/")[-2:]
    out_fn = "-".join(output_stuff) + ".h5"

    # Assemble the full path to the aws virtual directory
    out_full = os.path.join("/home/mrg/aws-esdata/frodo", out_fn)

    # if it is already present remotley, skip it.
    if os.path.isfile(out_full): 
        try:
            print "Testing: (%s)" % out_full
            test = tables.openFile(out_full)
            test.close()
            print "Looks good, so skipping", out_fn
            continue
        except:
            print "Test failed . . ."

    # Setup compression filters, and open file for storage
    filters = tables.Filters(complevel=3, complib='lzo', shuffle=True)
    h5 = tables.openFile(temp_file, "w", filters = filters)

    # Block the data into it
    # dir_to_h5(path, h5)
    try:
        dir_to_h5(path, h5)
    except:
        b0rked.append(path)

    # Flush and close to ensure that all data is on disk
    h5.flush()
    h5.close()

    # Block on copying unitl all done
    shutil.copy(temp_file, out_full)

print b0rked

    
