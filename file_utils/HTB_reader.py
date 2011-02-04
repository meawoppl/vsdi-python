import struct
from numpy import ( int8, int16, uint8, uint16, 
                    nbytes, fromfile )

# nCols = htb.nchannels
# nRows = htb.period

# So the value labled "func" actully describes the data type of the array in the DB
# The following dictionaries decodes that into python types, I also included the original
# Comments from the c header file . . . 
db_format_decode = {0:int16,         # signed 8 bit analog SUM  
                    1:int8,          # signed 8 bit analog APP  
                    2:uint16,        # unsigned counter USUM    
                    3:uint16,        # unsigned counter UAPP    
                    4:uint16,        # unsigned event ESUM      
                    5:uint16,        # unsigned event EAPP      
                    6:int16,         # signed 12 bit analog XSUM
                    7:int16 }        # signed 12 bit analog XAPP

db_sumappend_decode = {0:"sum", 2:"sum", 4:"sum", 6:"sum",
                       1:"app", 3:"app", 5:"app", 7:"app"}

db_type_decode = {0:"analog",  1:"analog", 
                  2:"counter", 3:"counter",
                  4:"event",   5:"event", 
                  6:"event",   7:"event" }

headerformat = ([("26s", "date"),
                 ("l", "ldate"), 
                 ("14s", "cfg_file"),
                 ("14s", "pro_file"),
                 ("52s", "UNUSED"),
                 ("L", "speed"),
                 ("L", "alloc"),
                 ("l", "offset"),
                 ("L", "period"),
                 ("L", "extension"),
                 ("H", "skip"),
                 ("H", "first_channel"),
                 ("H", "nchannels"),
                 ("H", "sweep_limit"),
                 ("L", "cancel_override"),
                 ("B", "func"),
                 ("s", "UNUSED2"),
                 ("H", "tag"),
                 ("H", "npages"),
                 ("L", "nsamples"),
                 ("H", "samples_per_page"),
                 ("H", "sweep"),
                 ("H", "next_page"),
                 ("H", "next_off"),
                 ("80s", "title"),
                 ("L", "speed_units"),
                 ("268s", "UNUSED3")])

# Assemble the above tuples into a reading format and list of varible names
# By sepcifying the endinaness, you force it to go w/o local compiler specified padding
format_string = "<"
var_names = []

for var_format, var_name in headerformat:
    format_string += var_format
    var_names += [var_name]

# Stupiud check on header size.  
header_size_in_bytes = struct.calcsize(format_string)
assert header_size_in_bytes == 512, "Internal sanity check failed!  Header should be 512 bytes per specification."

class openHTB:
    '''This is a class designed to read the .htb file format created by
    Reflective Computing (http://www.reflectivecomputing.com/) TEMPO
    Experiment Coordination Platform.
    
    The reader reads everything right off the bat.  That said if a .htb file
    is larger than the memory of the machine, this may need to be re-tooled.  
    Given that there are structures dedicated to EMM, I don't 
    think this will be an issue.  

    This code uses numpy for array storage, though any number of array like objects
    could be dropped in, including numeric, pyrex, etc.

    Signature: 
    htb = openHTB(filename)
    
    After loading the object will have three lists of interest, all keyed to the sequential
    databases encountered:

    htb.db_names contains the string name of the database represented
    htb.db_metas is a list ordered by the order of the appearing database in the .htb file
        Each entry contains a dictionary keyed on the names for varibles in the C-spec'ed
        header format that appears earlier in this file.
    htb.db_array is a list of numpy arrays ordered by appearance of the database in the 
    .htb file
    '''

    def __init__(self, fn, debug=False):
        # Debugging flag
        self.debug = debug
        # Save the filename and open the file
        self.filename = fn
        self.flo = open(self.filename, "rb")
        
        # Create the structures we will populate
        self.db_names = []
        self.db_frmts = []
        self.db_types = []
        self.db_smapp = []
        self.db_metas = []
        self.db_array = []
        
        # Populate all of them
        self._load_all()

    def dbprint(self, *args):
        if self.debug: print args

    def _load_header(self, offset=0):
        '''Loads the 512 byte header given at a certain offset in an htb file'''
        self.dbprint("Loading header of file (%s) at offset (%i)" % (self.filename, offset))

        header = {}
        self.flo.seek(offset)
        # Read the header and unpack based on format string
        header_data = self.flo.read(header_size_in_bytes)
    
        # If you seeked to the file-end looking for the next db:
        if len(header_data) == 0:
            return None
        elif len(header_data) < 512:
            raise RuntimeError("Unexpected end of htb file!")
        data_from_header = struct.unpack(format_string, header_data)

        # Assign the header stuff into a dictionary.
        # Some of these are used in loading image arrays
        for name, value in zip(var_names, data_from_header):
            # It it is a string remove extraneous null bytes
            if type(value) == str: value = value.replace("\x00", "")
            self.dbprint(name, "==>", value)
            header[name] = value        
        return header 
    
    def _load_array(self, data_offset, data_len, nchan, format):
        # Seek to the start of the actual data, skipping over the header
        
        self.flo.seek(data_offset)
        # Read the record the header length
        data_points = data_len / nbytes[format]
        aray_data = fromfile(self.flo, dtype=format, count = data_points)
        # Arrange it into a (time x channels) array
        self.db_array += [aray_data.reshape((-1, nchan)).T.squeeze()]

    def _load_all(self):
        '''This function loads the htb header at the 0th offset, 
        it then uses it to find the next one, repeat until file end'''
        current_offset = 0
        while True:
            self.dbprint( "*" * 25)
            # Try to load the next header.  If the file is empty or
            # has reached the end, this routine is done!
            meta_data = self._load_header(current_offset)
            if meta_data == None:
                break

            # Sanitize the name and add it to our list-o-names
            self.db_names += [meta_data["title"].strip("\x00")]

            # Find the length of the DB, and channel count
            db_length   = meta_data["alloc"]
            chan_count  = meta_data["nchannels"]
            data_points = meta_data["period"]

            # Compute the length of the actual data array in bytes, and its offset
            data_length = db_length - 512
            data_offset = current_offset + 512
            
            # This should decode the "func" value into a numpy data type
            # And database type, and whether is is summing or appending?
            frmt =    db_format_decode[ meta_data["func"] ]
            typ  =      db_type_decode[ meta_data["func"] ]
            sora = db_sumappend_decode[ meta_data["func"] ]

            # Populate the list with the various whatnot
            self.db_frmts += [ frmt ]
            self.db_types += [ typ ]
            self.db_smapp += [ sora ]
            self.db_metas += [ meta_data ]
            
            # Stupid checks . . . 
            assert (data_length % chan_count) == 0, "There should be an even number of data points for each channel"

            # Use the above datae to get the actual array data . . .
            self._load_array(data_offset, data_length, chan_count, frmt)
            current_offset += meta_data["alloc"] 


if __name__ == "__main__":
    import sys
    res = openHTB(sys.argv[1], debug=True)
