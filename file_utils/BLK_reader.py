# Functions to work with BLK imaging file
# By MRG Based on the the matlab code of YC
# 
# '''
# Original header comments:
# www.opt-imaging-inc.com/download/adobe/VdaqTasks3001.pdf
# Appendix A in VDaqTasks3001.pdf
# fnBLK is the BLK filename incliding pathname.
#
#
# YC at ES lab
# Created on Apr. 13, 2008
# Last modified on Apr. 14, 2008
# '''
# Created Aug. 19, 2010

import struct, numpy
from pylab import *

decode_data_type = {11:numpy.char,        # char           - 1 byte
                    12:numpy.ushort,      # unsigned short - 2 byte
                    13:numpy.long,        # long           - 4 byte 
                    14:numpy.float32}     # float          - 8 byte
# Blk file format . . .
headerformat = ([("l", "FileSize"),
                 ("l", "CheckSum_Header"),
                 ("l", "CheckSum_Data"),
                 ("l", "LenHeader"),
                 ("l", "VersionID"),
                 ("l", "FileType"),    
                 ("l", "FileSubtype"), 
                 ("l", "DataType"),    
                 ("l", "SizeOf"),	     
                 ("l", "FrameWidth"),
                 ("l", "FrameHeight"),
                 ("l", "NFramesPerStim"),
                 ("l", "NStimuli"),
                 ("l", "InitialXBinFactor"),
                 ("l", "InitialYBinFactor"),
                 ("l", "XBinFactor"),
                 ("l", "YBinFactor"),
                 ("32s", "acUserName"),
                 ("16s", "acRecordingDate"),
                 ("l", "X1ROI"),
                 ("l", "Y1ROI"),
                 ("l", "X2ROI"),
                 ("l", "Y2ROI"),                           
                 ("l", "StimOffs"),
                 ("l", "StimSize"),
                 ("l", "FrameOffs"),
                 ("l", "FrameSize"),
                 ("l", "RefOffs"),
                 ("l", "RefSize"),
                 ("l", "RefWidth"),
                 ("l", "RefHeight"),
                 # These are both actully 16 unsigned short
                 # This is length Equiv, but gloms them together
                 ("32s", "aushWhichBlocks"),
                 ("32s", "aushWhichFrames"),

                 ("f", "LoClip"),
                 ("f", "HiClip"),
                 ("l", "LoPass"),
                 ("l", "HiPass"),
                 ("64s", "acOperationsPerformed"),
                 ("f", "Magnification"),
                 ("H", "ushGain"),
                 ("H", "ushWavelength"),
                 ("l", "ExposureTime"),
                 ("l", "NRepetitions"), 
                 ("l", "AcquisitionDelay"), 
                 ("l", "InterStimInterval"),
                 ("16s", "acCreationDate"),
                 ("64s", "acDataFilename"),
                 ("256s", "acOraReserved"),
                 ("l", "IncludesRefFrame"), 
                 ("265s", "acListOfStimuli"),
                 ("l", "NFramesPerDataFrame"),
                 ("l", "NTrials"),
                 ("l", "ScaleFactor"),
                 ("f", "MeanAmpGain"),
                 ("f", "MeanAmpDC"),
                 ("B", "ucBegBaselineFrameNo"), 
                 ("B", "ucEndBaselineFrameNo"),
                 ("B", "ucBegActivityFrameNo"),
                 ("B", "ucEndActivityFrameNo"),
                 ("B", "ucDigitizerBits"), 
                 ("B", "ucActiveSystemID"),
                 ("B", "ucDummy2"),
                 ("B", "ucDummy3"),
                 ("l", "X1SuperPix"),
                 ("l", "Y1SuperPix"),
                 ("l", "X2SuperPix"),
                 ("l", "Y2SuperPix"),
                 ("f", "FrameDuration"),
                 ("l", "ValidFrames"),
                 ("224s", "acVdaqReserved"),
                 ("16s", "rTimeBlockStart"),
                 ("16s", "rTimeBlockEnd"),
                 ("224s", "acUser"),
                 ("256s", "acComment") ])

# Assemble the above tuples into a reading format and list of varible names
format_string = ""
var_names = []
for var_format, var_name in headerformat:
    format_string += var_format
    var_names += [var_name]
    
# Compute the length of the header based on string
header_size_in_bytes = struct.calcsize(format_string)

class openBLK:
    '''This is a class for reading BLK files as generated by
    TODO
    '''
    def __init__(self, filepath):
        self.filename = filepath
        self.header_loaded = False
        self.data_loaded = False
        # Dont call load header or data until that time 
        # b/c we may only want certain information

    def force_header(self):
        if self.header_loaded:
            return
        self.load_header()
        self.header_loaded = True

    def force_data(self):
        if self.data_loaded:
            return
        self.load_data()
        self.data_loaded = True

    def load_header(self):
        self.header = {}
        # Open the block file
        block_data_file = open(self.filename, "rb")
        # Read the header and unpack based on format string
        header_data = block_data_file.read(header_size_in_bytes)
        data_from_header = struct.unpack(format_string, header_data)

        # Assign the header stuff into a dictionary.
        # Some of these are used in loading image arrays
        for name, value in zip(var_names, data_from_header):
            self.header[name] = value        

        # Assertions about file structure i.e. sanity checks
        # print data_from_header[3] + 12
        # print header_size_in_bytes
        assert (header_size_in_bytes == (data_from_header[3] + 12)), "File header size mismatch"
        assert len(data_from_header) == len(var_names), "Incorrect variable count unpacked"
        
    def load_data(self):
        # Make sure the header is loaded up
        self.force_header()
        
        # Open datafile and skip to the offset of the data start
        block_data_file = open(self.filename, "rb")
        block_data_file.seek(self.header["LenHeader"])

        # Decode the data dtype (u)int8/16 etc.
        self.dtype = decode_data_type[self.header["DataType"]]

        # Load the raw bin and intrepate
        self.data_array = numpy.fromstring(block_data_file.read(), dtype=self.dtype)
        self.data_array = self.data_array.reshape((-1, 
                                                    self.header["FrameWidth"], 
                                                    self.header["FrameHeight"]))

        # Switch axes so image is shape (x, y, t)
        self.data_array = self.data_array.swapaxes(0,2)
    
    

    def print_header(self):
        # Make sure the header is loaded up
        self.force_header()
        # Print the vars in the header order
        for v in var_names:
            print v, "=", self.header[v]
        

def numpy_array_to_pgm(filename, nda, mx_value=None):
    # This is a terrible format.  That said this implementaiton is
    # based on the specification given at:
    # http://netpbm.sourceforge.net/doc/pgm.html
    # Check to see if it can be stored w/o loss in this format
    assert nda.ndim  == 2
    assert (nda.dtype == uint8) or (nda.dtype == uint16)
    
    # Maxes for uint 8/16
    if (nda.dtype == uint8 ) and (mx_value==None): mx_value = 255
    if (nda.dtype == uint16) and (mx_value==None): mx_value = 65535
    
    # Figure out values for header
    w, h = nda.shape

    # If we have gotten this far without error
    # Open the file 
    f = open(filename, "w")

    # Write the header
    header = "P5\n%i\n%i\n%i\n" % (w, h, mx_value)
    print header
    f.write(header)
    f.write(nda.tostring())

# xs, ys = mgrid[0:256, 0:256]
# xs = xs.astype(uint8)
# ys = ys.astype(uint8)
# numpy_array_to_pgm("test.pgm", ys)

# 1/0

if __name__ == "__main__":
    '''This is a not-bad unittest'''
    # TODO: nose unittesting
    bf = BlockFileReader("test-blockfile.BLK")
    bf.load_data()
    bf.print_header()
    
    vals = []
    max_brightness_of_all_frames = bf.data_array.max()

    for frame_number in range(bf.data_array.shape[2]):
        vals.append(bf.data_array[:,:,frame_number].mean())
        
        print bf.data_array[:,:,frame_number].dtype
        
        print bf.data_array.max(), bf.data_array.min()

        numpy_array_to_pgm("frames/%04i.pgm" % frame_number, 
                           bf.data_array[:,:,frame_number],
                           mx_value = max_brightness_of_all_frames)
    



