import os, numpy, pickle, hashlib

class DoOnce:
    def __init__(self, dec_f):
        self.f = dec_f
        print "In Init"

    def __call__(self, *args, **kwargs):
        try:
            doover = kwargs.pop("doover")
        except KeyError:
            doover = False

        md5 = hashlib.md5()
            
        # Figure out what we would save it as . . .
        # make a hash on the args and kwargs
        arg_strings  = [str(a) for a in args] 
        arg_strings += [ str(k) + "=" + str(v) for k, v in kwargs.iteritems() ]
        md5.update( ".".join(arg_strings) )
        cache_string = self.f.func_name + "-" + md5.hexdigest()  
        
        # the +4 is for the ".ext"
        if len(cache_string) + 4 >= 254:
            raise RuntimeError("Cache String Too Long")

        # Different file names for numpy arrays and pkls
        pkl_file = os.path.join(".doonce-cache/", cache_string + ".pkl")
        npy_file = os.path.join(".doonce-cache/", cache_string + ".npy")

        # Look for an already computed the result
        # We dont know if it will be a npy or pkl so check both
        if os.path.isfile(pkl_file) and doover == False:
            print "Previous PKL found:", pkl_file
            return pickle.load(open(pkl_file))
        if os.path.isfile(npy_file) and doover == False:
            print "Previous NPY found:", npy_file
            return numpy.load(npy_file)
        print "No previous result found, or redoing. . ."

        # I guess we didnt compute the result, or we are recomputing it . . .
        print "Computing Result . . ."
        res = self.f(*args, **kwargs)
        
        # If it is an ndarray, save it that way, otherwise . . . pickle it!
        if isinstance(res, numpy.ndarray):
            print "Saving numpy array", npy_file
            flo = open(npy_file, "w")
            numpy.save(flo, res)
        else:
            print "Saving pickle array", pkl_file
            flo = open(pkl_file, "w")
            pickle.dump(res, flo)
        # Flush and close the file object
        flo.flush()
        flo.close()

        return res
