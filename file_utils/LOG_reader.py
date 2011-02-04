from pyparsing import *

def str_int_float_list(token):
    # If it is a string
    token = token.strip()
    if token.startswith("'") and token.endswith("'"):
        return token

    # List
    if token.startswith("[") and token.endswith("]"):
        return [ str_int_float_list(ele) for ele in token[1:-1].split(",") ]

    # Float or int reamin
    if "." in token:
        return float(token)

    if str(int(token)) == token:
        return int(token)
    
    raise ValueError("Token (%s) is not recognised" % token)


# Prereqs for value
# pm = Literal("+") | Literal("-")
# flt_value = Optional(pm) + Word(nums) + Literal('.') + Optional(Word(nums))
# int_value = Optional(pm) + Word(nums)
# str_value = Literal("'") + SkipTo("'", include = True)
# base_vals = (flt_value  ^ int_value ^ str_value)
# list_value = Literal("[") + base_vals + ZeroOrMore(Literal(",") + base_vals) + Literal("]")
# value = flt_value  ^ int_value ^ str_value ^ list_value

# The above is clever and more accurate, but validation 
# now happens in the 'str_int_float_list' function

# Something equals something
token = Word(alphanums + "._():,")
eq_si = Literal("=")
value = SkipTo(LineEnd(), include=True) 
# TODO: This should be more specific.  
# It is ok for now (requires stripping /n etc)

token_eq_value = Combine(token) + Suppress(eq_si) + Combine(value)

# Make the header tokens have their own parse action

# These are parts of the conditions table . . .
astrix_line = Word("*")
slash_line = Word("/")

# TODO Breakup stuff inside (or ignore b/c it is all in the token data)
asterix_enclosed_data = Suppress( astrix_line + SkipTo(astrix_line, include=True) )

class openLOG:
    '''This is a class designed to deal with .log files created
    by the tempo experimental control system.  The instaciation 
    takes a filename as an argument, and creates three resulting
    dictionaries of information:
       'header_dict': Contains data from the header. 
       'footer_dict': Contains data from the footer.
       'trial_dict': Which is contiains sub-dictionaries keyed to the trial number.
    '''
    def __init__(self, filepath):
        # We define our own parts of the parser during invocation
        # This is so we can slave the appropriate parse-actions and
        # populate dictionaries with the right datae

        # Define the header token-value pair to populate the header dictionary
        header_tv = token_eq_value.copy()
        header_tv.setParseAction(self.add_header_token)
        header = OneOrMore(header_tv)

        # Same as above for a specific experiment token value
        exp_tv = token_eq_value.copy()
        exp_tv.setParseAction(self.add_experiment_token)
        experiment = OneOrMore(exp_tv)
        experiments = OneOrMore(experiment + Suppress(slash_line) ) + experiment

        # Same as above for footer token-value pairs
        foot_tv = token_eq_value.copy()
        foot_tv.setParseAction(self.add_footer_token)
        footer = OneOrMore(foot_tv)

        # Definition of how to parse a log file.
        log_file = header + asterix_enclosed_data + experiments + asterix_enclosed_data + footer

        # Make blank dictionaries that get populated by the parse
        self.header_dict = {}
        self.Current_Trial = None
        self.trial_dict = {}
        self.footer_dict = {}

        # Run the g-d parse already
        log_file.parseFile(filepath)

    def add_header_token(self, s, loc, toks):
        '''This function adds a key-value pair to a class attached dictionary 
        corresponding to data found in the log file header.'''
        key, value = toks
        self.header_dict[key] =  str_int_float_list(value)

    def add_footer_token(self, s, loc, toks):
        '''This function adds a key-value pair to a class attached dictionary 
        corresponding to data found in the log file footer.'''
        key, value = toks
        self.footer_dict[key] =  str_int_float_list(value)

    def add_experiment_token(self, s, loc, toks):
        # This needs to be a bit more complicated.  
        # If the token is the experiment number, we use that to 
        # set the currenty experiment we are parsing
        key, value = toks

        # We need to make a nested dictionary to make sure everything is kosher
        if key == "TrialNum":
            self.trial_dict[int(value)] = {}
            self.Current_Trial = self.trial_dict[int(value)]
            self.Current_Trial[key] = str_int_float_list(value)
        else:
            self.Current_Trial[key] = str_int_float_list(value)

if __name__ == "__main__":
    # Test against all the sample Frodo log files
    from pprint import pprint
    import os, sys
        
    if sys.argv[1] == "-t":
        logs_dir = "sample-logs"
        for fi in os.listdir(logs_dir):
            path = os.path.join(logs_dir, fi)
            
            print "Parsing File:", path
            logfile_reader = openLOG(path)
            sys.stdout.flush()
            print '\tDone'
    elif sys.argv[1]=="-h":
        lg = openLOG(sys.argv[2])
        pprint(lg.header_dict)

    elif sys.argv[1]=="-f":
        lg = openLOG(sys.argv[2])
        pprint(lg.footer_dict)

    elif sys.argv[1]=="-t":
        lg = openLOG(sys.argv[2])
        pprint(lg.trial_dict)
