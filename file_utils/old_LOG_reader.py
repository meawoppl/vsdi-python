from pyparsing import *

def str_int_float_list(token):
    token = token.strip()
    return_token = None
    for fun in [str, float, int]:
        try:
            if str(fun(token)) in token:
                return_token = fun(token)
        except ValueError:
            pass
    if return_token == None: 
        raise ValueError("Token (%s) is not recognised" % token)

    return return_token

# The above is clever and more accurate, but validation 
# now happens in the 'str_int_float_list' function

# Something equals something
token = Word(alphanums + "._():,-[]")
eq_si = Literal("=")
wht   = White(ws=" ")
value = Word(alphanums + "._():,-[]")
# TODO: This should be more specific.  
# It is ok for now (requires stripping /n etc)
token_eq_value = Combine(token) + Suppress(wht) + Combine(SkipTo(LineEnd(), include=True))

# Make the header tokens have their own parse action
# These are parts of the conditions table . . .
astrix_line = Word("*") + SkipTo( LineEnd(), include=True)
slash_line = Word("/")  + SkipTo( LineEnd(), include=True)

# TODO Breakup stuff inside (or ignore b/c it is all in the token data)
asterix_enclosed_data = Suppress( astrix_line + SkipTo(astrix_line, include=True) )

foot_start = astrix_line + Literal("Begin of File Footer") + astrix_line
foot_end   = astrix_line + Literal("End of File Footer") + astrix_line

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
        self.header_tv = token_eq_value.copy()
        self.header_tv.setParseAction(self.add_header_token)
        self.header = OneOrMore(self.header_tv)

        # Same as above for a specific experiment token value
        self.exp_tv = token_eq_value.copy()
        self.exp_tv.setParseAction(self.add_experiment_token)
        self.experiment = self.exp_tv + OneOrMore(self.exp_tv)
        self.experiments = delimitedList(self.experiment, delim=astrix_line)

        # Same as above for footer token-value pairs
        self.foot_tv = token_eq_value.copy()
        self.foot_tv.setParseAction(self.add_footer_token)
        self.footer = OneOrMore(self.foot_tv)

        self.old_header = (self.header + asterix_enclosed_data + self.header + 
                      asterix_enclosed_data + asterix_enclosed_data + astrix_line)
        self.old_footer = Suppress( foot_start + SkipTo(foot_end, include=True) )
        # self.old_footer = (asterix_enclosed_data + asterix_enclosed_data + self.footer + asterix_enclosed_data)
        # Definition of how to parse a log file.
        self.log_file = self.old_header + self.experiments + self.old_footer

        # Make blank dictionaries that get populated by the parse
        self.header_dict = {}
        self.Current_Trial = None
        self.trial_dict = {}
        self.footer_dict = {}

        # Run the g-d parse already
        self.log_file.parseFile(filepath)

    def add_header_token(self, s, loc, toks):
        '''This function adds a key-value pair to a class attached dictionary 
        corresponding to data found in the log file header.'''
        key, value = toks
        self.header_dict[key] =  str_int_float_list(value)
        print "Header:", key, ":", value

    def add_footer_token(self, s, loc, toks):
        '''This function adds a key-value pair to a class attached dictionary 
        corresponding to data found in the log file footer.'''
        key, value = toks
        self.footer_dict[key] =  str_int_float_list(value)
        print "Footer:", key, ":", value

    def add_experiment_token(self, s, loc, toks):
        # This needs to be a bit more complicated.  
        # If the token is the experiment number, we use that to 
        # set the currenty experiment we are parsing
        key, value = toks
        print "Exp:", key, ":", value
        # We need to make a nested dictionary to make sure everything is kosher
        if (key == "TrialNum") or (key == "CountBlockTotal"):
            self.trial_dict[int(value)] = {}
            self.Current_Trial = self.trial_dict[int(value)]
            self.Current_Trial[key] = str_int_float_list(value)
        else:
            self.Current_Trial[key] = str_int_float_list(value)

if __name__ == "__main__":
    # Test against all the sample Frodo log files
    import os, sys

    path = "old-log-ex.log"
    logfile_reader = openLOG(path)
