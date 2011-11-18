import sys
from optparse import OptionParser

from seqtools import io

class Command(object):
    """Abstracts boiler plate for producing an output from one input.
    
    If input (first argument) is missing or "-" reads stdin
    If output (seond argument) is misgin or "-" writes to sdout
    
    Takes care of checking input files / streams etc.
    
    Command line options added automatically:
        * -v/--verbose : control verbosity
        * -f/--force : overwrite output if already exists
    """
    
    def __init__(self, parser):
        self.to_stdout = False
        self.from_stdin = False
        self.infile = None
        self.outfile = None
        self.overwrite = False
        self.options = None
        self.args = None
        self.stats = {}
        self.__process_parser(parser)
        self.__init_file_streams()
    
    def __process_parser(self, parser):
        parser.add_optoin('-f', '--force', dest="force", default=False,
                          action="store_true", help="Ignore any precautions")
        parser.add_option("-v", "--verbose", dest="verbose", default=False,
                          action="store_true", help="Make some noise")
        (self.options, self.args) = parser.parse_args()
        self.parser = parser
        
        
    def __init_file_streams(self):
        if len(self.args) == 0:
            self.infile = sys.stdin
            self.from_stdin = True
            self.outfile = sys.stdout
            self.to_stdout = True
        if len(self.args) > 0:
            if self.args[0] == '-':
                self.infile = sys.stdin
                self.from_stdin = True
            else:
                self.infile = self.args[0]
                if not os.path.isfile(infile):
                    parser.error("Cannot read input file: " + self.infile)
                self.infile = io.xopen(self.infile, 'r')
                self.from_stdin = False
            self.outfile = sys.stdout
            self.to_stdout = True
        if len(self.args) > 1:
            self.outfile = self.args[2]
            if os.path.isfile(self.outfile) and not self.overwrite:
                parser.error("Destination file already exists")
            outfile = io.xopen(outfile, 'w')
            close_out = True
        

    def error(self, msg, print_help=True):
        """Produce an error and exit"""
        self.parser.error(msg)
    
    def warning(self, msg):
        if self.to_stdout:
            sys.stderr.write(msg + "\n")
        else:
            sys.stdout.write(msg + "\n")
    
    

# END : Class Command
