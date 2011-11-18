import os, sys, time
from optparse import OptionParser

from seqtools import io

class Command(object):
    """Abstracts boiler plate for command line scripts.
    
    Includes a handy abstraction to deal with input/output files/streams.
      + If input argument is expected, set `in_arg` to its position in the
        argument list. If this argument is missing or is "-", then stdin is read.
      + Ditto for output/stdout
    
    Command line options added automatically:
        * -v/--verbose : control verbosity
        * -f/--force : overwrite output if already exists
    """
    
    def __init__(self, name, parser, in_arg=0, out_arg=1):
        """Create Command Line program
        
        set in_arg,out_arg to None if you do not want to process an
        input/output file
        
        """
        self.name = name
        self.to_stdout = False
        self.from_stdin = False
        self.infile = None
        self.outfile = None
        self.overwrite = False
        self.options = None
        self.args = None
        self.stats = []
        self.__process_parser(parser)
        self.__init_file_streams(in_arg, out_arg)
    
    def __process_parser(self, parser):
        parser.add_option('-f', '--force', dest="force", default=False,
                          action="store_true", help="Ignore any precautions")
        parser.add_option("-v", "--verbose", dest="verbose", default=False,
                          action="store_true", help="Make some noise")
        (self.options, self.args) = parser.parse_args()
        self.parser = parser
        
        
    def __init_file_streams(self, in_arg=0, out_arg=1):
        if in_arg is not None:
            if len(self.args) - 1 < in_arg or self.args[in_arg] == '-':
                self.infile = sys.stdin
                self.from_stdin = True
            else:
                self.infile = self.args[in_arg]
                if not os.path.isfile(self.infile):
                    parser.error("Cannot read input file: " + self.infile)
                self.infile = io.xopen(self.infile, 'r')
                self.from_stdin = False

        if out_arg is not None:
            if len(self.args) - 1 < out_arg or self.args[out_arg] == '-':
                self.outfile = sys.stdout
                self.to_stdout = True
            else:
                self.outfile = self.args[out_arg]
                if os.path.isfile(self.outfile) and not self.options.force:
                    self.error("Outfile already exists, " \
                               "use -f/--force to overwite")
                self.outfile = io.xopen(self.outfile, 'w')
                self.to_stdout = False

    def error(self, msg, print_help=False):
        """Produce an error and exit"""
        self.parser.error(msg)
    
    def warning(self, msg):
        if self.to_stdout:
            sys.stderr.write(msg + "\n")
        else:
            sys.stdout.write(msg + "\n")

    def add_stat(self, name, value):
        self.stats.append((name, value))

    def run(self, what):
        t0 = time.time()
        what(self)
        elapsed = time.time() - t0

        if self.infile is not None and not self.from_stdin:
            self.infile.close()
        if self.outfile is not None and not self.to_stdout:
            self.outfile.close()
        
        report = sys.stderr if self.to_stdout else sys.stdout
        
        report.write('\n')
        report.write('========== ' + self.name + ' Finished ==========\n')
        report.write('Elapsed time:  %.2f seconds\n' % elapsed)
        for (name,val) in self.stats:
            report.write('%s: %s\n' % (name, val))
        report.write("\n")
    
# END : Class Command
