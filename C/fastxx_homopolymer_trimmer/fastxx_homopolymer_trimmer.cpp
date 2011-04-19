// This tool depends on the fastx toolkit.
// Make sure libfastx and gtextutils is in your -I(include) path when compiling.
// 
// Use it to strip stretches of homopolymers at the end of a sequence.
// 
// To install:
//   * Ensure that gtestutils and fastx is installed
//     libgtextutils.* will be installed into /usr/local/lib
//     
//   + Create /usr/local/include/fastx
//   + Copy /path/to/fastx/src/*.h to /usr/local/include/fastx/
//   + Copy /path/to/fastx/src/libfastx/libfastx.* in /usr/local/lib
#include <err.h>
#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include "fastx/fastx.h"
#include "fastx/fastx_args.h"

#include "sequtils/strings.h"

#define MAX_TRASHFILE_LEN 1000

const char* usage=
"usage: fastxx_homopolymer_trimmer [-h] [-v] [-b A] [-l N] [-t trashfile] " \
"[-i INFILE] [-o OUTFILE]\n" \
"by Steve Lianoglou (slianoglou@gmail.com)\n" \
"\n" \
"   [-h]         = Help me help you.\n" \
"   [-v]         = verbose: print short summary of input/output counts\n"      \
"   [-b N]       = The basepair-run to strip from 3' end of read [A]\n"        \
"   [-t outfile] = The path to the file to save 'trashed' reads in. If the\n"  \
"                  filename ends in *.gz, the file will be compressed.\n"      \
"   [-l N]       = discard sequences shorter than N nucleotids [5]\n"          \
"   [-i INFILE]  = FASTA/Tabular input file. default is STDIN.\n"              \
"   [-o OUTFILE] = FASTA/Tabular output file. default is STDOUT.\n"            \
"\n";


// Working parameters
char NT = 'A';
unsigned int min_length = 5;
int save_trash = 0;
char trashfile[MAX_TRASHFILE_LEN] = "homopolymer.trash.fastq";
int compress_trash = 0;

FASTX fastx;
FASTX trashx;

int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg) {
  switch(optc) {
    case 'b':
      if (optarg == NULL) {
        errx(1, "[-b] parameter requires an argument value");
      }
      if (strlen(optarg) != 1) {
        errx(1, "[-b] can only be a character of length 1, eg: A, C, G, or T");
      }
      NT = optarg[0];
      break;
    
    case 'l':
      if (optarg == NULL) {
        errx(1, "[-l] parameter requires a value");
      }
      min_length = strtoul(optarg, NULL, 10);
      if (min_length < 5) {
        errx(1, "Minimum length of reads cannot be less than 5.");
      }
      break;
    
    case 't':
      if (optarg == NULL) {
        errx(1, "[-t] parameter needs name of file to save trash to.");
      }
      // TODO: Check trashfile is valid filename AND it doesn't already exist
      if (strlen(optarg) > MAX_TRASHFILE_LEN) {
        errx(1, "[-t] trashfile name is too long");
      }
      
      save_trash = 1;
      strncpy(trashfile, optarg, sizeof(trashfile) - 1);
      
      // Compress it?
      if (strncmp(".gz", substring(trashfile, strlen(trashfile) - 3, 3), 3) == 0) {
        fprintf(get_report_file(), "Compressing: %s\n", trashfile);
        compress_trash = 1;
      }
      break;
    
    default:
      errx(1,"Internal error: unknown option '%c'",optc);
  }
  return 1;
}

unsigned int homopolymer_start(const std::string& read, char nt);
unsigned int homopolymer_start(const std::string& read, char nt) {
  unsigned int run_start = read.length();
  unsigned int i;
  
  for (i = read.length(); i > 0;) {
    if (read[i - 1] != nt) {
      break;
    }
    run_start = i--;
  }
  if (i < read.length()) {
    run_start--;
  }
  return run_start;
}

int main(int argc, char* argv[]) {
  int skip_read;
  unsigned int axe_at;
  unsigned int reads_count = 0;
  unsigned int trashed_count = 0;
  unsigned int trimmed_count = 0;
  unsigned int bases_trimmed = 0;
  std::string read;
  
  // parse_commandline(argc, argv);
  fastx_parse_cmdline(argc, argv, "b:l:t:", parse_program_args);
  
  fastx_init_reader(&fastx, get_input_filename(), FASTA_OR_FASTQ, ALLOW_N,
                    REQUIRE_UPPERCASE, get_fastq_ascii_quality_offset());
  fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT,
                    compress_output_flag());
  
  // Setup the trash?
  if (save_trash) {
    fastx_init_reader(&trashx, get_input_filename(), FASTA_OR_FASTQ, ALLOW_N,
      REQUIRE_UPPERCASE, get_fastq_ascii_quality_offset());
    fastx_init_writer(&trashx, trashfile, OUTPUT_SAME_AS_INPUT,
      compress_trash);
  }
  
  while (fastx_read_next_record(&fastx)) {
    reads_count++;
    skip_read = 0;
    read = std::string(fastx.nucleotides);
    axe_at = homopolymer_start(read, NT);
    
    if (axe_at < min_length) {
      skip_read = 1;
    } else if (axe_at < read.length()) {
      // Trim the read at this position
      fastx.nucleotides[axe_at] = 0;
      trimmed_count++;
      bases_trimmed += read.length() - axe_at;
    }
    
    if (skip_read) {
      trashed_count++;
      if (save_trash) {
        fastx_read_next_record(&trashx);
        fastx_write_record(&trashx);
      }
    } else {
      fastx_write_record(&fastx);
    }
  }
  
  // Print verbose report
  if (verbose_flag()) {
    fprintf(get_report_file(), "Number of reads processed: %d\n", reads_count);
    fprintf(get_report_file(), "Number of reads trashed: %d\n", trashed_count);
    fprintf(get_report_file(), "Number of reads trimmed: %d\n", trimmed_count);
    fprintf(get_report_file(), "Number of bases trimmed: %d\n", bases_trimmed);
  }
  
  return 0;
}
