/*
  This tool depends on the fastx toolkit.
  Make sure libfastx and gtextutils is in your -I(include) path when compiling.
  
  Use it to strip stretches of homopolymers at the end of a sequence.
*/
#include <err.h>
#include <getopt.h>
#include <string.h>
#include <algorithm>
#include <cstdlib>
#include <ios>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <ostream>
#include <fstream>

#include "config.h"

#include "fastx.h"
#include "fastx_args.h"

using namespace std;

const char* usage=
"usage: fastxx_homopolymer_trimmer [-h] [-v] [-b A] [-l N] [-i INFILE] [-o OUTFILE]\n" \
"by Steve Lianoglou (slianoglou@gmail.com)\n" \
"\n" \
"   [-h]         = Help me help you.\n" \
"   [-v]         = verbose: print short summary of input/output counts\n" \
"   [-b N]       = The basepair to hunt for and remove from 3' end of read\n" \
"                  Defaults to A\n" \
"   [-l N]       = discard sequences shorter than N nucleotids (default is 5)\n" \
"   [-i INFILE]  = FASTA/Tabular input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Tabular output file. default is STDOUT.\n" \
"\n";

// Working parameters
char NT = 'A';
int min_length = 5;
int verbose = 0;

FASTX fastx;

int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg)
{
  switch(optc) {
    case 'b':
      if (optarg == NULL) {
        errx(1, "[-b] parameter requires an argument value");
      }
      if (strlen(optarg) != 1) {
        errx(1, "[-b] can only be a character of length 1 (A, C, G, or T)");
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
    
    // case 'v':
    //   verbose = 1;
    //   break;
    
    default:
      errx(1,"Internal error: unknown option '%c'",optc);
  }
  return 1;
}

// int parse_commandline(int argc, char* argv[]) {
//   fastx_parse_cmdline(argc, argv, "b:l:");
//   return 1;
// }

int homopolymer_start(const std::string& read, char nt);
int homopolymer_start(const std::string& read, char nt) {
  int run_start = -1;
  int i;
  
  for (i = read.length(); i >= 0;) {
    if (read[i - 1] != nt) {
      break;
    }
    run_start = i--;
  }
  return run_start;
}

int main(int argc, char* argv[]) {
  int axe_at, skip_read;
  unsigned int reads_count = 0;
  unsigned int trashed_count = 0;
  unsigned int trimmed_count = 0;
  unsigned int bases_trimmed = 0;
  std::string read;
  
  // parse_commandline(argc, argv);
  fastx_parse_cmdline(argc, argv, "b:l:", parse_program_args);
  
  fastx_init_reader(&fastx, get_input_filename(), FASTA_OR_FASTQ, ALLOW_N,
                    REQUIRE_UPPERCASE, get_fastq_ascii_quality_offset());
  fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT,
                    compress_output_flag());
  
  while (fastx_read_next_record(&fastx)) {
    reads_count++;
    skip_read = 0;
    read = std::string(fastx.nucleotides);
    axe_at = homopolymer_start(read, NT);
    
    if (axe_at != -1 && axe_at < min_length) {
      skip_read = 1;
    } else if (axe_at < read.length()){
      // Trim the read at this position
      fastx.nucleotides[axe_at] = 0;
      trimmed_count++;
      bases_trimmed += read.length() - axe_at;
    }
    
    if (skip_read) {
      trashed_count++;
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
