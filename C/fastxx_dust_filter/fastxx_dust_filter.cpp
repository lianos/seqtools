// This tool depends on the fastx toolkit.
// Make sure libfastx and gtextutils is in your -I(include) path when compiling.
//
// Use it to strip stretches of homopolymers at the end of a sequence.
//
// To install:
//   * Ensure that gtextutils and fastx is installed
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
"usage: fastxx_dust_filter [-h] [-v] [-b A] [-l N] [-t trashfile] " \
"[-i INFILE] [-o OUTFILE]\n" \
"by Steve Lianoglou (slianoglou@gmail.com)\n" \
"\n" \
"   [-h]         = Help me help you.\n" \
"   [-v]         = verbose: print short summary of input/output counts\n"      \
"   [-d FLOAT]   = (float) value for dust score [4.0]\n"                   \
"   [-t outfile] = The path to the file to save 'trashed' reads in. If the\n"  \
"                  filename ends in *.gz, the file will be compressed.\n"      \
"   [-l N]       = discard sequences shorter than N nucleotids [5]\n"          \
"   [-i INFILE]  = FASTA/Tabular input file. default is STDIN.\n"              \
"   [-o OUTFILE] = FASTA/Tabular output file. default is STDOUT.\n"            \
"\n";


// Working parameters
double dust_threshold = 4.0f;
unsigned int min_length = 5;
int save_trash = 0;
char trashfile[MAX_TRASHFILE_LEN] = "dust.trash.fastq";
int compress_trash = 0;
int trim_front = 0;

FASTX fastx;
FASTX trashx;

// Dust scoring scheme lifted from `sga`:
//   https://github.com/jts/sga/blob/master/src/Util/Util.cpp
//
// As explained in Morgulis A. J Comp Bio:
//   A fast and symmetric DUST implementation to Mask Low-Complexity DNA
//   Sequences.
double calculate_dust_score(const std::string& read);
double calculate_dust_score(const std::string& read) {
  std::map<std::string, int> scoreMap;
  
  // Slide a 3-mer window over the sequence and insert the sequences into the map
  for (size_t i = 0; i < read.size() - 3; ++i) {
    std::string triMer = read.substr(i, 3);
    scoreMap[triMer]++;
  }

  // Calculate the score by summing the square of every element in the map
  double sum = 0;
  std::map<std::string, int>::iterator iter = scoreMap.begin();
  for (; iter != scoreMap.end(); ++iter) {
    int tc = iter->second;
    double score = (double)(tc * (tc - 1)) / 2.0f;
    sum += score;
  }
  
  return sum / (seq.size() - 2);
}


int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg) {
  switch(optc) {
    case 'd':
      if (optarg == NULL) {
        errx(1, "[-d] parameter requires an argument value");
      }
      if (strlen(optarg) != 1) {
        errx(1, "[-b] can only be a character of length 1, eg: A, C, G, or T");
      }
      dust_threshold = strtod(optarg, NULL);
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

int main(int argc, char* argv[]) {
  int dusted;
  int too_short;
  double dust_score;
  unsigned int reads_count = 0;
  unsigned int read_length = 0;
  unsigned int trashed_count = 0;
  unsigned int trimmed_count = 0;
  unsigned int i;

  // unsigned int (*posf)(const std::string&, char) = NULL;
  std::string read;

  // parse_commandline(argc, argv);
  fastx_parse_cmdline(argc, argv, "d:l:t:", parse_program_args);

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
    if (save_trash) {
      fastx_read_next_record(&trashx);
    }
    reads_count++;
    dusted = 0;
    
    read = std::string(fastx.nucleotides);
    read_length = read.length();
    too_short = read_length < min_length;
    
    if (too_short) {
      dusted = 1;
    } else {
      dust_score = calculate_dust_score(read);
      dusted = dust_score >= dust_threshold;
    }
    
    if (dusted) {
      trashed_count++;
      if (save_trash) {
        fastx_write_record(&trashx);
      }
    } else {
      fastx_write_record(&fastx);
    }
  }
  
  // Print verbose report
  if (verbose_flag()) {
    fprintf(get_report_file(), "Number of reads processed: %d\n", reads_count);
    fprintf(get_report_file(), "Number of reads dusted: %d\n", trashed_count);
  }
  
  return 0;
}
