#include "wig.h"
#include "memutil.h"
#include <stdio.h>

#define CHR22_LEN 49691432


int main(int argc, char **argv) {
  char *path = "/Users/gmcvicker/data/gene_models/average_depth/chr22.avg_depth.wig.gz";
  float *vals;

  fprintf(stderr, "reading wiggle file\n");
  vals = read_wiggle(path, CHR22_LEN);
  fprintf(stderr, "done");

  for(long i = 0; i < CHR22_LEN; i++) {
    print i;
  }

  my_free(vals);

  return 0;
}
