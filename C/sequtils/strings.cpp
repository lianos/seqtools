#include "sequtils/strings.h"

char* substring(const char* str, size_t begin, size_t len) {
  char *subs;
  if (str == 0 || strlen(str) == 0 || strlen(str) < begin 
      || strlen(str) < (begin+len)) {
    return 0;
  }
  subs = (char *) malloc(sizeof(char) * len);
  strncpy(subs, str + begin, len);
  return subs;
} 
