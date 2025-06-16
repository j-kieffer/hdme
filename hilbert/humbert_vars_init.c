
#include <stdlib.h>
#include "hilbert.h"

char** humbert_vars_init(void)
{
  char** vars;
  vars = malloc(2 * sizeof(char*));
  vars[0] = malloc(2 * sizeof(char));
  vars[1] = malloc(2 * sizeof(char));
  return vars;
}
