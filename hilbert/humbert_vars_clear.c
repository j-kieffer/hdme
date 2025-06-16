
#include <stdlib.h>
#include "hilbert.h"

void humbert_vars_clear(char** vars)
{
  free(vars[0]);
  free(vars[1]);
  free(vars);
}
