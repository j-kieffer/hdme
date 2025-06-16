#include "hdme_data.h"
#include <string.h>

void hdme_data_vars_set(char** vars, const char* name, slong k)
{
  strcpy(vars[k], name);
}
