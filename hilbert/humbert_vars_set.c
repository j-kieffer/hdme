
#include <string.h>
#include "hilbert.h"

void humbert_vars_set(char** vars, slong delta)
{
  if ((delta == 5) ||
      (delta == 13) ||
      (delta == 17) ||
      (delta == 53) ||
      (delta == 56) ||
      (delta == 57) ||
      (delta == 60) ||
      (delta == 61)
      )
    {
      strcpy(vars[0], "g");
      strcpy(vars[1], "h");
    }
  else if ((delta == 12) ||
	   (delta == 40) ||
	   (delta == 85)
	   )
    {
      strcpy(vars[0], "e");
      strcpy(vars[1], "f");
    }
  else if (delta == 24)
    {
      strcpy(vars[0], "a");
      strcpy(vars[1], "d");
    }
  else if ((delta == 28) ||
	   (delta == 29) ||
	   (delta == 37)
	   )
    {
      strcpy(vars[0], "f");
      strcpy(vars[1], "g");
    }
  else if (delta == 93)
    {
      strcpy(vars[0], "m");
      strcpy(vars[1], "n");
    }
  else
    {
      strcpy(vars[0], "r");
      strcpy(vars[1], "s");
    }
}
