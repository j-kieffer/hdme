#include <stdlib.h>

#include "igusa.h"


int main()
{
  acb_ptr roots;
  acb_ptr reorder;
  slong k, perm;

  flint_printf("thomae_reorder....");
  fflush(stdout);

  roots = _acb_vec_init(6);
  reorder = _acb_vec_init(6);
  
  for (k = 0; k < 6; k++)
    {
      acb_set_si(&roots[k], k+1);
    }
  for (perm = 0; perm < 720; perm++)
    {
      thomae_reorder(reorder, roots, perm);
      /* flint_printf("perm = %wd: ", perm);
      for (k = 0; k < 6; k++) acb_printd(&reorder[k], 10);
      flint_printf("\n"); */
    }
    
  _acb_vec_clear(roots, 6);
  _acb_vec_clear(reorder, 6);
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
