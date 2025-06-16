#include <stdlib.h>

#include "modular.h"


int main()
{
	slong iter;

	flint_printf("siegel_2step_direct_isog_Q....");
	fflush(stdout);

	/* Up to iter < 4 */
	for (iter = 0; iter < 6; iter++)
	{
		fmpz* I1;
		fmpz* I2;
		fmpz* all_I;
		slong ell;
		slong max_nb_roots = 2;
		slong nb_roots = 1;
		slong k;
		int res;

		I1 = _fmpz_vec_init(4);
		I2 = _fmpz_vec_init(4);
		all_I = _fmpz_vec_init(4 * max_nb_roots);

		switch(iter) {
			case 0:
				// something in the class https://beta.lmfdb.org/Genus2Curve/Q/63749/a/
				fmpz_set_si(&I1[0], 19864);
				fmpz_set_si(&I1[1], 9722500);
				fmpz_set_si(&I1[2], 61093064375);
				fmpz_set_si(&I1[3], 62254882812500);
				ell = 2;
				fmpz_set_si(&I2[0], 1204816);
				fmpz_set_si(&I2[1], 46060937500);
				fmpz_set_si(&I2[2], 10446228721484375);
				// -744748353958129882812500 == -2^40 * 677344682079  - 31376986196)
				fmpz_set_si(&I2[3], -677344682079);
				fmpz_mul_2exp(&I2[3], &I2[3], 40);
				fmpz_add_si(&I2[3], &I2[3], -31376986196);
				break;

			case 1:
				// something in the class https://beta.lmfdb.org/Genus2Curve/Q/63749/a/
				fmpz_set_si(&I1[0], 19864);
				fmpz_set_si(&I1[1], 9722500);
				fmpz_set_si(&I1[2], 61093064375);
				fmpz_set_si(&I1[3], 62254882812500);
				ell = 2;
				// https://beta.lmfdb.org/Genus2Curve/Q/63749/a/63749/1
				fmpz_set_si(&I2[0], 19024);
				fmpz_set_si(&I2[1], 2911276);
				fmpz_set_si(&I2[2], 17003250887);
				fmpz_set_si(&I2[3], 254996);
				break;

			case 2:
				/* https://beta.lmfdb.org/Genus2Curve/Q/249/a/249/1,
https://beta.lmfdb.org/Genus2Curve/Q/249/a/6723/1 */
				fmpz_set_si(&I1[0], 108);
				fmpz_set_si(&I1[1], 57);
				fmpz_set_si(&I1[2], 2259);
				fmpz_set_si(&I1[3], -31872);
				ell = 2;
				fmpz_set_si(&I2[0], 1932);
				fmpz_set_si(&I2[1], 87897);
				fmpz_set_si(&I2[2], 65765571);
				fmpz_set_si(&I2[3], 860544);
				break;

			case 3:
				/* https://beta.lmfdb.org/Genus2Curve/Q/277/a/277/1,
https://beta.lmfdb.org/Genus2Curve/Q/277/a/277/2 */
				fmpz_set_si(&I1[0], 64);
				fmpz_set_si(&I1[1], 352);
				fmpz_set_si(&I1[2], 9552);
				fmpz_set_si(&I1[3], -1108);
				ell = 3;
				fmpz_set_si(&I2[0], 4480);
				fmpz_set_si(&I2[1], 1370512);
				fmpz_set_si(&I2[2], 1511819744);
				fmpz_set_si(&I2[3], -1108);
				break;

			case 4:
				/* https://beta.lmfdb.org/Genus2Curve/Q/523/a/523/1,
https://beta.lmfdb.org/Genus2Curve/Q/523/a/523/2 */
				fmpz_set_si(&I1[0], 120);
				fmpz_set_si(&I1[1], -540);
				fmpz_set_si(&I1[2], -29169);
				fmpz_set_si(&I1[3], -2092);
				ell = 5;
				fmpz_set_si(&I2[0], 332400);
				fmpz_set_si(&I2[1], 10084860);
				fmpz_set_si(&I2[2], 1107044456391);
				fmpz_set_si(&I2[3], -2092);
				break;

			case 5:
				/* https://beta.lmfdb.org/Genus2Curve/Q/295/a/295/1,
https://beta.lmfdb.org/Genus2Curve/Q/295/a/295/2 */
				fmpz_set_si(&I1[0], 108);
				fmpz_set_si(&I1[1], -39);
				fmpz_set_si(&I1[2], 20835);
				fmpz_set_si(&I1[3], 37760);
				ell = 7;
				fmpz_set_si(&I2[0], 198804);
				fmpz_set_si(&I2[1], 305807001);
				fmpz_set_si(&I2[2], 18482629056189);
				fmpz_set_si(&I2[3], -37760);
				break;



			default:
				flint_printf("iteration = %wd has not been initialized\n", iter);
				flint_abort();
		}

		igusa_from_IC_fmpz(I1, I1);
		igusa_from_IC_fmpz(I2, I2);

		siegel_2step_direct_isog_Q(&nb_roots, all_I, I1, ell);

		if (nb_roots == 0)
		{
			flint_printf("FAIL (roots)\n");
			flint_printf("nb_roots = %wd\n", nb_roots);
			fflush(stdout);
			flint_abort();
		}

		res = 0;
		for (k = 0; k < nb_roots; k++)
		{
			if (_fmpz_vec_equal(I2, &all_I[4*k], 4)) res = 1;
		}
		if (!res)
		{
			flint_printf("FAIL (values)\n");
			fflush(stdout);
			flint_abort();
		}

		_fmpz_vec_clear(I1, 4);
		_fmpz_vec_clear(I2, 4);
		_fmpz_vec_clear(all_I, 4 * max_nb_roots);
	}

	flint_cleanup();
	flint_printf("PASS\n");
	return EXIT_SUCCESS;
}
