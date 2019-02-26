#include <fplll/sieve/sieve_gauss.h>
#include <fplll/lll.h>
#include <fplll.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <cstring>

using namespace std;
using namespace fplll;

int main(int argc,char **argv)
{
	int dim;sscanf(argv[1],"%d",&dim);
// put the basis into matrix A
	ZZ_mat<mpz_t> A;
	char buf[40] = {0};
	char buf2[40] = {0};
	sprintf(buf2,"GaussSieveszd_time_test_out_%d.txt",dim);
	ofstream  Out(buf2);
		sprintf(buf,"./svp_basis40-80/dim%d.txt",dim);
		ifstream  finBasis(buf);
		finBasis >> A;
		int status;

		status = lll_reduction(A);
		if( status == 0)
			cout << "LLL success" << endl;
		Out << A << endl;

		//use bkz to solve and record the time 
		clock_t start_time = clock();

		GaussSieve<mpz_t, FP_NR<double>> gsieve( A, 2, 1, 0 );
		Z_NR<mpz_t> goal_norm;
		goal_norm = 0;
		gsieve.sieve(goal_norm);
		NumVect<Z_NR<mpz_t>> v = gsieve.return_first();


		clock_t end_time = clock();

		Out << v << endl;

		Out << "The run time to solve " << dim  << "-dimension  is:" << (end_time - start_time) << "s!" << endl;
	
	return 0;
}
