#include <fplll/sieve/sieve_gauss.h>
#include <fplll/lll.h>
#include <fplll.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <cstring>
#define BLOCK 3
using namespace std;
using namespace fplll;

void DG_Sieve(int dim,int ini_dim,int round,int n,int split_no)
{
	if (dim>40)
	{
		char buf[40] = {0};
		Z_NR<mpz_t> goal_norm;  goal_norm = 0;
		//划分格基，调用三个函数，读入子块列表，输出自身列表
		//划分当前块
		//split   todo!!!!!!!!!!!!!!!!!!!!!!!!
		sprintf(buf,"bash  ./split.sh %d %d %d %d %d",dim,ini_dim,round,n,split_no);
		system(buf); 
		int split_number=1;
		DG_Sieve(dim-n,ini_dim,round+1,n,split_number++);
		DG_Sieve(dim-n,ini_dim,round+1,n,split_number++);
		DG_Sieve(dim-n,ini_dim,round+1,n,split_number++);
		//hebing   todo!!!!!!!!!!!!!!!!!!!!!!
		sprintf(buf,"bash  ./merge.sh %d %d %d %d %d",dim,ini_dim,round,n,split_no);
		system(buf);
		//当前块初始化
		sprintf(buf,"./svp_basis40-80/dim%d_%d_%d.txt",ini_dim,round,split_no);	ifstream  finBasis(buf);
		ZZ_mat<mpz_t> A;
		finBasis >> A;
		finBasis.close();
		GaussSieve<mpz_t, FP_NR<double>> gsieve(A,2,1,0);
		//读入当前块列表
		sprintf(buf,"./A%d_%d_%d.txt",ini_dim,round,split_no); ifstream yhq(buf);
		if(!yhq)     {cout<<"A.txt error2"<<round<<split_no<<endl; exit(0);}
		gsieve.scanf_list(yhq);
		yhq.close();
		if (round!=0)
		{
			gsieve.sieve(goal_norm);
			//输出当前块列表
			sprintf(buf,"./A%d_%d_%d.txt",ini_dim,round,split_no); 
			ofstream outfile(buf);
			if(!outfile) {cout<<"A.txt outfile error"<<round<<split_no<<endl; exit(0);}
			gsieve.print_list(outfile);
			outfile.close();
		}
		else
		{
			sprintf(buf,"GaussSieveshb_time_test_out_%d.txt",dim);	ofstream  Out(buf);
			clock_t start_time = clock();
			gsieve.sieve(goal_norm);
			NumVect<Z_NR<mpz_t>> v = gsieve.return_first();
			clock_t end_time = clock();
			Out << v << endl;
			Out << "The run time to solve " << dim << "-dimension  is:" << (end_time - start_time) /CLOCKS_PER_SEC <<" s  "<< (end_time - start_time) << "s!" << endl;
			Out.close();
		}
	}
	else
	{
		char buf[40] = {0};
		Z_NR<mpz_t> goal_norm;  goal_norm = 0;
		//当前块初始化
		sprintf(buf,"./svp_basis40-80/dim%d_%d_%d.txt",ini_dim,round,split_no);	ifstream  finBasis(buf);
		ZZ_mat<mpz_t> A;
		finBasis >> A;
		finBasis.close();
		GaussSieve<mpz_t, FP_NR<double>> gsieve(A,2,1,0);
		//筛
		gsieve.sieve(goal_norm);
		//输出当前块列表
		sprintf(buf,"./A%d_%d_%d.txt",ini_dim,round,split_no);
		ofstream outfile(buf);
		if(!outfile) {cout<<"A.txt outfile error"<<round<<split_no<<endl; exit(0);}
		gsieve.print_list(outfile);
		outfile.close();
	}
}
	
int main(int argc,char **argv)
{
	int dim,n;
	sscanf(argv[1],"%d",&dim);
	sscanf(argv[2],"%d",&n);
	ZZ_mat<mpz_t> *A_block=new ZZ_mat<mpz_t>[BLOCK];
	Z_NR<mpz_t>* haha;
	char buf[40] = {0};
	char buf2[40] = {0};
	int i,j;
	ZZ_mat<mpz_t> A;
    sprintf(buf2,"./svp_basis40-80/dim%d.txt",dim);
    ifstream  finalBasis(buf2);
    finalBasis >> A;
    int status;
    status = lll_reduction(A);if( status == 0)  cout << "LLL success" << endl;
    sprintf(buf,"./svp_basis40-80/dim%d_0_0.txt",dim);
    ofstream qhyfile(buf);
    qhyfile << A;
    qhyfile.close();
	DG_Sieve(dim,dim,0,n,0);
	return 0;
}

