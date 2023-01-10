#include <stdio.h>
#include "CKKS_lineartransform.h"

#define logN 10
#define N (1<<logN)
#define L 9
#define DNUM 3
#define K (L/DNUM)
#define Delta (1ULL<<40)
#define NumDiags 3


double norm_square_exp(const double zr[N/2], const double zi[N/2]){
	double sum=0;
	for(int i=0; i<N/2; i++) sum+=zr[i]*zr[i]+zi[i]*zi[i];
	return (sum/(N/2));
}

double norm_max(const double zr[N/2], const double zi[N/2]){
	double max=0;
	for(int i=0; i<N/2; i++){ double abs=sqrt(zr[i]*zr[i]+zi[i]*zi[i]);
		if(abs>max) max=abs;
	}
	return max;
}

void main(){
	uint64_t q[L]; 	find_RNS_primes<L>(  Delta,N,q);
	uint64_t p[K]; 	find_RNS_primes<K>( q[L-1],N,p);
	//-----------------------------------------------
	int h=3; int s[N]; uint64_t pk_hat[2][L][N];
	keygen<N>(h,s);
	keygen<N,L>(s,q,pk_hat);
	int ss[N]; conv<N>(s,s,ss);

	uint64_t evk_hat[DNUM][2][DNUM*K+K][N];
	swkgen<N,L,DNUM>(ss,s,q,p,evk_hat);
	//-----------------------------------------------
	double zr[N/2], zi[N/2];
	for(int i=0;i<N/2;i++){
		zr[i]=((double)rand())/RAND_MAX;
		zi[i]=((double)rand())/RAND_MAX;
	}
	//-----------------------------------------------
	uint64_t pt[L][N];
	encode<N,logN,L>(zr,zi,Delta,q,pt);
	uint64_t ct_hat[2][L][N];
	enc<N,L>(pt,s,q,ct_hat);
	
	//-----------------------------------------------
	SparseComplexMatrix<N/2, 3> A;
	A.shift[0] = 0;
	A.shift[1] = 1;
	A.shift[2] = N/2 - 1;
	
	for(int i=0; i<3  ; i++)
	for(int j=0; j<N/2; j++){
		A.diagr[i][j] = ((double)rand())/RAND_MAX;
		A.diagi[i][j] = ((double)rand())/RAND_MAX;
	}

	uint64_t rkey_hat[3][DNUM][2][DNUM*K+K][N];
	rkey_gen<N, L, DNUM, K, NumDiags>(A, q, p, s, rkey_hat);

	uint64_t res_hat[2][L][N];
	lineartransform<N, logN, L, DNUM, K, NumDiags>(A, Delta, q, p, rkey_hat, ct_hat, res_hat);

	uint64_t res_rs_hat[2][L-1][N];
	RS_hat<N, L>(q, res_hat[0], res_rs_hat[0]);
	RS_hat<N, L>(q, res_hat[1], res_rs_hat[1]);
	//-----------------------------------------------
	uint64_t pt1[L-1][N];
	dec<N, L-1>(res_rs_hat, s, q, pt1);

	double er[N/2], ei[N/2];
	decode<N, logN, L-1>(pt1, Delta, q, er, ei);
	for(int i=0; i<N/2; i++){
		int im1 = (i==0) ? (i-1 + N/2) : (i-1);
		int ip1 = (i==N/2-1) ? (i+1 - N/2) : (i+1);
		er[i] -= A.diagr[0][i] * zr[i] + A.diagr[1][i]*zr[ip1] + A.diagr[2][i]*zr[im1];
		ei[i] -= A.diagr[0][i] * zi[i] + A.diagr[1][i]*zi[ip1] + A.diagr[2][i]*zi[im1];
		er[i] += A.diagi[0][i] * zi[i] + A.diagi[1][i]*zi[ip1] + A.diagi[2][i]*zi[im1];
		ei[i] -= A.diagi[0][i] * zr[i] + A.diagi[1][i]*zr[ip1] + A.diagi[2][i]*zr[im1];
	}

	//-----------------------------------------------
	printf("%e\n",norm_square_exp(er,ei));
	

}