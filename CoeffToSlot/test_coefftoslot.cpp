#include <stdio.h>
#include "CKKS_CoeffToSlot.h"

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

void main(){
	uint64_t q[L]; 	find_RNS_primes<L>(  Delta,N,q);
	uint64_t p[K]; 	find_RNS_primes<K>( q[L-1],N,p);
	//-----------------------------------------------
	int h=3; int s[N]; uint64_t pk_hat[2][L][N];
	keygen<N>(h,s);
	keygen<N,L>(s,q,pk_hat);
	//-----------------------------------------------
	double zr[N/2], zi[N/2];
	for(int i=0;i<N/2;i++){
		zr[i]=((double)rand())/RAND_MAX/sqrt(2);
		zi[i]=((double)rand())/RAND_MAX/sqrt(2);
	}
	//-----------------------------------------------
	uint64_t pt[L][N];
	encode<N,logN,L>(zr,zi,Delta,q,pt);
	uint64_t ct_hat[2][L][N];
	enc<N,L>(pt,s,q,ct_hat);
	
	//-----------------------------------------------
	SparseComplexMatrix<N/2, 27> A[3];
	splitU0R_logN_10( A );

	uint64_t  evk_hat       [DNUM][2][DNUM*K+K][N];
	uint64_t ckey_hat       [DNUM][2][DNUM*K+K][N];
	uint64_t rkey_hat[3][27][DNUM][2][DNUM*K+K][N];
	swkgen_logN_10<L,DNUM,K>(q,p,s,A,evk_hat,ckey_hat,rkey_hat);
	//-----------------------------------------------
	uint64_t ct_cts_hat[2][2][L-3][N];
	CoeffToSlot_logN_10<L,DNUM,K>(q,p,Delta,A,ckey_hat,rkey_hat,ct_hat,ct_cts_hat);
	//-----------------------------------------------
	uint64_t pt_cts[2][L-3][N];
	double wr[2][N/2], wi[2][N/2];
	dec<N, L-3>(ct_cts_hat[0], s, q, pt_cts[0]);
	dec<N, L-3>(ct_cts_hat[1], s, q, pt_cts[1]);
	decode<N, logN, L-3>(pt_cts[0], Delta, q, wr[0], wi[0]);
	decode<N, logN, L-3>(pt_cts[1], Delta, q, wr[1], wi[1]);
	//-----------------------------------------------
	double pt_over_Delta[N];
	ifft<N,logN>(zr,zi,pt_over_Delta);
	double pt_over_Delta_bitReversed[N];
	bitReverse<N/2>(pt_over_Delta    , pt_over_Delta_bitReversed    );
	bitReverse<N/2>(pt_over_Delta+N/2, pt_over_Delta_bitReversed+N/2);
	double er[2][N/2], ei[2][N/2];
	for(int i=0;i<N/2;i++){
		er[0][i] = wr[0][i]-pt_over_Delta_bitReversed[i];
		ei[0][i] = wi[0][i]-0;
		er[1][i] = wr[1][i]-pt_over_Delta_bitReversed[i+N/2];
		ei[1][i] = wi[1][i]-0;
	}
	printf("%e,%e\n",norm_square_exp(er[0],ei[0]),
		             norm_square_exp(er[1],ei[1]));
}