#include <stdio.h>
#include "CKKS_ks_wih_dnum.h"

#define logN 10
#define N (1<<logN)
#define L 9
#define DNUM 3
#define K (L/DNUM)
#define Delta (1ULL<<40)


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
	int h=3; int sfr[N], sto[N]; uint64_t pk_hat[2][L][N];
	keygen<N>(h,sfr);
	keygen<N>(h,sto);
	keygen<N,L>(sfr,q,pk_hat);

	uint64_t swk_hat[DNUM][2][L+K][N];
	swkgen<N,L,DNUM>(sfr,sto,q,p,swk_hat);
	//-----------------------------------------------
	double zr[N/2], zi[N/2];
	for(int i=0;i<N/2;i++){
		zr[i]=((double)rand())/RAND_MAX/sqrt(2);
		zi[i]=((double)rand())/RAND_MAX/sqrt(2);
	}
	//-----------------------------------------------
	uint64_t pt[L-2][N];
	encode<N,logN,L-2>(zr,zi,Delta,q,pt);
	uint64_t ct_hat[2][L-2][N];
	enc<N,L-2>(pt,sfr,q,ct_hat);
	//-----------------------------------------------
	uint64_t out_hat[2][L-2][N];
	ks<N,L-2,DNUM,K>(q,p,swk_hat,ct_hat,out_hat);
	//-----------------------------------------------
	dec<N, L-2>(out_hat, sto, q, pt);
	double er[N/2],ei[N/2];
	decode<N, logN, L-2>(pt, Delta, q, er, ei);
	for(int i=0;i<N/2;i++){er[i]-=zr[i]; ei[i]-=zi[i];}
	//-----------------------------------------------
	double sig_sqr = 24;
	printf("%e, %e, %e\n",norm_square_exp(er,ei),
		((double)N)*N/Delta/Delta*sig_sqr/2 +
		((double)N)  /Delta/Delta*(1./12+sig_sqr*(1+h)), norm_max(er,ei));

}