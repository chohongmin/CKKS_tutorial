#include <stdio.h>
#include "CKKS_ks_wih_dnum.h"
#include "CKKS_poly.h"

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
		zi[i]=0;
	}
	//-----------------------------------------------
	uint64_t pt[L][N];
	encode<N,logN,L>(zr,zi,Delta,q,pt);
	uint64_t ct_hat[2][L][N];
	enc<N,L>(pt,s,q,ct_hat);
	//-----------------------------------------------
	uint64_t T2_hat[2][L-1][N];
	Chebyshev_doubling<N,L,DNUM,K>(q,p,Delta,ct_hat,evk_hat,T2_hat);
	uint64_t T4_hat[2][L-2][N];
	Chebyshev_doubling<N,L-1,DNUM,K>(q,p,Delta,T2_hat,evk_hat,T4_hat);
	uint64_t T8_hat[2][L - 3][N];
	Chebyshev_doubling<N, L - 2, DNUM, K>(q, p, Delta, T4_hat, evk_hat, T8_hat);
	uint64_t T16_hat[2][L - 4][N];
	Chebyshev_doubling<N, L - 3, DNUM, K>(q, p, Delta, T8_hat, evk_hat, T16_hat);
	uint64_t T32_hat[2][L - 5][N];
	Chebyshev_doubling<N, L - 4, DNUM, K>(q, p, Delta, T16_hat, evk_hat, T32_hat);
	//-----------------------------------------------
	double c[64];
	for (int i = 0; i < 64; i++) c[i] = 0.1 * (i + 1);

	uint64_t res_hat[2][L-6][N];
	eval_poly_deg64<N,L,DNUM,K>(q,p,c,Delta,ct_hat,T2_hat,T4_hat,T8_hat, T16_hat, T32_hat, evk_hat,res_hat);
	//-----------------------------------------------
	dec<N, L-6>(res_hat, s, q, pt);
	double er[N/2],ei[N/2];
	decode<N, logN, L-6>(pt, Delta, q, er, ei);
	for(int i=0;i<N/2;i++){
		double T1 = zr[i];
		double T2 = 2*T1*T1-1;
		double T4 = 2*T2*T2-1;
		double T8 = 2 * T4 * T4 - 1;
		double T16 = 2 * T8 * T8 - 1;
		double T32 = 2 * T16 * T16 - 1;


		er[i] -= c[0] + c[1] * T1 + T2 * (c[2] + c[3] * T1)
			+ T4 * (c[4] + c[5] * T1 + T2 * (c[6] + c[7] * T1))
			+ T8 * (c[8] + c[9] * T1 + T2 * (c[10] + c[11] * T1)
				    + T4 * (c[12] + c[13] * T1 + T2 * (c[14] + c[15] * T1)))
			+ T16 * (c[16] + c[17] * T1 + T2 * (c[18] + c[21] * T1)
					+ T4 * (c[20] + c[21] * T1 + T2 * (c[22] + c[23] * T1))
					+ T8 * (c[24] + c[25] * T1 + T2 * (c[26] + c[27] * T1)
							+ T4 * (c[28] + c[29] * T1 + T2 * (c[30] + c[31] * T1))))
			+ T32 * (c[32] + c[33] * T1 + T2 * (c[34] + c[35] * T1)
					+ T4 * (c[36] + c[37] * T1 + T2 * (c[38] + c[39] * T1))
					+ T8 * (c[40] + c[41] * T1 + T2 * (c[42] + c[43] * T1)
							+ T4 * (c[44] + c[45] * T1 + T2 * (c[14+32] + c[15+32] * T1)))
							+ T16 * (c[16+32] + c[17+32] * T1 + T2 * (c[18+32] + c[21+32] * T1)
									+ T4 * (c[20+32] + c[21+32] * T1 + T2 * (c[22+32] + c[23+32] * T1))
									+ T8 * (c[24+32] + c[25+32] * T1 + T2 * (c[26+32] + c[27+32] * T1)
											+ T4 * (c[28+32] + c[29+32] * T1 + T2 * (c[30+32] + c[31+32] * T1)))));
		ei[i]-=0;
	}
	//-----------------------------------------------
	double sig_sqr = 24;
	printf("%e\n",norm_square_exp(er,ei));

}