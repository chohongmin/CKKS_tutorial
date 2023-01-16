#include <stdio.h>
#include "CKKS_CoeffToSlot.h"
#include "CKKS_poly.h"
#include "CKKS_EvalMod.h"

#define logN 10
#define N (1<<logN)
#define L 24
#define DNUM 4
#define K (L/DNUM)
#define Delta (1ULL<<50)


double norm_square_exp(const double zr[N/2], const double zi[N/2]){
	double sum=0;
	for(int i=0; i<N/2; i++) sum+=zr[i]*zr[i]+zi[i]*zi[i];
	return (sum/(N/2));
}

//
void main(){
	//---------------------------------------------------------------------
	// Initial setting
	//---------------------------------------------------------------------
	uint64_t q[L];	uint64_t p[K];
	{
		uint64_t q50[19]; find_RNS_primes<19>(1ULL<<50,N,q50);
		uint64_t q60[ 9]; find_RNS_primes< 9>(1ULL<<60,N,q60);
		q[0] = q60[0]; q[21] = q50[0]; p[5] = q60[8];

		for(int i= 1; i<=12; i++) q[i] = q50[i];
		for(int i=14; i<=20; i++) q[i] = q60[1+i-14];
		q[23]=q50[13];
		for(int i= 0; i<= 4; i++) p[i] = q50[14+i];
		find_RNS_primes<1>((uint64_t)((1ULL<<50)*12),N,&q[22]);
		find_RNS_primes<1>(            1ULL<<40     ,N,&q[13]);
	}
	int h=64; int s[N];
	keygen<N>(h,s);
	double zr[N/2], zi[N/2];
	for(int i=0;i<N/2;i++){
		zr[i]=((double)rand())/RAND_MAX/sqrt(2);
		zi[i]=((double)rand())/RAND_MAX/sqrt(2);
	}
	uint64_t pt[1][N];
	encode<N,logN,1>(zr,zi,Delta,q,pt);
	uint64_t ct_bot_hat[2][1][N];
	enc<N,1>(pt,s,q,ct_bot_hat);
	//---------------------------------------------------------------------
	// ModUp
	//---------------------------------------------------------------------
	uint64_t ct_hat[2][L][N];
	for(int i=0;i<2;i++)
	for(int j=0;j<N;j++)
		ct_hat[i][0][j] = ct_bot_hat[i][0][j];
	intt<N,1>(q,ct_hat[0]);
	intt<N,1>(q,ct_hat[1]);
	modUP<N, L>(q, 1, ct_hat[0]);
	modUP<N, L>(q, 1, ct_hat[1]);
	ntt<N,L>(q, ct_hat[0]);
	ntt<N,L>(q, ct_hat[1]);
	//----------------------------------------------------
	// debug : pt = ifft(z), pt+qI = ifft(decode(dec(ct))
	//----------------------------------------------------
	{
		double pt_real[N]; ifft<N,logN>(zr,zi,pt_real);
		uint64_t pt[L][N]; dec<N,L>(ct_hat,s,q,pt);
		double wr[N/2],wi[N/2]; decode<N,logN,L>(pt,Delta,q,wr,wi);
		double pt_qI_real[N]; ifft<N,logN>(wr,wi,pt_qI_real);
		double Imax=0;
		for(int i=0;i<N;i++){
			double I = (pt_qI_real[i]*Delta-pt_real[i]*Delta)/q[0];
			if( fabs(I)>Imax)
				Imax = fabs(I);
		}
		printf("%f\n", Imax);
	}
	//---------------------------------------------------------------------
	// CoeffToSlot
	//---------------------------------------------------------------------
	SparseComplexMatrix<N/2, 27> A[3];
	splitU0R_logN_10( A );
	uint64_t  evk_hat       [DNUM][2][DNUM*K+K][N];
	uint64_t ckey_hat       [DNUM][2][DNUM*K+K][N];
	uint64_t rkey_hat[3][27][DNUM][2][DNUM*K+K][N];
	swkgen_logN_10<L,DNUM,K>(q,p,s,A,evk_hat,ckey_hat,rkey_hat);
	uint64_t ct_cts_hat[2][2][L-3][N];
	CoeffToSlot_logN_10<L,DNUM,K>(q,p,Delta,A,ckey_hat,rkey_hat,ct_hat,ct_cts_hat);
	//-----------------------------------------------
	// debug
	//-----------------------------------------------
	double x                    [2][N/2];
	double pt_over_q_bitReversed[N];
	{
		uint64_t pt_cts[2][L-3][N];
		double wr[2][N/2], wi[2][N/2];
		dec<N, L-3>(ct_cts_hat[0], s, q, pt_cts[0]);
		dec<N, L-3>(ct_cts_hat[1], s, q, pt_cts[1]);
		decode<N, logN, L-3>(pt_cts[0], q[0], q, wr[0], wi[0]);
		decode<N, logN, L-3>(pt_cts[1], q[0], q, wr[1], wi[1]);
		double pt_over_Delta[N];  ifft<N,logN>(zr,zi,pt_over_Delta);
		double pt_over_Delta_bitReversed[N];
		bitReverse<N/2>(pt_over_Delta    , pt_over_Delta_bitReversed    );
		bitReverse<N/2>(pt_over_Delta+N/2, pt_over_Delta_bitReversed+N/2);
		double er[2][N/2], ei[2][N/2];
		for(int i=0;i<N/2;i++){
			er[0][i] = wr[0][i]-pt_over_Delta_bitReversed[i    ]*Delta/q[0]/12;
			ei[0][i] = wi[0][i]-0;
			er[1][i] = wr[1][i]-pt_over_Delta_bitReversed[i+N/2]*Delta/q[0]/12;
			ei[1][i] = wi[1][i]-0;

			pt_over_q_bitReversed[i    ]=pt_over_Delta_bitReversed[i    ]/(1<<10);
			pt_over_q_bitReversed[i+N/2]=pt_over_Delta_bitReversed[i+N/2]/(1<<10);
			x[0][i]=wr[0][i]*12;
			x[1][i]=wr[1][i]*12;
		}
	}
	//---------------------------------------------------------------------
	// EvalMod
	//---------------------------------------------------------------------
	uint64_t ct_hat_evalmod[2][2][L-10][N];
	
	EvalMod_h_64<N,L-3,DNUM,K>(q,p,q[0],evk_hat,ct_cts_hat[0],ct_hat_evalmod[0]);
	EvalMod_h_64<N,L-3,DNUM,K>(q,p,q[0],evk_hat,ct_cts_hat[1],ct_hat_evalmod[1]);

	//-----------------------------------------------
	// debug
	//-----------------------------------------------
	{
		uint64_t pt_evalmod[2][L-10][N];
		double wr[2][N/2], wi[2][N/2];
		dec<N, L-10>(ct_hat_evalmod[0], s, q, pt_evalmod[0]);
		dec<N, L-10>(ct_hat_evalmod[1], s, q, pt_evalmod[1]);
		decode<N, logN, L-10>(pt_evalmod[0], q[0], q, wr[0], wi[0]);
		decode<N, logN, L-10>(pt_evalmod[1], q[0], q, wr[1], wi[1]);
		double er[2][N/2], ei[2][N/2];
		for(int i=0;i<N/2;i++){
			er[0][i] = wr[0][i]-pt_over_q_bitReversed[i    ];
			ei[0][i] = wi[0][i]-0;
			er[1][i] = wr[1][i]-pt_over_q_bitReversed[i+N/2];
			ei[1][i] = wi[1][i]-0;
		}
		printf("%e,  %e\n",norm_square_exp(er[0],ei[0]),
			               norm_square_exp(er[1],ei[1]));
	}
}
