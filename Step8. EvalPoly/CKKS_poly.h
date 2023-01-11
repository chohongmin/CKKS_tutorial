#pragma once
#include <stdint.h>
#include <inttypes.h>
#include <cassert>
#include "DFT.h"
#include "NumberTheory.h"
#include "BigInt.h"

//
template<int N, int L, int DNUM, int K>
void Chebyshev_doubling( const uint64_t q[L],
						 const uint64_t p[K], uint64_t Delta,
						 const uint64_t T[2][L][N],
						 const uint64_t evk_hat[DNUM][2][DNUM*K+K][N],
							   uint64_t Tnext[2][L-1][N]){
	uint64_t temp[2][L][N];
	mul<N,L,DNUM,K>(q, p, evk_hat, T, T, temp);
	for(int i=0; i<L; i++)
	for(int j=0; j<N; j++){
		temp[0][i][j] = (2*temp[0][i][j]) % q[i];
		temp[1][i][j] = (2*temp[1][i][j]) % q[i];
	}
	RS_hat<N,L>(q, temp[0], Tnext[0]);
	RS_hat<N,L>(q, temp[1], Tnext[1]);

	///*
	uint64_t scale = Delta;//uint64_t((((double)Delta)*Delta)/q[L-1]+0.5);
	for(int i=0; i<L-1; i++)
	for(int j=0; j<N  ; j++){
		Tnext[0][i][j] = (Tnext[0][i][j] + q[i] - (scale%q[i])) % q[i];
	}//*/
}

//
template<int N, int L>
void eval_poly_deg2( const uint64_t q[L],
					 const double c[2], uint64_t Delta,
					 const uint64_t      T1[2][L][N],
						   uint64_t res_hat[2][L-1][N]){
	
	uint64_t L0[2][L][N];
	for(int i=0; i<2; i++){
	for(int j=0; j<L; j++){
		for(int k=1; k<N; k++)
			L0[i][j][k] = 0;
		if(c[i] >= 0)
			L0[i][j][0] = uint64_t (c[i]*Delta+.5)%q[j];
		else
			L0[i][j][0] = (q[j] - (uint64_t(-c[i]*Delta+.5))%q[j]);
		
	}}

	//----------------------------------
	// debug
	//----------------------------------
	//double tempr[N/2],tempi[N/2];
	//decode<1<<10,10,8>(L0[0],Delta,q,tempr,tempi); double c0=tempr[0];
	//decode<1<<10,10,8>(L0[1],Delta,q,tempr,tempi); double c1=tempr[0];

	ntt<N,L>(q,L0[0]);
	ntt<N,L>(q,L0[1]);
	
	{
		uint64_t temp[2][L][N];
		for(int j=0; j<L; j++)
			for(int k=0; k<N; k++){
				temp[0][j][k] = (mul_mod(L0[0][j][k], Delta % q[j], q[j]) + mul_mod(L0[1][j][k], T1[0][j][k], q[j]))%q[j];
				temp[1][j][k] =                                             mul_mod(L0[1][j][k], T1[1][j][k], q[j]);
			}
		RS_hat<N, L>(q, temp[0], res_hat[0]);
		RS_hat<N, L>(q, temp[1], res_hat[1]);
	}
}

//
template<int N, int L, int DNUM, int K>
void eval_poly_merge( const uint64_t q[L],
					  const uint64_t p[K], uint64_t Delta,
					  const uint64_t T[2][L][N],
					  const uint64_t evk_hat[DNUM][2][DNUM*K+K][N],
					  const uint64_t ct1_hat[2][L][N],
					  const uint64_t ct2_hat[2][L][N],
						    uint64_t ct_res_hat[2][L-1][N]){
	
	uint64_t temp[2][L][N];
	mul<N, L,DNUM,K>(q, p, evk_hat, T, ct2_hat, temp);
	for(int j=0; j<L; j++)
	for(int k=0; k<N; k++){
		temp[0][j][k] = (temp[0][j][k] + mul_mod(ct1_hat[0][j][k], Delta%q[j], q[j])) % q[j];
		temp[1][j][k] = (temp[1][j][k] + mul_mod(ct1_hat[1][j][k], Delta%q[j], q[j])) % q[j];
	}
	RS_hat<N, L>(q, temp[0], ct_res_hat[0]);
	RS_hat<N, L>(q, temp[1], ct_res_hat[1]);
}



//
template<int N, int L, int DNUM, int K>
void eval_poly_deg4( const uint64_t q[L],
					 const uint64_t p[K],
					 const double c[4], uint64_t Delta,
					 const uint64_t T1[2][L  ][N],
					 const uint64_t T2[2][L-1][N],
					 const uint64_t evk_hat[DNUM][2][DNUM*K+K][N],
						   uint64_t res_hat[2][L-2][N]){
	uint64_t L1[2][2][L-1][N];
	eval_poly_deg2<N,L>(q, c  , Delta, T1, L1[0]);
	eval_poly_deg2<N,L>(q, c+2, Delta, T1, L1[1]);

	eval_poly_merge<N, L-1,DNUM,K>(q, p, Delta, T2, evk_hat, L1[0], L1[1], res_hat);
}

//
template<int N, int L, int DNUM, int K>
void eval_poly_deg8( const uint64_t q[L],
					 const uint64_t p[K],
					 const double c[8], uint64_t Delta,
					 const uint64_t T1[2][L  ][N],
					 const uint64_t T2[2][L-1][N],
					 const uint64_t T4[2][L-2][N],
					 const uint64_t evk_hat[DNUM][2][DNUM*K+K][N],
						   uint64_t res_hat[2][L-3][N]){
	uint64_t L2[2][2][L-2][N];
	eval_poly_deg4<N,L,DNUM,K>(q, p, c  , Delta, T1, T2, evk_hat, L2[0]);
	eval_poly_deg4<N,L,DNUM,K>(q, p, c+4, Delta, T1, T2, evk_hat, L2[1]);

	eval_poly_merge<N, L-2,DNUM,K>(q, p, Delta, T4, evk_hat, L2[0], L2[1], res_hat);
}

//
template<int N, int L, int DNUM, int K>
void eval_poly_deg16(const uint64_t q[L],
					 const uint64_t p[K],
					 const double   c[16], uint64_t Delta,
					 const uint64_t T1[2][L][N],
					 const uint64_t T2[2][L - 1][N],
					 const uint64_t T4[2][L - 2][N],
					 const uint64_t T8[2][L - 3][N],
					 const uint64_t evk_hat[DNUM][2][K*DNUM + K][N],
						   uint64_t res_hat[2][L - 4][N]) {
	uint64_t L3[2][2][L - 3][N];
	eval_poly_deg8<N, L, DNUM, K>(q, p, c    , Delta, T1, T2, T4, evk_hat, L3[0]);
	eval_poly_deg8<N, L, DNUM, K>(q, p, c + 8, Delta, T1, T2, T4, evk_hat, L3[1]);

	eval_poly_merge<N, L - 3, DNUM, K>(q, p, Delta, T8, evk_hat, L3[0], L3[1], res_hat);
}


//
template<int N, int L, int DNUM, int K>
void eval_poly_deg32(const uint64_t q[L],
					 const uint64_t p[K],
					 const double   c[32], uint64_t Delta,
					 const uint64_t T1 [2][L][N],
					 const uint64_t T2 [2][L - 1][N],
					 const uint64_t T4 [2][L - 2][N],
					 const uint64_t T8 [2][L - 3][N],
					 const uint64_t T16[2][L - 4][N],
					 const uint64_t evk_hat[DNUM][2][K*DNUM+K][N],
					       uint64_t res_hat[2][L - 5][N]) {
	uint64_t L4[2][2][L - 4][N];
	eval_poly_deg16<N, L, DNUM, K>(q, p, c   , Delta, T1, T2, T4, T8, evk_hat, L4[0]);
	eval_poly_deg16<N, L, DNUM, K>(q, p, c+16, Delta, T1, T2, T4, T8, evk_hat, L4[1]);

	eval_poly_merge<N, L - 4, DNUM, K>(q, p, Delta, T16, evk_hat, L4[0], L4[1], res_hat);
}


//
template<int N, int L, int DNUM, int K>
void eval_poly_deg64(const uint64_t q[L],
					 const uint64_t p[K],
					 const double   c[64], uint64_t Delta,
					 const uint64_t T1[2][L][N],
					 const uint64_t T2[2][L - 1][N],
					 const uint64_t T4[2][L - 2][N],
					 const uint64_t T8[2][L - 3][N],
					 const uint64_t T16[2][L - 4][N],
					 const uint64_t T32[2][L - 5][N],
					 const uint64_t evk_hat[DNUM][2][K*DNUM+K][N],
					       uint64_t res_hat[2][L - 6][N]) {
	uint64_t L5[2][2][L - 5][N];
	eval_poly_deg32<N, L, DNUM, K>(q, p, c     , Delta, T1, T2, T4, T8, T16, evk_hat, L5[0]);
	eval_poly_deg32<N, L, DNUM, K>(q, p, c + 32, Delta, T1, T2, T4, T8, T16, evk_hat, L5[1]);

	eval_poly_merge<N, L - 5, DNUM, K>(q, p, Delta, T32, evk_hat, L5[0], L5[1], res_hat);
}


//
template<int N, int L, int DNUM, int K>
void eval_poly_deg128(const uint64_t q[L],
					  const uint64_t p[K],
					  const double   c[128], uint64_t Delta,
					  const uint64_t T1[2][L][N],
					  const uint64_t T2[2][L - 1][N],
					  const uint64_t T4[2][L - 2][N],
					  const uint64_t T8[2][L - 3][N],
					  const uint64_t T16[2][L - 4][N],
					  const uint64_t T32[2][L - 5][N],
					  const uint64_t T64[2][L - 6][N],
					  const uint64_t evk_hat[DNUM][2][K*DNUM+K][N],
						    uint64_t res_hat[2][L - 7][N]) {
	uint64_t L6[2][2][L - 6][N];
	eval_poly_deg64<N, L, DNUM, K>(q, p, c     , Delta, T1, T2, T4, T8, T16, T32, evk_hat, L6[0]);
	eval_poly_deg64<N, L, DNUM, K>(q, p, c + 64, Delta, T1, T2, T4, T8, T16, T32, evk_hat, L6[1]);

	eval_poly_merge<N, L - 6, DNUM, K>(q, p, Delta, T64, evk_hat, L6[0], L6[1], res_hat);
}


//
template<int N, int L, int DNUM, int K>
void eval_poly_deg256(const uint64_t q[L],
					  const uint64_t p[K],
					  const double   c[256], uint64_t Delta,
					  const uint64_t T1  [2][L][N],
					  const uint64_t T2  [2][L - 1][N],
					  const uint64_t T4  [2][L - 2][N],
					  const uint64_t T8  [2][L - 3][N],
					  const uint64_t T16 [2][L - 4][N],
					  const uint64_t T32 [2][L - 5][N],
					  const uint64_t T64 [2][L - 6][N],
				      const uint64_t T128[2][L - 7][N],
					  const uint64_t evk_hat[DNUM][2][K*DNUM+K][N],
					        uint64_t res_hat[2][L - 8][N]) {
	uint64_t L7[2][2][L - 7][N];
	eval_poly_deg128<N, L, DNUM, K>(q, p, c      , Delta, T1, T2, T4, T8, T16, T32, T64, evk_hat, L7[0]);
	eval_poly_deg128<N, L, DNUM, K>(q, p, c + 128, Delta, T1, T2, T4, T8, T16, T32, T64, evk_hat, L7[1]);

	eval_poly_merge<N, L - 7, DNUM, K>(q, p, Delta, T128, evk_hat, L7[0], L7[1], res_hat);
}

//
template<int N>
void convert_poly_to_binarytreeform(double *c, int N) {
	if (N == 2) return;

	for (int i = 0; i < N / 2; i++) {
		c[i] -= c[N / 2 + i];
		c[N / 2 + i] *= 2;
	}

	convert_poly_to_binarytreeform<N>(c, N / 2);
	convert_poly_to_binarytreeform<N>(c + N / 2, N / 2);
}
