#pragma once

#include "CKKS_basic.h"
/////////////////////////////////////////////////////////////////////////////
// Key Switch with DNUM
/////////////////////////////////////////////////////////////////////////////

template<int N, int L, int K>
void modUP(const uint64_t q[L],
		   const uint64_t p[K],
		   const uint64_t a[L][N],
		         uint64_t a_tilde[L+K][N]) {
	uint64_t Qmod[L][K];
	for(int j=0; j<L; j++)
	for(int k=0; k<K; k++) {
		Qmod[j][k] = 1;
		for (int i = 0; i < L; i++)
			if (i != j) 
				Qmod[j][k] = mul_mod(Qmod[j][k], q[i] % p[k], p[k]);
	}

	uint64_t invQmod[L];
	for (int j = 0; j < L; j++) {
		invQmod[j] = 1;
		for (int i = 0; i < L; i++)
			if (i != j)
				invQmod[j] = mul_mod(invQmod[j], inv_mod(q[i], q[j]), q[j]);
	}

	uint64_t QmodP[K];
	for (int k = 0; k < K; k++) {
		QmodP[k] = 1;
		for (int j = 0; j < L; j++)
			QmodP[k] = mul_mod(QmodP[k], q[j] % p[k], p[k]);
	}

	for (int i = 0; i < N; i++) {
		uint64_t b[L];
		int count = 0;
		for (int j = 0; j < L; j++) {
			b[j] = mul_mod(a[j][i], invQmod[j], q[j]);
			if (2 * b[j] >= q[j]) count++;
		}

		for (int k = 0; k < L; k++)
			a_tilde[k][i] = a[k][i];

		for (int k = 0; k < K; k++) {
			a_tilde[L + k][i] = 0;
			for (int j = 0; j < L; j++)
				a_tilde[L + k][i] += mul_mod(b[j] % p[k], Qmod[j][k], p[k]);

			if (count > 0)
				a_tilde[L + k][i] += mul_mod(QmodP[k], p[k] - count, p[k]);
			a_tilde[L + k][i] = a_tilde[L + k][i] % p[k];
		}
	}
}

//
template<int L, int DNUM>
void gadget_g(const uint64_t q[L],
			  const uint64_t p[L / DNUM],
					uint64_t g[DNUM][L + (L / DNUM)]){
	const int K = L / DNUM; assert(L % DNUM == 0);

	for(int d = 0; d < DNUM; d++)
		for(int i = 0; i < L + K; i++){
			if(d*K <= i && i < (d+1)*K) g[d][i] = 1;
			else g[d][i] = 0;
		}
}


//
template<int N, int L, int DNUM>
void gadget_ginv(const uint64_t q[L],
				 const uint64_t p[L / DNUM],
				 const uint64_t a[L][N],
					   uint64_t ginva[DNUM][L + (L / DNUM)][N]) {
	const int K = L / DNUM;
	
	uint64_t QP[L+K];
	for(int i=0; i<L; i++) QP[i] = q[i];
	for(int i=0; i<K; i++) QP[L+i] = p[i];

	for(int d=0; d<DNUM; d++){
		for(int i=0; i<K; i++){
			QP[i] = QP[d*K+i];
			QP[d*K+i] = q[i];
		}

		uint64_t a_Dd[K][N];
		for(int i=0; i<K; i++)
		for(int j=0; j<N; j++)
			a_Dd[i][j] = a[d*K+i][j];

		modUP<N, K, L>(QP, QP+K, a_Dd, ginva[d]);

		for(int i=0; i<K; i++)
		for(int j=0; j<N; j++){
			uint64_t temp = ginva[d][i][j];
			ginva[d][i][j] = ginva[d][d*K+i][j];
			ginva[d][d*K+i][j] = temp;
		}

		for(int i=0; i<K; i++){
			uint64_t temp = QP[i];
			QP[i] = QP[d*K+i];
			QP[d*K+i] = temp;
		}

	}
}

//
template<int N, int L, int DNUM>
void swkgen(const int sfr[N],
			const int sto[N],
			const uint64_t q[L],
			const uint64_t p[L / DNUM],
				  uint64_t swk_hat[DNUM][2][L + (L / DNUM)][N]) {
	const int K = L / DNUM;
	uint64_t g[DNUM][L + K];
	gadget_g<L, DNUM>(q, p, g);

	for (int n = 0; n < DNUM; n++) {
		uint64_t pt[L + K][N];
		for (int j = 0; j < L; j++) {
			uint64_t P = 1;
			for (int k = 0; k < K; k++)
				P = mul_mod(P, p[k] % q[j], q[j]);
			uint64_t Pg = mul_mod(P, g[n][j], q[j]);
			for (int i = 0; i < N; i++)
				pt[j][i] = mul_mod((q[j] + sfr[i]) % q[j], Pg, q[j]);
		}
		for (int j = 0; j < K; j++)
		for (int i = 0; i < N; i++)
			pt[j + L][i] = 0;
		uint64_t qp[L + K];
		for (int j = 0; j < L; j++) qp[    j] = q[j];
		for (int j = 0; j < K; j++) qp[L + j] = p[j];

		enc<N, L + K>(pt, sto, qp, swk_hat[n]);
	}
}

//
template<int N, int L, int DNUM>
void ks(const uint64_t q[L],
		const uint64_t p[L / DNUM],
		const uint64_t swk_hat[DNUM][2][L + (L / DNUM)][N],
		const uint64_t ct_hat [2][L][N],
			  uint64_t out_hat[2][L][N]) {
	uint64_t a[L][N];
	for (int i = 0; i < L; i++)
	for (int j = 0; j < N; j++)
		a[i][j] = ct_hat[1][i][j];
	intt<N, L>(q, a);
	
	const int K = L / DNUM;
	uint64_t ginva[DNUM][L + K][N];
	gadget_ginv<N, L, DNUM>(q, p, a, ginva);
	
	uint64_t qp[L + K];
	for (int j = 0; j < L; j++) qp[  j] = q[j];
	for (int j = 0; j < K; j++) qp[L+j] = p[j];
	
	uint64_t temp   [2][L + K][N];
	uint64_t temp_rs[2][L    ][N];
	for (int d = 0; d < DNUM; d++) {
		ntt<N, L + K>(qp, ginva[d]);
		for(int i=0; i<L+K; i++)
		for(int j=0; j<  N; j++) {
			temp[0][i][j] = mul_mod(ginva[d][i][j], swk_hat[d][0][i][j], qp[i]);
			temp[1][i][j] = mul_mod(ginva[d][i][j], swk_hat[d][1][i][j], qp[i]);
		}

		RS_hat<N, L + K, L>(qp, temp, temp_rs);
		for (int i = 0; i < L; i++)
		for (int j = 0; j < N; j++){
			if (d == 0) {
				out_hat[0][i][j] = temp_rs[0][i][j];
				out_hat[1][i][j] = temp_rs[1][i][j];
			}
			else {
				out_hat[0][i][j] = (out_hat[0][i][j] + temp_rs[0][i][j]) % q[i];
				out_hat[1][i][j] = (out_hat[1][i][j] + temp_rs[1][i][j]) % q[i];
			}
		}
	}

	for (int i = 0; i < L; i++)
	for (int j = 0; j < N; j++)
		out_hat[0][i][j] = (out_hat[0][i][j] + ct_hat[0][i][j]) % q[i];
}
/////////////////////////////////////////////////////////////////////////////