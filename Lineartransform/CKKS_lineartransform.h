#pragma once

#include "CKKS_basic.h"
#include "CKKS_ks_wih_dnum.h"

template< int N, int L>
void rot( const uint64_t q[L], const uint64_t m   [L][N], int r,
	                                 uint64_t mrot[L][N]){
	if(r==0){
		for(int i=0;i<L;i++)
		for(int j=0;j<N;j++)
			mrot[i][j]=m[i][j];
	}
	else{
		uint64_t temp[L][N];
		for(int i=0;i<L;i++)
		for(int j=0;j<N;j++){
			uint64_t Aj=(5*j)/N;
			uint64_t rj=(5*j)%N;
			temp[i][rj]=(Aj%2==0)?m[i][j]:((q[i]-m[i][j])%q[i]);
		}
		rot<N,L>(q, temp,r-1,mrot);
	}
}
//
template<int N>
void rot( const int s[N], int r, int srot[N] ){
	if(r==0){
		for(int j=0;j<N;j++)
			srot[j]=s[j];
	}
	else{
		int temp[N];
		for(int j=0;j<N;j++){
			int Aj=(5*j)/N;
			int rj=(5*j)%N;
			temp[rj]=(Aj%2==0)?s[j]:-s[j];
		}
		rot<N>(temp,r-1,srot);
	}
}
//
template< int N, int L, int DNUM, int K, int NumDiags >
void rkey_gen( const SparseComplexMatrix<N/2,NumDiags>& A,
	           const uint64_t q[L],
	           const uint64_t p[K],
	           const  int     s[N],
	                 uint64_t rkey_hat[NumDiags][DNUM][2][DNUM*K+K][N]){
	for(int i=0;i<NumDiags;i++){
		int srot[N];
		rot<N>(s,A.shift[i],srot);
		swkgen<N,L,DNUM>(srot,s,q,p,rkey_hat[i]);
	}
}

//
template<int N, int logN, int L, int DNUM, int K, int NumDiags>
void lineartransform(const SparseComplexMatrix<N/2, NumDiags> & A,
				     const uint64_t Delta,
					 const uint64_t q[L],
					 const uint64_t p[K],
					 const uint64_t rkey_hat[NumDiags][DNUM][2][DNUM*K+K][N],
					 const uint64_t ct_hat[2][L][N],
						   uint64_t res_hat[2][L][N]){
	uint64_t ct[2][L][N];
	for(int i=0; i<2; i++)
	for(int j=0; j<L; j++)
	for(int k=0; k<N; k++){
		ct[i][j][k] = ct_hat[i][j][k];
		res_hat[i][j][k] = 0;
	}
	intt<N, L>(q, ct[0]);
	intt<N, L>(q, ct[1]);

	for(int d=0; d<NumDiags; d++){
		uint64_t pt[L][N];
		encode<N, logN, L>(A.diagr[d], A.diagi[d], Delta, q, pt);

		uint64_t temp[2][L][N];
		rot<N, L>(q, ct[0], A.shift[d], temp[0]);
		ntt<N, L>(q, temp[0]);

		rot<N, L>(q, ct[1], A.shift[d], temp[1]);
		ntt<N, L>(q, temp[1]);

		uint64_t ct_rot_hat[2][L][N];
		ks<N, L, DNUM, K>(q, p, rkey_hat[d], temp, ct_rot_hat);

		ntt<N, L>(q, pt);
		for(int i=0; i<2; i++)
		for(int j=0; j<L; j++)
		for(int k=0; k<N; k++)
			res_hat[i][j][k] = (res_hat[i][j][k] + mul_mod(pt[j][k], ct_rot_hat[i][j][k], q[j]))%q[j];
	}
}
	


