#pragma once

#include "CKKS_ks_wih_dnum.h"
#include "CKKS_lineartransform.h"

//
template< int N, int NumDiagsA, int NumDiagsB >
void matmul( const SparseComplexMatrix<N,NumDiagsA>& A,
	         const SparseComplexMatrix<N,NumDiagsB>& B,
	               SparseComplexMatrix<N,NumDiagsA*NumDiagsB>& C ){
	for(int iA=0;iA<NumDiagsA;iA++)
	for(int iB=0;iB<NumDiagsB;iB++){
		int shiftA = A.shift[iA];
		int shiftB = B.shift[iB];
		int iC=iA*NumDiagsB+iB;
		C.shift[iC]=(shiftA+shiftB)%N;
		for(int i=0;i<N;i++){
			int k=(i+shiftA)%N;
			C.diagr[iC][i]=A.diagr[iA][i]*B.diagr[iB][k]-A.diagi[iA][i]*B.diagi[iB][k];
			C.diagi[iC][i]=A.diagr[iA][i]*B.diagi[iB][k]+A.diagi[iA][i]*B.diagr[iB][k];
		}
	}
}

//
template<int N>
void conj(const int s[N], int sconj[N]){
	sconj[0]=s[0];
	for(int i=1;i<N;i++)
		sconj[N-i]=-s[i];
}
//
template< int N, int L >
void conj( const uint64_t q    [L],
	       const uint64_t a    [L][N],
	             uint64_t aconj[L][N]){
	for(int i=0;i<L;i++){
		aconj[i][0]=a[i][0];
		for(int j=1;j<N;j++)
			aconj[i][N-j]=(q[i]-a[i][j])%q[i];
	}
}
//
void splitU0R_logN_10( SparseComplexMatrix<1<<9,3*3*3> A[3]){
	SparseComplexMatrix<1<<9,3> E[9];
	splitU0R<10>(E);
	for(int d=0;d<3;d++){
		SparseComplexMatrix<1<<9,9> temp;
		matmul<1<<9,3,3>(E[3*d],E[3*d+1],temp);
		matmul<1<<9,9,3>(temp  ,E[3*d+2],A[d]);
	}
}

//
template <int N, int NumDiags>
void mat_conjtranspose( const SparseComplexMatrix<N,NumDiags>& A,
	                          SparseComplexMatrix<N,NumDiags>& Aconjtranspose){
	for(int d=0;d<NumDiags;d++){
		Aconjtranspose.shift[d]=A.shift[d];
		for(int i=0;i<N;i++){
			Aconjtranspose.diagr[d][i]=0;
			Aconjtranspose.diagi[d][i]=0;
		}
	}
	//
	for(int d=0;d<NumDiags;d++){
		int d_=0;
		for(;d_<NumDiags;d_++)
			if(Aconjtranspose.shift[d_]==(N-A.shift[d])%N)
				break;
		assert(d_!=NumDiags);
		//
		int s = A.shift[d];
		for(int i=0;i<N;i++){
			Aconjtranspose.diagr[d_][(i+s)%N] += A.diagr[d][i];
			Aconjtranspose.diagi[d_][(i+s)%N] -= A.diagi[d][i];
		}
	}
}

//
template< int L, int DNUM, int K>
void swkgen_logN_10( const uint64_t q[L],
	                 const uint64_t p[K], const int s[1<<10],
	                 const SparseComplexMatrix< (1<<9),27> A[3],
	                       uint64_t  evk_hat       [DNUM][2][DNUM*K+K][1<<10],
	                       uint64_t ckey_hat       [DNUM][2][DNUM*K+K][1<<10],
	                       uint64_t rkey_hat[3][27][DNUM][2][DNUM*K+K][1<<10]){
	int ss[1<<10]; conv<1<<10>(s,s,ss);
	swkgen<1<<10,L,DNUM>(ss,s,q,p,evk_hat);
	int sconj[1<<10]; conj<1<<10>(s,sconj);
	swkgen<1<<10,L,DNUM>(sconj,s,q,p,ckey_hat);
	for(int d=0;d<3;d++)
		rkey_gen<1<<10,L,DNUM,K,27>(A[d],q,p,s,rkey_hat[d]);
}

//
template< int L, int DNUM, int K>
void CoeffToSlot_logN_10( const uint64_t q[L],
	                      const uint64_t p[K], uint64_t Delta,
	                      const SparseComplexMatrix< 1<<9, 27> A[3],
						  const	uint64_t ckey_hat       [DNUM][2][DNUM*K+K][1<<10],
	                      const uint64_t rkey_hat[3][27][DNUM][2][DNUM*K+K][1<<10],
					      const uint64_t     ct_hat   [2][L  ][1<<10],
	                            uint64_t ct_cts_hat[2][2][L-3][1<<10]){
	const int N=1<<10;
	SparseComplexMatrix<1<<9,27> Aconjtranspose[3];
	for(int n=0;n<3;n++){
		mat_conjtranspose<1<<9,27>(A[n],Aconjtranspose[n]);
		for(int d=0;d<27;d++)
		for(int i=0;i<N/2;i++){
			Aconjtranspose[n].diagr[d][i] /= (n==0)?16:8;
			Aconjtranspose[n].diagi[d][i] /= (n==0)?16:8;
		}
	}
	//
	uint64_t ct1_hat[2][L][N];

	for(int i=0;i<2;i++)
	for(int j=0;j<L;j++)
	for(int k=0;k<N;k++) ct1_hat[i][j][k] = ct_hat[i][j][k];

	for(int n=0;n<3;n++){
		uint64_t temp[2][L][N];
		lineartransform<N,10,L,DNUM,K,27>(Aconjtranspose[n],Delta,q,p,rkey_hat[n],ct1_hat,temp);
		for(int i=0;i<2;i++)
		for(int j=0;j<L;j++)
		for(int k=0;k<N;k++) ct1_hat[i][j][k] = temp[i][j][k];
	}

	//
	uint64_t ct2_hat[2][L][N];
	{
		uint64_t temp[2][L][N];
		intt<N,L>(q,ct1_hat[0]);
		intt<N,L>(q,ct1_hat[1]);
		conj<N,L>(q,ct1_hat[0],temp[0]);
		conj<N,L>(q,ct1_hat[1],temp[1]);
		ntt<N,L>(q,ct1_hat[0]); ntt<N,L>(q,temp[0]);
		ntt<N,L>(q,ct1_hat[1]); ntt<N,L>(q,temp[1]);
		ks<N,L,DNUM,K>(q,p,ckey_hat,temp,ct2_hat);
	}

	//
	for(int i=0;i<2;i++)
	for(int j=0;j<L;j++)
	for(int k=0;k<N;k++)
		ct1_hat[i][j][k] = (ct1_hat[i][j][k]+ct2_hat[i][j][k])%q[j];

	RS_hat<1<<10,L,L-3>(q,ct1_hat,ct_cts_hat[0]);

	//
	uint64_t ptmi_hat[L][N];
	for(int i=0;i<L;i++){
		for(int j=0;j<N;j++)
			ptmi_hat[i][j]=0;
		ptmi_hat[i][N/2]=q[i]-1;
	}
	ntt<N,L>(q,ptmi_hat);
	//
	for(int i=0;i<2;i++)
	for(int j=0;j<L;j++)
	for(int k=0;k<N;k++)
		ct1_hat[i][j][k] = mul_mod(ptmi_hat[j][k],(ct1_hat[i][j][k]+
			               mul_mod(q[j]-2,ct2_hat[i][j][k],q[j]))%q[j],q[j]);

	RS_hat<1<<10,L,L-3>(q,ct1_hat,ct_cts_hat[1]);
}