#pragma once

#include <stdint.h>
#include <inttypes.h>
#include <unordered_set>
#include <chrono>
#include <iostream>
#include "BigInt.h"

uint64_t   mul_mod( uint64_t a, uint64_t b, uint64_t q);
uint64_t power_mod( uint64_t a, uint64_t n, uint64_t q );
uint64_t   inv_mod( uint64_t a, uint64_t q );
uint64_t rand_uint64();
bool isprime(uint64_t n);
void findPrimeFactors( uint64_t n, std::unordered_set<uint64_t>& s );
uint64_t findPrimitive(uint64_t n);

//
template<int L>
void find_RNS_primes( uint64_t Delta, int N, uint64_t q[L]){
	if( Delta==(1ULL<<50) && N<=(1<<17) && L<=30 ){
		uint64_t q_save[30]={1125899915231233ULL,1125899927027713ULL,1125899930959873ULL,1125899935940609ULL,1125899941445633ULL,1125899963990017ULL,
							 1125899976572929ULL,1125899978145793ULL,1125899989155841ULL,1125899996233729ULL,1125899996495873ULL,1125899998068737ULL,
							 1125899998855169ULL,1125900015108097ULL,1125900015370241ULL,1125900027166721ULL,1125900043419649ULL,1125900047351809ULL,
							 1125900052856833ULL,1125900053905409ULL,1125900055478273ULL,1125900063866881ULL,1125900070944769ULL,1125900076711937ULL,
							 1125900078284801ULL,1125900079595521ULL,1125900083003393ULL,1125900093751297ULL,1125900094013441ULL,1125900097159169ULL};
		for(int i=0;i<L;i++)
			q[i]=q_save[i];
	}
	else if( Delta==(1ULL<<60) && N<=(1<<17) && L<=30 ){
		uint64_t q_save[30]={1152921504616808449ULL,1152921504618381313ULL,1152921504622575617ULL,1152921504634109953ULL,1152921504643809281ULL,
							 1152921504650100737ULL,1152921504663994369ULL,1152921504666615809ULL,1152921504672645121ULL,1152921504687849473ULL,
							 1152921504690208769ULL,1152921504690995201ULL,1152921504694927361ULL,1152921504695451649ULL,1152921504695713793ULL,
							 1152921504697286657ULL,1152921504717733889ULL,1152921504722190337ULL,1152921504730316801ULL,1152921504730841089ULL,
							 1152921504738181121ULL,1152921504741851137ULL,1152921504746569729ULL,1152921504756006913ULL,1152921504760201217ULL,
							 1152921504760987649ULL,1152921504767016961ULL,1152921504770424833ULL,1152921504783794177ULL,1152921504786677761ULL};
		for(int i=0;i<L;i++)
			q[i]=q_save[i];
	}
	else{
		int count=0;
		uint64_t a=Delta;
		while(count<L){
			a++;
			if((a-1)%(2*N)==0 && isprime(a)){
				q[count]=a;
				count++;
			}
		}
	}
}

//-------------------------------------------------------------------------------
// bitReverse
//-------------------------------------------------------------------------------
template< int N >
void bitReverse( const uint64_t a[N], uint64_t b[N] ){
	int logN=0;
	while(N>(1<<logN))
		logN++;
	for(int i=0;i<N;i++){
		int j=0;
		for(int k=0;k<logN;k++)
			j+=((i>>k)&1)<<(logN-1-k);
		b[j]=a[i];
	}
}

//
template<int N>
void ntt(uint64_t a[N], uint64_t p ){
	uint64_t psi=findPrimitive(p)%p;
	psi = power_mod(psi,(p-1)/(2*N),p);
	uint64_t Psi[N], PsiRev[N]; Psi[0]=1;
	for(int i=1;i<N;i++)
		Psi[i]=mul_mod(psi,Psi[i-1],p);
	bitReverse<N>(Psi,PsiRev);
	for(int m=1,t=N/2;m<N;m*=2,t/=2){
		for(int i=0;i<m;i++){
			uint64_t* ae=a+2*i*t;
			uint64_t* ao=ae+t;
			for(int k=0;k<t;k++){
				uint64_t U=mul_mod(ao[k],PsiRev[m+i],p);
				ao[k]=(ae[k]+p-U)%p;
				ae[k]=(ae[k]  +U)%p;
			}
		}
	}
}

//
template<int N>
void intt(uint64_t a[N], uint64_t p ){
	uint64_t psi=findPrimitive(p)%p; 
	psi = power_mod(psi,(p-1)/(2*N),p); 
	psi=inv_mod(psi,p);
	uint64_t Psi[N], PsiRev[N]; Psi[0]=1;
	for(int i=1;i<N;i++)
		Psi[i]=mul_mod(psi,Psi[i-1],p);
	bitReverse<N>(Psi,PsiRev);
	for(int t=1,m=N/2;m>0;t*=2,m/=2){
		for(int i=0;i<m;i++){
			uint64_t* ae=a+2*i*t;
			uint64_t* ao=ae+t;
			for(int k=0;k<t;k++){
				uint64_t U=(ae[k]+p-ao[k])%p;
				ae[k]=(ae[k]+ao[k])%p;
				ao[k]=mul_mod(U,PsiRev[m+i],p);
			}
		}
	}
	uint64_t Ninv=inv_mod(N,p);
	for(int i=0;i<N;i++) a[i]=mul_mod(a[i],Ninv,p);
}

//
template<int L> void icrt(const uint64_t q[L],
	                      const uint64_t a[L], BigInt& A ){
	BigInt Q(L); Q[0]=1;
	for(int i=0;i<L;i++){
		BigInt temp(Q); temp.mul(q[i],Q);}
	//Q.print("Q");
	A.resize(L+1);
	for(int i=0;i<L;i++){
		BigInt Qoverq; uint64_t r;
		Q.div(q[i],Qoverq,r);
		//Qoverq.print("Qoverq");
		BigInt temp1; 
		Qoverq.mul( mul_mod(a[i],inv_mod(Qoverq.mod(q[i]),q[i]),q[i]),temp1);
		BigInt temp2(A); temp1.add(temp2,A);
	}
	//A.print("A");
	while(A>=Q  ){
		BigInt temp; 
		A.sub(Q  ,temp); 
		A=temp;
	}
}

//
template<int N, int L> 
void ntt(const uint64_t q[L],
		 const uint64_t p[L], uint64_t a[2*L][N]){
	for(int i=0; i<L; i++){
		ntt<N>(a[i], q[i]);
		ntt<N>(a[L+i], p[i]);
	}
}

//
template<int N, int L> 
void ntt(const uint64_t q[L], uint64_t a[L][N]){
	for(int i=0; i<L; i++) ntt<N>(a[i], q[i]);
}

//
template<int N, int L> 
void intt(const uint64_t q[L],
		  const uint64_t p[L], uint64_t a[2*L][N]){
	for(int i=0; i<L; i++){
		intt<N>(a[i], q[i]);
		intt<N>(a[L+i], p[i]);
	}
}

//
template<int N, int L> 
void intt(const uint64_t q[L], uint64_t a[L][N]){
	for(int i=0; i<L; i++) intt<N>(a[i], q[i]);
}


//
template<int N, int L>
void rot(const uint64_t q[L], const uint64_t a[L][N],
									uint64_t a_rot[L][N]){
	for(int i=0; i<L; i++)
	for(int j=0; j<N; j++){
		int Aj = (5*j) / N;
		int rj = (5*j) % N;

		a_rot[i][rj] = (Aj % 2 == 0) ? a[i][j] : ((q[i]-a[i][j])%q[i]);
	}
}

//
template<int N, int L>
void rot_hat(const uint64_t a_hat[L][N], uint64_t a_hat_rot[L][N]){
	uint64_t Normal[N]; 
	uint64_t Reverse[N];

	for(int i=0; i<N; i++) Normal[i]=i;

	bitReverse<N>(Normal, Reverse);

	for(int j=0; j<N; j++){
		int bRj = Reverse[j];
		int rj = (10*bRj+5)%(2*N);
		for(int i=0; i<L; i++)
			a_hat_rot[i][j] = a_hat[i][Reverse[(rj-1)/2]];
	}
}

// ADDED by KSH 1024
// rot pt by r
template<int N, int L>
void rot(const uint64_t q[L], const uint64_t a[L][N],
									uint64_t a_rot[L][N],
									const int r){
	// pow = (5^r) % (2*N)
	int pow = 1;
	for(int i = 0; i < r; ++i) {
		pow = (5 * pow) % (2*N);
	}

	for(int l = 0; l < L; l++)
	for(int i = 0, j = 0; i < N; ++i, j = (j + pow) % (2*N)) {
		if(j < N)
			a_rot[l][j] = a[l][i];
		else
			a_rot[l][j-N] = (q[l] - a[l][i])%q[l];
	}
}

// ADDED by KSH 1024
// rot_hat pt by r
template<int N, int L>
void rot_hat(const uint64_t a_hat[L][N], uint64_t a_hat_rot[L][N], const int r){
	uint64_t Normal[N]; 
	uint64_t Reverse[N];

	for(int i=0; i<N; i++) Normal[i]=i;

	bitReverse<N>(Normal, Reverse);

	// pow = (5^r) % (2*N)
	int pow = 1;
	for(int i = 0; i < r; ++i) {
		pow = (5 * pow) % (2*N);
	}

	for(int j=0; j<N; j++){
		int bRj = Reverse[j];
		int rj_orig = (2*bRj+1)%(2*N);
		int rj = ((uint64_t) rj_orig * (uint64_t) (pow)) % (2*N);
		for(int i=0; i<L; i++)
			a_hat_rot[i][j] = a_hat[i][Reverse[(rj-1)/2]];
	}
}

// ADDED by KSH 1024
// rot_hat ct by r
template<int N, int L>
void rot_hat(
	const uint64_t ct_hat[2][L][N],
	uint64_t ct_rot_hat[2][L][N],
	const int r) {
	rot_hat<N,L>(ct_hat[0], ct_rot_hat[0], r);
	rot_hat<N,L>(ct_hat[1], ct_rot_hat[1], r);
}

// ADDED by KSH 1024
// rot s by r
template<int N>
void rot( const int s[N], int s_rot[N], const int r){
	// pow = (5^r) % (2*N)
	int pow = 1;
	for(int i = 0; i < r; ++i) {
		pow = (5 * pow) % (2*N);
	}

	for(int i = 0, j = 0; i < N; ++i, j = (j + pow) % (2*N)) {
		if(j < N)
			s_rot[j] = s[i];
		else
			s_rot[j-N] = -s[i];
	}
}

// ADDED by KSH 1024
// conj s
template<int N>
void conj( const int s[N], int s_conj[N]){
	s_conj[0] = s[0];
	for(int i=1; i<N; i++) s_conj[i] = -s[N-i];
}

// ADDED by KSH 1024
// conj pt
template<int N, int L>
void conj(
	const uint64_t q[L],
	const uint64_t pt[L][N],
	uint64_t pt_conj[L][N]){
	for(int i = 0; i < L; ++i) {
		pt_conj[i][0] = pt[i][0];
		for(int j=1; j<N; j++){
			pt_conj[i][j] = q[i] - pt[i][N-j];
		}
	}
}

// ADDED by KSH 1024
// conj_hat pt
template<int N, int L>
void conj_hat(
	const uint64_t q[L],
	const uint64_t pt_hat[L][N],
			uint64_t pt_conj_hat[L][N]){
	uint64_t Normal[N]; 
	uint64_t Reverse[N];

	for(int i=0; i<N; i++) Normal[i]=i;

	bitReverse<N>(Normal, Reverse);

	for(int j=0; j<N; j++){
		int bRj = Reverse[j];
		int rj_orig = (2*bRj+1)%(2*N);
		int rj = 2*N - rj_orig;
		for(int i=0; i<L; i++) {
			pt_conj_hat[i][j] = pt_hat[i][Reverse[(rj-1)/2]];
		}
	}
}

// ADDED by KSH 1024
// conj_hat ct
template<int N, int L>
void conj_hat(const uint64_t q[L],
		  const uint64_t ct[2][L][N],
			uint64_t ct_conj[2][L][N]){
	conj_hat<N,L>(q, ct[0], ct_conj[0]);
	conj_hat<N,L>(q, ct[1], ct_conj[1]);
}







