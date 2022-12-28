# CKKS_tutorial
This code provides a tutorial code for the CKKS-HEAAN crypto algorithm in RNS.
### BigInt
This folder contains an implementation of Big integers, self-contained and independent to any library. 128 bit multiplication should involve assembly routine, and it contains BigInt.asm, BigInt.h, BigInt.cpp. The assembly code is platform-dependent, should compile with Microsoft Macro assembler, that is built in usual VC++ compilers.

### DFT
DFT.h is an implementation of dft/idft and fft/ifft.
DFT matrix is a product of (logN-1) number of Sparsematrices, and FFT is the matrix multiplications of the sparse matrices.
The header includes SparseComplexMatrix and the decomposition of DFT matrix into FFT matrices.

### Convert_poly_to_binarytreeForm 
<img width="512" alt="readme_convert_poly" src="https://user-images.githubusercontent.com/121416455/209746234-60ee7df1-adf1-410c-aa72-86d9ec1134ea.png">
