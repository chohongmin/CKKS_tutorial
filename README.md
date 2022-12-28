# CKKS_tutorial
This code provides a tutorial code for the CKKS-HEAAN crypto algorithm in RNS.
### BigInt
This folder contains an implementation of Big integers, self-contained and independent to any library. 128 bit multiplication should involve assembly routine, and it contains BigInt.asm, BigInt.h, BigInt.cpp. The assembly code is platform-dependent, should compile with Microsoft Macro assembler, that is built in usual VC++ compilers.

### DFT
DFT.h is an implementation of dft/idft and fft/ifft.
DFT matrix is a product of (logN-1) number of Sparsematrices, and FFT is the matrix multiplications of the sparse matrices.
The header includes SparseComplexMatrix and the decomposition of DFT matrix into FFT matrices.

### NumberTheory ( depends on BigInt )
Contains some basic routines in number theory such as modular multiplication and inversion, finding primitive roots, and NTT/INTT.

### Convert_poly_to_binarytreeForm 
Converts a polynomial in the Chebyshev form to the form of Binary-tree evaluation.
