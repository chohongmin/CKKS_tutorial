# CKKS_tutorial
This code provides a tutorial code for the CKKS-HEAAN crypto algorithm in RNS.
### 1. BigInt
This folder contains an implementation of Big integers, self-contained and independent to any library. 128 bit multiplication should involve assembly routine, and it contains BigInt.asm, BigInt.h, BigInt.cpp. The assembly code is platform-dependent, should compile with Microsoft Macro assembler, that is built in usual VC++ compilers.

### 2. DFT
DFT.h is an implementation of dft/idft and fft/ifft.
DFT matrix is a product of (logN-1) number of Sparsematrices, and FFT is the matrix multiplications of the sparse matrices.
The header includes SparseComplexMatrix and the decomposition of DFT matrix into FFT matrices.

### 3. NumberTheory 
Contains some basic routines in number theory such as modular multiplication and inversion, finding primitive roots, and NTT/INTT.

### 4. CKKS basic
Implements CKKS encodoing/decoding and encryption/decrytion together with Rescaling.

### 5. KS_with_DNUM
Key-switching algorithm with the gadget decomposition toolkit

### 6. Lineartransform
Implements the ciphertext operation that corresponds to a linear transform on data. The implementations involves rotation and key-switching.

### Convert_poly_to_binarytreeForm 
Converts a polynomial in the Chebyshev form to the form of Binary-tree evaluation.
