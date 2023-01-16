# CKKS_tutorial
This code provides a tutorial code for the CKKS-HEAAN crypto algorithm in RNS. The provided code is not an optimized one but a word-to-word simple implementation just for tutorial purpose.

### Step 1. BigInt
This folder contains an implementation of Big integers, self-contained and independent to any library. 128 bit multiplication should involve assembly routine, and it contains BigInt.asm, BigInt.h, BigInt.cpp. The assembly code is platform-dependent, should compile with Microsoft Macro assembler, that is built in usual VC++ compilers.

### Step 2. DFT
DFT.h is an implementation of dft/idft and fft/ifft.
DFT matrix is a product of (logN-1) number of Sparsematrices, and FFT is the matrix multiplications of the sparse matrices.
The header includes SparseComplexMatrix and the decomposition of DFT matrix into FFT matrices.

### Step 3. NumberTheory 
Contains some basic routines in number theory such as modular multiplication and inversion, finding primitive roots, and NTT/INTT.

### Step 4. CKKS basic
Implements CKKS encodoing/decoding and encryption/decrytion together with Rescaling.

### Step 5. KS_with_DNUM
Key-switching algorithm with the gadget decomposition toolkit

### Step 6. Lineartransform
Implements the ciphertext operation that corresponds to a linear transform on data. The implementations involves rotation and key-switching.

### Step 7. CoeffToSlot
In this step, the CoeffToSlot algorithm in the case of logN=10 is implemented. FFT matrices are grouped by three matrices, of which rotation keys are generated. Then the CoeffToSlot is carried out.

### Step 8. EvalPoly
Converts a polynomial in the Chebyshev form to the form of Binary-tree evaluation, and evaluate the polynomial in the cases of degree=2,4,8,...,256.

### Step 9. EvalMod
Approximation of EvalMod by a polynomial of degree 128.

### Step 10. SlotToCoeff
The SlotToCoeff algorithm in the case of logN=10 is implemented. 

### All Combined
All the codes are combined to implement a basic bootstrapping algorithm.

