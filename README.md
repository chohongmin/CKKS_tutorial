# CKKS_tutorial
This code provides a tutorial code for the CKKS-HEAAN crypto algorithm in RNS.

Folder "Convert_poly_to_binarytreeForm" 
: it containts a code that converts polynomials in the form sum(c[i]*T_i(x),i,0,N-1) to binary-tree evlauation form 
sum(c[i]*T_i(x),i,0,N/2-1)+T_{N/2}*sum(c[N/2+i]*T_i(x),i,0,N/2-1)
