# Quantum Code Distance

The repository contains a program written in GAP language for calculating distance of a quantum code over arbitrary field, using its stabilizer matrix and exploiting a random algorithm.

The file "RandomDistance.gap" provides ReadGeneratorMatrix function for reading a matrix of stabilizer group from file, DoRandDist for calculating distance of the code, SaveNewMTX and SaveAddMTX for creating a new file with matrix in matrix market or inserting a code parameters into existing file.

All additional details about the code functionality can be found in the last section of an overleaf project:
https://www.overleaf.com/read/fcyfjxwytsrj



Files "5qubitInteger.mtx" and "5qubitComplex.mtx" are test files with matricies for 5 qubit code stored in a full and compressed forms.
