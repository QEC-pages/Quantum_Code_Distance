# Quantum Code Distance

The repository contains a program written in GAP language for calculating the distance of a quantum code over an arbitrary field, using the matrix corresponding to the stabilizer group and exploiting a random algorithm.

The file "RandomDistance.gap" provides ReadGeneratorMatrix function for reading the matrix corresponding to the stabilizer group from the file, DoRandDist for calculating the distance of the code, SaveNewMTX and SaveAddMTX for creating a new file with parameters and matrix in Matrix Market format or inserting a code parameters into the existing file.

All additional details about the code functionality can be found in the last section of an overleaf project:

https://www.overleaf.com/read/fcyfjxwytsrj



Files "5qubitInteger.mtx" and "5qubitComplex.mtx" are test files with matricies corresponding to 5 qubit code stored in a full and compressed forms.
