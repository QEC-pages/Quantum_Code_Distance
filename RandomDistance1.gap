# Functions for finding upper bound on the distance of a quantum code over field F


# Package for some function, used in code
LoadPackage("guava");
LoadPackage("float");


# Reading generator matrix G from file. It has to be in a specific MatrixFormat with permuted columns (details in pdf file).
ReadGeneratorMatrix := function(StrPath)
    # StrPath - path to the file with matrix

    local input, data, pair, F, rowsG, colsG, G, i, iComment;
    # local variables:
    # input - text, read from file and converted to string
    # data - string of numbers converted to list of them
    # pair - parameter, which can be either "integer" (normal storage) or "complex" (compress storage)
    # F - field, over which we are working
    # rowsG - number of rows in G
    # colsG - number of columns in G
    # G - generator matrix
    # i - for "for" loop
    # pair - indicates if we store matrix in the compressed form a+i*b (normal form if pair=integer and compressed form if pair=complex),
    # iComment - line, where comment section ends

    # Read from file
    input := ReadAll(InputTextFile(StrPath));;
    # Splitting lines
    data := SplitString(input, "\n");;

    # Splitting string in a line
    for i in [1..Length(data)] do
	data[i] := SplitString(data[i], " ");;
    od;

    # Find pair
    pair := data[1, 5];;

    # Find field
    F := EvalString(data[2, 3]);;

    # We search where top comment section ends:
    iComment := 2;
    while (data[iComment + 1, 1] = "%") do
	iComment := iComment + 1;;
    od;

    # Parameters (dimensions and number of non-zero elemments) of G
    rowsG := Int(data[iComment + 1, 1]);;
    colsG := Int(data[iComment + 1, 2]);;

    # We intoruce empty matrix
    G := NullMat(rowsG, colsG, F);;

    # Then we fill G with the elements from data
    if (pair = "integer") then
	if IsPrime(Size(F)) then
	    for i in [(iComment + 2)..Length(data)] do
		if (data[i, 1] <> "%") then
		    if (Int(data[i, 3]) = -1) then
			G[Int(data[i, 1]), Int(data[i, 2])] := Zero(F);;
		    else
			G[Int(data[i, 1]), Int(data[i, 2])] := Int(data[i, 3]) * One(F);;
		    fi;
		fi;
	    od;
	else
	    for i in [(iComment + 2)..Length(data)] do
                if (data[i, 1] <> "%") then
                    if (Int(data[i, 3]) = -1) then
                        G[Int(data[i, 1]), Int(data[i, 2])] := Zero(F);;
                    else
			G[Int(data[i, 1]), Int(data[i, 2])] := PrimitiveElement(F) ^ Int(data[i, 3]);;
                    fi;
                fi;
            od;
        fi;
    elif (pair = "complex") then
	if IsPrime(Size(F)) then
            for i in [(iComment + 2)..Length(data)] do
                if (data[i, 1] <> "%") then
                    if (Int(data[i, 3]) = -1) then
                        G[Int(data[i, 1]), 2 * Int(data[i, 2]) - 1] := Zero(F);;
                    else
			G[Int(data[i, 1]), 2 * Int(data[i, 2]) - 1] := Int(data[i, 3]) * One(F);;
                    fi;

		    if (Int(data[i, 4]) = -1) then
                        G[Int(data[i, 1]), 2 * Int(data[i, 2])] := Zero(F);;
                    else
                        G[Int(data[i, 1]), 2 * Int(data[i, 2])] := Int(data[i, 4]) * One(F);;
                    fi;
                fi;
            od;
        else
            for i in [(iComment + 2)..Length(data)] do
                if (data[i, 1] <> "%") then
                    if (Int(data[i, 3]) = -1) then
                        G[Int(data[i, 1]), 2 * Int(data[i, 2]) - 1] := Zero(F);;
                    else
                        G[Int(data[i, 1]), 2 * Int(data[i, 2]) - 1] := PrimitiveElement(F) ^ Int(data[i, 3]);;
                    fi;

		    if (Int(data[i, 4]) = -1) then
                        G[Int(data[i, 1]), 2 * Int(data[i, 2])] := Zero(F);;
                    else
                        G[Int(data[i, 1]), 2 * Int(data[i, 2])] := PrimitiveElement(F) ^ Int(data[i, 4]);;
                    fi;
                fi;
            od;
        fi;
    else
	Error("\n \n", "The pair != \"integer\" or \"complex\"!", "\n \n");
    fi;

    return G;
end;


# Function: Symplectic weight of a vector
#
# Output: symplectic weight of a vector
# Input: vec - vector in the form (a1,b1,a2,b2,...)
SymplVecWeight := function(vec, F)

    local sympvecweight, i, len; # local variables: sympvecweight - the symplectic weight, i - for "for" loop, len - length of vec

    len := Length(vec);;
    sympvecweight := 0;;

    for i in [1..(len/2)] do
	# for each odd and even pairs (a1,b1), (a2,b2), etc we check if ai or bi is nonzero. If it is, then operator acts nontrivially on
	# a qubit i. Thus we increase symplectic weight by 1
	if (vec[2*i-1] <> Zero(F) or vec[2*i] <> Zero(F)) then
            sympvecweight := sympvecweight + 1;;
        fi;
    od;

    return sympvecweight;
end;


# Function: Check matrix introducing.
# We assume G is in the form (A1,B1,A2,B2,...), then we want H to has the form (-B1,A1,-B2,A2,...).
# We operate with transposed G, since it is much easier to take and change a particular row, then a particular column.
CheckMatrix := function(G, F, debug)

    local dims, i, H; # local variables: dims - dimensions of G matrix, i - for "for" loop, H - check matrix 

    dims := DimensionsMat(G);;

    # Checking that G has even number of columns
    if (not IsInt(dims[2] / 2)) then
	Error("\n \n", "Generator Matrix G has odd number of columns!", "\n \n");
    fi;

    # Introducing check matrix
    H := TransposedMatMutable(G);;
    for i in [1..(dims[2]/2)] do
        H[2*i] := (-1) * H[2*i];; #H = (A1,-B1,A2,-B2,...)^T
        H := Permuted(H, (2*i-1,2*i));; # H = (-B1,A1,-B2,A2,...)^T. Permutation of rows in transposed matrix (so, permutation of columns in initial).
    od;

    # Checking the result on ortogonality GH^T = 0, if 2-nd bit of debug equals 1. Right now we have H^T stored in H. In .tex file I showed why this condition can be checked using the same form GH^T,
    # even if the columns are switched. NullMat - zero matrix.
    if (debug[2] = 1) then
	if (G * H = NullMat(dims[1], dims[1], F)) then
	    Print("Ortogonality GH^T condition is met!", "\n");
	else
	    Error("\n \n", "Problem with ortogonality GH^T!", "\n \n");
	fi;
    fi;

    # We return to the normal form: H = (-B1,A1,-B2,A2,...)
    H := TransposedMatMutable(H);;

    return H;
end;


# Probability calculation
AverageCalc := function(Multiplicities)
    return Sum(Multiplicities)/Length(Multiplicities);
end;


# Entropy function
EntrFunc := function(x)
    return (-x*Log2(x) - (1-x)*Log2(1-x));
end;


# Function: a simplified version of a Random Information Set algorithm. See tex file for more details about concepts.
#
# Input:
# G - input stabilizer matrix. It has to consist of two blocks (A|B). We represent them as G=(A1,B1,A2,B2,...), where Ai and Bi - columns of A and B respectively.
# num - number of information sets to construct. In the algorithm we have to many times create random permutation and then look at the weights. There are several terminators to stop this process:
# we get too small (e.g. uninteresting) distance, we make too much iterations, one particular weight appears large enough number of times, and many others. num controls the number of iterations.
# It should be exponentially large (num ~ exp(n)).
# mindist - the algorithm stops when distance below its value is found (set it to 1 if you want actual distance).
# F - field to use (e.g., GF(9); it has to match the field of the matrices).
# debug - parameter for additional output. Its binary representation will be found:
# bit = 0 - no debugging, bit = 1 - debugging.
# 1-st bit corresponds to the displaying of additional information in the end (program prints first and last found vectors corresponding to the found upper bound on distance,
# if their lengths<=100; prints how many times the first vector corresponding to the upper bound on the distance was found; prints overall number of found vectors corresponding to 
# found upper bound on the distance; prints found upper bound on the distance; prints number of iterations made before the completion of the program;
# prints the probability of existence of a vector with lower weight).
# 2-nd bit corresponds to the checking of ortogonality GH^T = 0 in CheckMatrix function.
# 3-rd bit corresponds to the displaying of the upper bound on the distance found so far and the overall number of times vectors were found with this weight. It prints this information 
# if the the number of iteration rounds is 200*i times for some integer i.
# 4-th bit corresponds to the ... !!!
# Output: parameters of the code as a list of strings.
DoRandDist := function(arg)
    # The function receives arg which is a list [G, num, mindist, F, debug] or [G, num, mindist, F]. This method of defining the function allows us to get changeable number of variables, which is
    # useful, since we do not want to pass debug parameter everytime to the function. If we do not pass debug parameter, we assume it to be 0 (no debugging).

    local G, num, mindist, F, debug, CodeWords, Multiplicities, TempPos, dims, H, i, j, W, V, dimsW, rows, cols, DistBound, VecCount, per, W1, W2, TempVec, TempWeight, S1, Theta, FR, FFail, FTheta,
	R, FRPGE, delta, l1, l2;
    # G, num, mindist, F, debug - are input data
    # CodeWords - if debug[1] = 1, then we record the first 100 different codewords with the lowest weight found so far
    # Multiplicities - number of times codewords from CodeWords were found
    # TempPos - temporal variable corresponding to the position of TempVec in CodeWords
    # dims - dimensions of matrix G
    # H - check matrix
    # i and j - for "for" loop
    # W and V - are vector spaces ortogonal to H and G correspondingly
    # dimsW - dimensions of W (W presented as a matrix, with row being basis vectors)
    # rows and cols are parts of dimW
    # DistBound - upper bound, we are looking for
    # VecCount - number of times vectors with a lowest weight were found so far.
    # per - is for permutations, which we are using in the code
    # W1, W2, TempVec, TempWeight - temporal variables in the code of two matrices, vector and weight
    # S1 - temporal variable for calculation of erasure threshold (it equals to the sum of i_k - index of the last column in information set)
    # Theta - erasure threshold
    # FR, FRPGE, FTheta - complexity exponents
    # FFail - variable, for which 2^(-cols*FFail) gives probability of NOT finding a codeword with lowest weight in a single shot
    # R - rate (lim (k/n) for n-> inf)
    # delta - relative distance (lim (d/n) for n-> inf)
    # l1, l2 - auxiliary variables for taking permutation on pairs

    # Initial parameters
    G := arg[1];;
    num := arg[2];;
    mindist := arg[3];;
    F := arg[4];;

    # Debug parameter
    debug := [0,0,0,0];;
    if (Length(arg) = 5) then
	debug := debug + CoefficientsQadic(arg[5], 2);;
    fi;

    # Dimensions of generator matrix
    dims := DimensionsMat(G);;

    # Creation of check-matrix
    H := CheckMatrix(G, F, debug);;

    # Below we are getting vector spaces W and V ortogonal to the columns of H and G correspondingly.
    # TransposedMatMutable(H) - creates mutable transpose matrix of H.
    # NullspaceMat(R) - returns a list of row vectors that form a basis of the vector space of solutions to the equation vec*R=0.
    W := NullspaceMat(TransposedMatMutable(H));;
    V := NullspaceMat(TransposedMatMutable(G));;

    # There we found dimentions of vector space W (how many vectors in a basis and their length).
    dimsW := DimensionsMat(W);;
    rows := dimsW[1];;
    cols := dimsW[2];;

    # There we define the main variable - bound on the distance of a code. At first, we don't know anything about code,
    # so as a bound we take the value that is bigger than the maximum distance this code could have (((the number of columns in H divided by 2)+1); because number of qudits is col/2).
    # VecCount - number of times vectors with minimal distance were found so far. Initially the value is 0.
    # FirstVecCount - number of times the first vector corresponding to the found upper bound on the distance was found. Initially the value is 0.
    DistBound := Int(cols / 2) + 1;;
    CodeWords := [];;

    # If debug[4] = 1 then we introduce auxiliary variable S1 for the calculation of the erasure threshold
    if (debug[4] = 1) then
	S1 := 0;
    fi;

    # The main part of algorithm.
    for i in [1..num] do
	# We start by creating random permutation for columns in W.
	l1 := ListPerm(Random(SymmetricGroup(cols/2)), cols/2);;  # this is a random permutation of length cols/2, written as a list
	l2 := [];
	# We extend the permutation, so it works now on pairs
	for i in [1..cols/2] do
	    Append(l2, [2*l1[i]-1, 2*l1[i]]);  # l2 contains the permutation we want as a list
	od;
	per := PermList(l2); # this is a permutation of length 2n moving pairs

        # per := Random(SymmetricGroup(cols));; # - old code

	# Perform that permutation.
        W1 := PermutedCols(W, per);;
	# Reduced row echelon form
        W2 := TriangulizedMat(W1);;
	# Inverse permutation
        W2 := PermutedCols(W2,Inverse(per));;

        for j in [1..rows] do
	    # We take one of the sample vectors for this iteration. It supposed to be low-weight.
            TempVec := W2[j];;
	    # Symplectic weight of this sample vector.
            TempWeight := SymplVecWeight(TempVec, F);;

	    # We check if this sample vector satisfies the conditions (this will mean that it corresponds to the logical operator).
	    # First, rough check:
            if (TempWeight > 0) and (TempWeight <= DistBound) then
                if (WeightVecFFE(V * TempVec) > 0) then # linear independence from rows of G
                    if (TempWeight < DistBound) then
                        DistBound := TempWeight;; # min weight found
                        VecCount := 1;; # reset the overall count of vectors of such weight
			
			# Recording all discovered codewords of minimum weight and their multiplicities
			if ((debug[1] = 1) or (debug[4] = 1)) then
			    CodeWords := [TempVec];;
			    Multiplicities := [1];;
			fi;

		    # If we already received such a weight (up to now - it is minimal), we want to update number of vectors, corresponding to it
                    elif (TempWeight = DistBound) then
                        VecCount := VecCount + 1;;

			# Recording of the first 100 different discovered codewords of minimum weight with their multiplicities
                        if ((debug[1] = 1) or (debug[4] = 1)) then
                            TempPos := Position(CodeWords, TempVec);
			    if ((TempPos = fail) and (Length(Multiplicities) < 100)) then
				Add(CodeWords, TempVec);
				Add(Multiplicities, 1);
			    elif (TempPos <> fail) then
				Multiplicities[TempPos] := Multiplicities[TempPos] + 1;;
			    fi;
			fi;

			# Debugging (3-rd bit of debug): program prints some details if the the number of passed iteration rounds is 200*i for integer i.
                        if ((debug[3] = 1) and (RemInt(i, 200) = 0) and (j = 1)) then
                            Print("Distance = ", DistBound, ". The number of times vectors with this weight were found so far: ", VecCount, ". The round = ", i, " of ", num, " iterations.", "\n");
                        fi;

                    fi;

		    # Specific terminator, if we don't want to obtain results for distance below a particular value. 
                    if (DistBound <= mindist) then # not interesting, exit immediately!
			Print("\n", "The found distance equals to or less then specified minimum distance!", "\n");
                        return -DistBound;
                    fi;

                fi;
            fi;
        od;

	if (debug[4] = 1) then
	    S1 := S1 + PositionNonZero(W2[RankMat(G)], 0);;
	fi;

    od;

    # We check that H*c = 0 for c being the first found vector with the lowest weight; that is c is a correct codeword
    if (WeightVecFFE(H * CodeWords[1]) > 0) then 
        Print("\n \n", "The first codeword is not orthogonal to the rows of check matrix H!");

	if (Length(CodeWords[1]) <= 100) then
	    Print("\n", "The improper vector is", "\n");
	    Display(CodeWords[1]);
	fi;

	Error("\n \n");
    fi;

    # If 1-th bit in debug = 1 program prints additional information about the found vector with the lowest weight
    if (debug[1] = 1) then
	if (Length(CodeWords[1]) <= 100) then
	    Print("First vector with the lowest weight found (pay attention to the permuted columns format of matrices and vectors!):", "\n");
	    Display(CodeWords[1]);
	fi;

	Print("The first vector with the lowest symplectic weight ", DistBound, " was found ", Multiplicities[1], " time(s)", "\n");
	Print(i," rounds of ", num," were made.", "\n");
	Print("The overall number of times vectors with the lowest weight were found: ", VecCount, "\n");
	Print("The probability of existence of a vector with lower weight after ", num, " steps: exp(-", Float(AverageCalc(Multiplicities)), ")", "\n");
    fi;

    if (debug[4] = 1) then
	Theta := Float(1 - S1 / (num * cols));;
	Print("Empirical value for erasure threshold: ", Theta, "\n");
	R := RankMat(G) / cols;;
	delta := DistBound / cols;;
	FR := EntrFunc(Float(delta)) - R * EntrFunc(Float(delta / R));;
	# P_Fail = 2^(-cols*F_Fail)^num = 2^(-cols*num*F_Fail); P_Fail = exp(-<n>); <n> = AverageCalc(Multiplicities)
	FFail := AverageCalc(Multiplicities) / (num * cols * Log(2.0));;
	# 2^(-cols*F_RPGE) + 2^(-cols*F_Fail) = 1
	FRPGE := -Log2(1-2^(-cols*FFail))/cols;;
	FTheta := EntrFunc(Float(delta)) - (1 - Theta) * EntrFunc(Float(delta / (1 - Theta)));;
	Print("Complexity exponents: F_R = ", FR, ", F_RPGE = ", FRPGE, ", F_(1-Theta) = ", FTheta, "\n");
	Print("We expect that F_R <= F_RPGE <= F_(1-Theta)", "\n");
	Print("F_Fail = ", FFail, "\n");
	Print("Residual difference of the squares: ", Float(Sum(Multiplicities, x->x^2) - Sum(Multiplicities)^2 / Length(Multiplicities)), "\n");
    fi;

    # Output
    Print("[[", dims[1], ",",dims[1] - RankMat(G), ",", DistBound, "]];", "  Field: ", F, "\n");
    return Concatenation("[[", String(dims[1]), ",", String(dims[1] - RankMat(G)), ",", String(DistBound), "]]");
end;


# Function: adds parameters received from DoRandDist as a third comment line of MatrixMarket file
# Input: StrPath - pathe to the file; params - parameters, received from DoRandDist
# Output: no output.
SaveAddMTX := function (StrPath, params)

    local input, data, i;
    # local variables:
    # input - text, read from file and converted to string
    # data - string of numbers converted to list of them
    # i - for "for" loop

    # Read from file
    input := ReadAll(InputTextFile(StrPath));;
    # Splitting lines
    data := SplitString(input, "\n");;

    # Add parameters as a third line
    Add(data, Concatenation("% Code parameters: ", params), 3);

    # We return "\n" in each line
    for i in [1..(Length(data)-1)] do
	Append(data[i], "\n");
    od;

    # We print all lines over the old ones in the file (with parameters line included)
    PrintTo(StrPath, Concatenation(data));

    Print("The line with code parameters was successfully added to the file!");

end;


# Function: exports a matrix in MatrixMarket format with parameters of the code in the third line of comment section.
# Comments:
# Input: name - name for the file; generator matrix G in the form (A1,B1,A2,B2,...)
# pair - parameter that represents if we should store matrix in a compressed form or not ("pair=integer" - normal storage; "pair=complex" - compressed storage);
# params - string with code parameters;
# F - field, we are using
# comment - list of strings for comment section (this variable can be omitted)
# Output: no output.
SaveNewMTX := function (arg)
    local name, G, params, pair, F, comment, filename, dims, rows, cols, nonzero, i, row, pos;
    # Function was written with ability to support changeable number of variables. We may or may not pass additional comments as the last variable. First six variables are passed in the function.
    # filename - name of the file
    # dims - dimensions of matrix G
    # rows and cols - number of rows and columns of the matrix G
    # nonzero - number of lines needed to store all nonzero elements of matrix with respect to pair parameter
    # i - for "for" loop
    # row and pos - temporal variables; first one will store rows of G, second is for positions of nonzero elements in this row

    # Initial variables
    name := arg[1];;
    G := arg[2];;
    params := arg[3];;
    pair := arg[4];;
    F := arg[5];;

    # We check pair parameter
    if (pair <> "integer") and (pair <> "complex") then
	Error("\n \n", "The \"pair\" parameter is not integer or complex!", "\n \n");
    fi;

    # Full file name
    filename := Concatenation(name, ".mtx");;

    # First two lines in Matrix Market
    PrintTo(filename, "%% MatrixMarket matrix coordinate ", pair, " general", "\n");
    AppendTo(filename,"% Field: ", F, "\n");

    # We write parameters of the code from DoRandDist
    AppendTo(filename, "% Parameters of the code: ", params, "\n");

    # We write comment section in the file, if it was passed into the function
    if (Length(arg) = 6) then
	comment := arg[6];;
	for i in [1..Length(comment)] do
	    AppendTo(filename, "% ", comment[i], "\n");
	od;
    fi;

    # Matrix dimensions
    dims := DimensionsMat(G);;
    rows := dims[1];;
    cols := dims[2];;

    # We count number of line needed for nonzero elements (it depends on the pair; for "pair=integer" it is simply number of nonzero elements in G;
    # for "pair=complex" it is number of pairs G[...,2i] and G[...,2i+1], where at least one of the elements is not 0; we use symplectic weight to implement the second case)
    nonzero := 0;
    if (pair = "integer") then
	for i in [1..rows] do
	    nonzero := nonzero + WeightVecFFE(G[i]);;
	od;
    else
	for i in [1..rows] do
            nonzero := nonzero + SymplVecWeight(G[i], F);;
        od;
    fi;

    # We write dimensions of the matrix and number of line containing nonzero elements
    AppendTo(filename, rows, " ", cols, " ", nonzero, "\n");

    # We append all nonzero elements and their positions into file with respect to pair parameter and field F.
    if (pair = "integer") then
        if IsPrime(Size(F)) then
	    for i in [1..rows] do
		row := G[i];;
		
		pos := PositionNonZero(row, 0);;
		while pos <= cols do
                    AppendTo(filename, i, " ", pos, " ", IntFFE(row[pos]), "\n");
                    pos := PositionNonZero(row, pos);;
		od;
	    od;
	else
	    for i in [1..rows] do
                row := G[i];;

                pos := PositionNonZero(row, 0);;
                while pos <= cols do
                    AppendTo(filename, i, " ", pos, " ", LogFFE(row[pos], PrimitiveElement(F)), "\n");
                    pos := PositionNonZero(row, pos);;
                od;
            od;
	fi;
    else
	if IsPrime(Size(F)) then
	    for i in [1..rows] do
                row := G[i];;

		pos := PositionNonZero(row, 0);;
		while pos <= cols do
		    # For Ai = 0
		    if IsInt(pos/2) then
			AppendTo(filename, i, " ", pos/2, " ", 0, " ", IntFFE(row[pos]), "\n");
			pos := PositionNonZero(row, pos);;
		    # For Ai != 0
		    else
			# Check if Bi = 0
			if (row[pos + 1] = Zero(F)) then
			    AppendTo(filename, i, " ", (pos+1)/2, " ", IntFFE(row[pos]), " ", 0, "\n");
			else
			    AppendTo(filename, i, " ", (pos+1)/2, " ", IntFFE(row[pos]), " ", IntFFE(row[pos + 1]), "\n");
			fi;

			pos := PositionNonZero(row, pos + 1);;
		    fi;
                od;
            od;
        else
	    for i in [1..rows] do
                row := G[i];;

                pos := PositionNonZero(row, 0);;
                while pos <= cols do
                    # For Ai = 0
                    if IsInt(pos/2) then
                        AppendTo(filename, i, " ", pos/2, " ", -1, " ", LogFFE(row[pos], PrimitiveElement(F)), "\n");
                        pos := PositionNonZero(row, pos);;
                    # For Ai != 0
                    else
                        # Check if Bi = 0
                        if (row[pos + 1] = Zero(F)) then
                            AppendTo(filename, i, " ", (pos+1)/2, " ", LogFFE(row[pos], PrimitiveElement(F)), " ", -1, "\n");
                        else
                            AppendTo(filename, i, " ", (pos+1)/2, " ", LogFFE(row[pos], PrimitiveElement(F)), " ", LogFFE(row[pos + 1], PrimitiveElement(F)), "\n");
                        fi;

                        pos := PositionNonZero(row, pos + 1);;
		    fi;
		od;
            od;
        fi;
    fi;

Print("The file was successfully created!");

end;
