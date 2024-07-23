(* ::Package:: *)

BeginPackage["PowerSeriesODESystemSolver`"];


ODESystemReducer::usage = "ODESystemReducer \
Returns the R and S matrices associated to the system of first-order differential equations.";


ODESystemReducerGeneral::usage = "ODESystemReducer \
Returns the R and S matrices associated to the system of first-order differential equations.
This method is valid any diagonalizable real R matrix. It is, therefore slower than EDOSystemReducerInteger.";


ODESystemReducerInteger::usage = "ODESystemReducer \
Returns the R and S matrices associated to the system of first-order differential equations.
This method specializes in the case when the R matrix only has integer eigenvalues.";


getPMatrixNoRec::usage =  "getPMatrixNoRec \
Returns the P matrix associated to the system of first-order differential equations.";


Begin["`Private`"];


(***
OrderJordanDecomposition\
Returns the JordanDecomposition system with the order where the non-integer eigenvalues are in the first columns and the last columns have the integer eigenvalues.\
The integer eigenvalues are ordered such that the highest eigenvalue is always in the last position of the diagonal.
***)

SortedJordanDecomposition[InitMatrix_]:=Block[{T,Q,JMatrix,diagJMatrix,swappedJMatrix,eigenVList,tempMatrix,currentPositionSwappedJMatrix,currentPositionJMatrix},
	{T,JMatrix}=JordanDecomposition[InitMatrix];
	
	(*Get distinct eigenvalues and dimension of the associated subspace*)
	diagJMatrix=Diagonal[JMatrix];
	eigenVList=Tally[diagJMatrix];
	
	swappedJMatrix=ConstantArray[0,Dimensions[JMatrix]];
	currentPositionSwappedJMatrix=1;
	(*This guarantees that the non-integer eigenvalues are in the first columns*)
	For[i=1,i<=Length[eigenVList],i++,
		If[!IntegerQ[eigenVList[[i,1]]],
			currentPositionJMatrix=FirstPosition[diagJMatrix,eigenVList[[i,1]]][[1]];
			swappedJMatrix[[currentPositionSwappedJMatrix;;currentPositionSwappedJMatrix+eigenVList[[i,2]]-1,
				currentPositionSwappedJMatrix;;currentPositionSwappedJMatrix+eigenVList[[i,2]]-1]] = 
				JMatrix[[currentPositionJMatrix;;currentPositionJMatrix+eigenVList[[i,2]]-1,
					currentPositionJMatrix;;currentPositionJMatrix+eigenVList[[i,2]]-1]];
			
			currentPositionSwappedJMatrix=currentPositionSwappedJMatrix+eigenVList[[i,2]];
		,
			Continue
		]
	];
	
	(*The next line guarantees that we have the correct order for the integer eigenvalues*)
	eigenVList=Sort[eigenVList,#1[[1]]<#2[[1]]&];
	(*Then we write the integer eigenvalues in the next columns. We write first the lower value eigenvalues*)
	For[i=1,i<=Length[eigenVList],i++,
		If[IntegerQ[eigenVList[[i,1]]],
			currentPositionJMatrix=FirstPosition[diagJMatrix,eigenVList[[i,1]]][[1]];
			swappedJMatrix[[currentPositionSwappedJMatrix;;currentPositionSwappedJMatrix+eigenVList[[i,2]]-1,
				currentPositionSwappedJMatrix;;currentPositionSwappedJMatrix+eigenVList[[i,2]]-1]] = 
				JMatrix[[currentPositionJMatrix;;currentPositionJMatrix+eigenVList[[i,2]]-1,
					currentPositionJMatrix;;currentPositionJMatrix+eigenVList[[i,2]]-1]];
			
			currentPositionSwappedJMatrix=currentPositionSwappedJMatrix+eigenVList[[i,2]];
		,
			Continue
		]
	];
	
	
	(****
	Now we need to get change of basis matrix eigenvectors swappedJMatrix= Q JMatrix Q^-1. 
	Since mathematica is consistent with the ordering of the JordanDecompositions we can use that.
	*****)
	
	Q=JordanDecomposition[swappedJMatrix][[1]];
	Q=T . Inverse[Q];
	
	Return[{Q,swappedJMatrix}]
]


ODESystemReducer[RMatrix_, \[CapitalTheta]Matrix_, parametrizationVariable_Symbol, SIMPLIFYTIMECONSTRAINT_:300] := Block[{ExpansionOrder,eigenVList,
	 r=parametrizationVariable},
	
	eigenVList=Diagonal[JordanDecomposition[RMatrix][[2]]];
	(*Check if RMatrix only has non integer eigenvalues. If so, return RMatrix and \[CapitalTheta]Matrix*)
	If[NoneTrue[eigenVList,IntegerQ],
		Return[ {RMatrix,\[CapitalTheta]Matrix,IdentityMatrix[Length[RMatrix]]} ]
		,
		Continue
	];
	
	(*Check if RMatrix only has integer eigenvalues. If so, call ODESystemReduceInteger and return the result*)
	If[AllTrue[eigenVList,IntegerQ],
		Return[ ODESystemReducerInteger[RMatrix,\[CapitalTheta]Matrix,r] ]
		,
		Continue
	];
	
	(*If RMatrix has both integer and non integer eigenvalues, call ODESystemReduceGeneral and return the result*)
	Return[ ODESystemReducerGeneral[RMatrix,\[CapitalTheta]Matrix,r] ]
];


(**
Method EDOSystemReducer has two implementations.
Since it uses recurrence, the first call TUMatrix is just the identity matrix.
The next calls TUMatrix must be specified. To avoid to have the user to provide the identity matrix, we have two implementations.
**)

ODESystemReducerGeneral[RMatrix_, \[CapitalTheta]Matrix_, parametrizationVariable_Symbol, SIMPLIFYTIMECONSTRAINT_:300] := Block[{U, InvU, Uprime, T, InvT, ExpansionOrder,
	tempSMatrix, transformed\[CapitalTheta]Matrix, tempTUMatrix,
	temp, r=parametrizationVariable, tempJordanMatrix, Value, TUMatrix = IdentityMatrix[Length[RMatrix]]},
	
	(*Construct matrix T containing the generalized eigenvectors of RMatrix.
	Always use since this guarantees that we get the right order of the eigenvalues: T^-1RT has lowest eigenvalue in last position of the diagonal*)
	{T,tempJordanMatrix}=SortedJordanDecomposition[RMatrix];
	
	(*Get the integer eigenvalues and set ExpansionOrder as maximum difference between those eigenvalues*)
	temp= Cases[Diagonal[tempJordanMatrix] ,_Integer];
	ExpansionOrder=Abs[Min[temp]-Max[temp]];
	
	If[ExpansionOrder==0,
		Return[{RMatrix,\[CapitalTheta]Matrix,TUMatrix}],
		Continue
	];
	
	(*Avoid denominators that mathematica adds. Do this only if not the last one, so after the previous line!*)
	Block[{LineN,ColumnN,tempDenominator,tempT},
		For[ColumnN=1,ColumnN<=Dimensions[T][[2]],ColumnN++,
			For[LineN=1,LineN<=Dimensions[T][[1]],LineN++,
				tempDenominator=Denominator[T[[LineN,ColumnN]]//Together];
				tempT=T//Transpose;
				tempT[[ColumnN]]=tempT[[ColumnN]]*tempDenominator//Simplify;
				T=tempT//Transpose;
			]
		]
	];
	
	(*Construct inverse matrix T. Do after the return above to save time*)
	InvT=Inverse[T];
	
	(*Construct U matrix*)
	U=IdentityMatrix[Length[RMatrix]];
	U[[Length[RMatrix],Length[RMatrix]]]=r;
	InvU=Inverse[U];
	
	
	(*Construct U' matrix*)
	Uprime=0*IdentityMatrix[Length[RMatrix]];
	Uprime[[Length[RMatrix],Length[RMatrix]]]=1;
	
	
	transformed\[CapitalTheta]Matrix=Series[\[CapitalTheta]Matrix,{r,0,ExpansionOrder}]//Normal;
	transformed\[CapitalTheta]Matrix=InvU . InvT . transformed\[CapitalTheta]Matrix . T . U;
	
	temp=1/r*Series[r*transformed\[CapitalTheta]Matrix,{r,0,0}]//Normal;
	
	tempSMatrix=r*Simplify[
		1/r*InvT . RMatrix . T-InvU . Uprime+temp,
		TimeConstraint->SIMPLIFYTIMECONSTRAINT
	];
	
	transformed\[CapitalTheta]Matrix=Simplify[InvU . InvT . \[CapitalTheta]Matrix . T . U-temp, TimeConstraint->SIMPLIFYTIMECONSTRAINT];
	
	{tempSMatrix,transformed\[CapitalTheta]Matrix,tempTUMatrix}=ODESystemReducerGeneral[tempSMatrix, transformed\[CapitalTheta]Matrix, TUMatrix . T . U, r]
];


ODESystemReducerGeneral[RMatrix_, \[CapitalTheta]Matrix_, TUMatrix_, parametrizationVariable_Symbol, SIMPLIFYTIMECONSTRAINT_:300] := Block[{U,InvU,Uprime,T,InvT,ExpansionOrder,
	tempSMatrix,transformed\[CapitalTheta]Matrix,tempTUMatrix,
	temp, r=parametrizationVariable, tempJordanMatrix, Value},
	
	(*Construct matrix T containing the generalized eigenvectors of RMatrix.
	Always use since this guarantees that we get the right order of the eigenvalues: T^-1RT has lowest eigenvalue in last position of the diagonal*)
	{T,tempJordanMatrix}=SortedJordanDecomposition[RMatrix];
	
	(*Get the integer eigenvalues and set ExpansionOrder as maximum difference between those eigenvalues*)
	temp= Cases[Diagonal[tempJordanMatrix] ,_Integer];
	ExpansionOrder=Abs[Min[temp]-Max[temp]];
	
	If[ExpansionOrder==0,
		Return[{RMatrix,\[CapitalTheta]Matrix,TUMatrix}],
		Continue
	];
	
	(*Avoid denominators that mathematica adds. Do this only if not the last one, so after the previous line!*)
	Block[{LineN,ColumnN,tempDenominator,tempT},
		For[ColumnN=1,ColumnN<=Dimensions[T][[2]],ColumnN++,
			For[LineN=1,LineN<=Dimensions[T][[1]],LineN++,
				tempDenominator=Denominator[T[[LineN,ColumnN]]//Together];
				tempT=T//Transpose;
				tempT[[ColumnN]]=tempT[[ColumnN]]*tempDenominator//Simplify;
				T=tempT//Transpose;
			]
		]
	];
	
	(*Construct inverse matrix T. Do after the return above to save time*)
	InvT=Inverse[T];
	
	(*Construct U matrix*)
	U=IdentityMatrix[Length[RMatrix]];
	U[[Length[RMatrix],Length[RMatrix]]]=r;
	InvU=Inverse[U];
	
	
	(*Construct U' matrix*)
	Uprime=0*IdentityMatrix[Length[RMatrix]];
	Uprime[[Length[RMatrix],Length[RMatrix]]]=1;
	
	
	transformed\[CapitalTheta]Matrix=Series[\[CapitalTheta]Matrix,{r,0,ExpansionOrder}]//Normal;
	transformed\[CapitalTheta]Matrix=InvU . InvT . transformed\[CapitalTheta]Matrix . T . U;
	
	temp=1/r*Series[r*transformed\[CapitalTheta]Matrix,{r,0,0}]//Normal;
	
	tempSMatrix=r*Simplify[
		1/r*InvT . RMatrix . T-InvU . Uprime+temp,
		TimeConstraint->SIMPLIFYTIMECONSTRAINT
	];
	
	transformed\[CapitalTheta]Matrix=Simplify[InvU . InvT . \[CapitalTheta]Matrix . T . U-temp,TimeConstraint->SIMPLIFYTIMECONSTRAINT];
	
	{tempSMatrix,transformed\[CapitalTheta]Matrix,tempTUMatrix}=ODESystemReducerGeneral[tempSMatrix, transformed\[CapitalTheta]Matrix, TUMatrix . T . U, r]
];


(**
Method EDOSystemReducer has two implementations.
Since it uses recurrence, the first call TUMatrix is just the identity matrix.
The next calls TUMatrix must be specified. To avoid to have the user to provide the identity matrix, we have two implementations.
**)



ODESystemReducerInteger[RMatrix_, \[CapitalTheta]Matrix_, parametrizationVariable_Symbol, SIMPLIFYTIMECONSTRAINT_:300] := Block[{U, InvU, Uprime, T, InvT, ExpansionOrder,
	tempSMatrix, transformed\[CapitalTheta]Matrix, tempTUMatrix,
	temp, r=parametrizationVariable, tempJordanMatrix, Value, TUMatrix = IdentityMatrix[Length[RMatrix]]},
	
	(*Construct matrix T containing the generalized eigenvectors of RMatrix.
	Always use since this guarantees that we get the right order of the eigenvalues: T^-1RT has lowest eigenvalue in the last position of the diagonal*)
	{T,tempJordanMatrix}=JordanDecomposition[RMatrix];
	
	(*Get maximum difference between eigenvalues*)
	temp=Diagonal[tempJordanMatrix];
	ExpansionOrder=Abs[Min[temp]-Max[temp]];
	
	If[ExpansionOrder==0,
		Return[{RMatrix,\[CapitalTheta]Matrix,TUMatrix}],
		Continue
	];
	
	
	(*Avoid denominators that mathematica adds. Do this only if not the last one, so after the previous line!*)
	Block[{LineN,ColumnN,tempDenominator,tempT},
		For[ColumnN=1,ColumnN<=Dimensions[T][[2]],ColumnN++,
			For[LineN=1,LineN<=Dimensions[T][[1]],LineN++,
				tempDenominator=Denominator[T[[LineN,ColumnN]]//Together];
				tempT=T//Transpose;
				tempT[[ColumnN]]=tempT[[ColumnN]]*tempDenominator//Simplify;
				T=tempT//Transpose;
			]
		]
	];
	
	(*Construct inverse matrix T. Do after the return above to save time*)
	InvT=Inverse[T];
	
	(*Construct U matrix*)
	U=IdentityMatrix[Length[RMatrix]];
	U[[Length[RMatrix],Length[RMatrix]]]=r;
	InvU=Inverse[U];
	
	
	(*Construct U' matrix*)
	Uprime=0*IdentityMatrix[Length[RMatrix]];
	Uprime[[Length[RMatrix],Length[RMatrix]]]=1;
	
	
	transformed\[CapitalTheta]Matrix=Series[\[CapitalTheta]Matrix,{r,0,ExpansionOrder}]//Normal;
	transformed\[CapitalTheta]Matrix=InvU . InvT . transformed\[CapitalTheta]Matrix . T . U;
	
	temp=1/r*Series[r*transformed\[CapitalTheta]Matrix,{r,0,0}]//Normal;
	
	tempSMatrix=r*Simplify[
		1/r*InvT . RMatrix . T-InvU . Uprime+temp,
		TimeConstraint->SIMPLIFYTIMECONSTRAINT
	];
	
	transformed\[CapitalTheta]Matrix=Simplify[InvU . InvT . \[CapitalTheta]Matrix . T . U-temp, TimeConstraint->SIMPLIFYTIMECONSTRAINT];
	
	{tempSMatrix,transformed\[CapitalTheta]Matrix,tempTUMatrix}=ODESystemReducerInteger[tempSMatrix, transformed\[CapitalTheta]Matrix, TUMatrix . T . U, r]
];


ODESystemReducerInteger[RMatrix_, \[CapitalTheta]Matrix_, TUMatrix_, parametrizationVariable_Symbol, SIMPLIFYTIMECONSTRAINT_:300] := Block[{U,InvU,Uprime,T,InvT,ExpansionOrder,
	tempSMatrix,transformed\[CapitalTheta]Matrix,tempTUMatrix,
	temp, r=parametrizationVariable, tempJordanMatrix, Value},
	
	(*Construct matrix T containing the generalized eigenvectors of RMatrix.
	Always use since this guarantees that we get the right order of the eigenvalues: T^-1RT has lowest eigenvalue in 1,1-position*)
	{T,tempJordanMatrix}=JordanDecomposition[RMatrix];
	
	(*Get maximum difference between eigenvalues*)
	temp=Diagonal[tempJordanMatrix];
	ExpansionOrder=Abs[Min[temp]-Max[temp]];
	
	If[ExpansionOrder==0,
		Return[{RMatrix,\[CapitalTheta]Matrix,TUMatrix}],
		Continue
	];
	
	(*Avoid denominators that mathematica adds. Do this only if not the last one, so after the previous line!*)
	Block[{LineN,ColumnN,tempDenominator,tempT},
		For[ColumnN=1,ColumnN<=Dimensions[T][[2]],ColumnN++,
			For[LineN=1,LineN<=Dimensions[T][[1]],LineN++,
				tempDenominator=Denominator[T[[LineN,ColumnN]]//Together];
				tempT=T//Transpose;
				tempT[[ColumnN]]=tempT[[ColumnN]]*tempDenominator//Simplify;
				T=tempT//Transpose;
			]
		]
	];
	
	(*Construct inverse matrix T. Do after the return above to save time*)
	InvT=Inverse[T];
	
	(*Construct U matrix*)
	U=IdentityMatrix[Length[RMatrix]];
	U[[Length[RMatrix],Length[RMatrix]]]=r;
	InvU=Inverse[U];
	
	
	(*Construct U' matrix*)
	Uprime=0*IdentityMatrix[Length[RMatrix]];
	Uprime[[Length[RMatrix],Length[RMatrix]]]=1;
	
	
	transformed\[CapitalTheta]Matrix=Series[\[CapitalTheta]Matrix,{r,0,ExpansionOrder}]//Normal;
	transformed\[CapitalTheta]Matrix=InvU . InvT . transformed\[CapitalTheta]Matrix . T . U;
	
	temp=1/r*Series[r*transformed\[CapitalTheta]Matrix,{r,0,0}]//Normal;
	
	tempSMatrix=r*Simplify[
		1/r*InvT . RMatrix . T-InvU . Uprime+temp,
		TimeConstraint->SIMPLIFYTIMECONSTRAINT
	];
	
	(*transformed\[CapitalTheta]Matrix=transformed\[CapitalTheta]Matrix-temp//Simplify;*)
	transformed\[CapitalTheta]Matrix=Simplify[InvU . InvT . \[CapitalTheta]Matrix . T . U-temp,TimeConstraint->SIMPLIFYTIMECONSTRAINT];
	
	{tempSMatrix,transformed\[CapitalTheta]Matrix,tempTUMatrix}=ODESystemReducerInteger[tempSMatrix, transformed\[CapitalTheta]Matrix, TUMatrix . T . U, r]
];


getPMatrixNoRec[RMatrix_, \[CapitalTheta]Matrix_, SeriesMaxOrder_, parametrizationVariable_Symbol] := Block[
{tempPMatrix, \[CapitalTheta]MatrixCoefficients, dummyFunc, CMatrix, tempMatrix,variablesList, Order, r=parametrizationVariable},
	
	(*Do this to avoid to append to list*)
	tempPMatrix= Table[
		Array[dummyFunc,Dimensions[RMatrix]],
		{n,0,SeriesMaxOrder}
	];
	
	(*Pre-compute series coefficients for \[CapitalTheta] matrix*)
	\[CapitalTheta]MatrixCoefficients=Table[SeriesCoefficient[\[CapitalTheta]Matrix, {r,0,n}], {n,0,SeriesMaxOrder-1}];
	
	
	tempPMatrix[[1]]=IdentityMatrix[Length[RMatrix]];
	variablesList = Array[dummyFunc, Dimensions[RMatrix]]//Flatten;
	For[Order=1, Order<=SeriesMaxOrder, Order++,
		CMatrix=Sum[\[CapitalTheta]MatrixCoefficients[[Order-j]] . tempPMatrix[[j+1]], {j,0,Order-1}];
		
		tempMatrix= Array[dummyFunc, Dimensions[RMatrix]];
		
		
		tempMatrix= Array[dummyFunc, Dimensions[RMatrix]]/.Solve[
			tempMatrix . ((Order)*IdentityMatrix[Dimensions[RMatrix]]+RMatrix) - RMatrix . tempMatrix == CMatrix, variablesList
		][[1]]//Simplify;
		
		tempPMatrix[[Order+1]]=tempMatrix
	];
	
	tempPMatrix=tempPMatrix*Table[r^n, {n, 0, SeriesMaxOrder}];
	tempPMatrix//Total
]


End[];


EndPackage[];
