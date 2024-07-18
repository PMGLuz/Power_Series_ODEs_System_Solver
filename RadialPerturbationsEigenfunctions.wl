(* ::Package:: *)

BeginPackage["RadialPerturbationsEigenfunctions`"];


EDOSystemReducer::usage = "EDOSystemReducer \
Returns the R and S matrices associated to the system of first-order differential equations.";


getPMatrixNoRec::usage =  "getPMatrixNoRec \
Returns the P matrix associated to the system of first-order differential equations.";


Begin["`Private`"];


(**
Method EDOSystemReducer has two implementations.
Since it uses recurrence, the first call TUMatrix is just the identity matrix.
The next calls TUMatrix must be specified. To avoid to have the user to provide the identity matrix, we have two implementations.
**)
EDOSystemReducer[RMatrix_, \[CapitalTheta]Matrix_, parametrizationVariable_Symbol, SIMPLIFYTIMECONSTRAINT_:300] := Block[{U, InvU, Uprime, T, InvT, ExpansionOrder,
	tempSMatrix, transformed\[CapitalTheta]Matrix, tempTUMatrix,
	temp, r=parametrizationVariable, tempJordanMatrix, Value, TUMatrix = IdentityMatrix[Length[RMatrix]]},
	
	(*Construct matrix T containing the generalized eigenvectors of RMatrix.
	Always use since this guarantees that we get the right order of the eigenvalues: T^-1RT has lowest eigenvalue in 1,1-position*)
	{T,tempJordanMatrix}=JordanDecomposition[RMatrix];
	
	(*Get maximum difference between eigenvalues*)
	ExpansionOrder=Abs[Min[Diagonal[tempJordanMatrix]]-Max[Diagonal[tempJordanMatrix]]];
	
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
	transformed\[CapitalTheta]Matrix=Simplify[InvU . InvT . \[CapitalTheta]Matrix . T . U-temp, TimeConstraint->SIMPLIFYTIMECONSTRAINT];
	
	{tempSMatrix,transformed\[CapitalTheta]Matrix,tempTUMatrix}=EDOSystemReducer[tempSMatrix, transformed\[CapitalTheta]Matrix, TUMatrix . T . U, r]
];


EDOSystemReducer[RMatrix_, \[CapitalTheta]Matrix_, TUMatrix_, parametrizationVariable_Symbol, SIMPLIFYTIMECONSTRAINT_:300] := Block[{U,InvU,Uprime,T,InvT,ExpansionOrder,
	tempSMatrix,transformed\[CapitalTheta]Matrix,tempTUMatrix,
	temp, r=parametrizationVariable, tempJordanMatrix, Value},
	
	(*Construct matrix T containing the generalized eigenvectors of RMatrix.
	Always use since this guarantees that we get the right order of the eigenvalues: T^-1RT has lowest eigenvalue in 1,1-position*)
	{T,tempJordanMatrix}=JordanDecomposition[RMatrix];
	
	(*Get maximum difference between eigenvalues*)
	ExpansionOrder=Abs[Min[Diagonal[tempJordanMatrix]]-Max[Diagonal[tempJordanMatrix]]];
	
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
	
	{tempSMatrix,transformed\[CapitalTheta]Matrix,tempTUMatrix}=EDOSystemReducer[tempSMatrix, transformed\[CapitalTheta]Matrix, TUMatrix . T . U, r]
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
