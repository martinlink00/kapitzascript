#!/usr/bin/env wolframscript
(* ::Package:: *)

(*

Notes for usage:

- Effective potential terms will be calculated up to Subscript[V, 2*highestcoeff] and the boolean incond can be switched on if initial conditions are to be calculated.
- Before running the script, the kernel needs to be reset in order to delete all global parameters from memory.
- Note, that results are dumpsaved into the subfolder "exports", in order to save RAM and in order to not lose results if Mathematica crashes which it tends to do.

*)




(**************************************************************************)




total=
Timing[


(*INIT*)

ClearSystemCache[];
Off[FrontEndObject::notavail];

highestcoeffresp = ReadList["params.txt"][[1]];

If[Not[DirectoryQ["exports"]],
  Print["Creating exports directory."];CreateDirectory["exports"];SetDirectory["exports"];CreateDirectory["logs"];,
  SetDirectory["exports"];];



(*MISCDEF*)

logfilename = "pot_"<>DateString["ISODateTime"]<>".log";

LogPrint[x_,type_] := Module[{},
	Print[x];
	outp = "["<>DateString["ISODateTime"]<>"]\t["<>type<>"]\t[Maximum RAM: "<>ToString[Round[N[MaxMemoryUsed[]/(10^9) *1.0],0.01]]<>" GB ("<>ToString[Round[N[MaxMemoryUsed[]/MemoryAvailable[] *100.0],0.01]]<>"%)]\t"<>x;
	WriteLine[FileNameJoin[{"logs",logfilename}],outp];];

myByteCount[symbolName_String]:=
	Replace[
		ToHeldExpression[symbolName],
		Hold[x__]:>If[MemberQ[Attributes[x],
		Protected|ReadProtected],
		Sequence@@{},
		{ByteCount[Through[{OwnValues,DownValues,UpValues,SubValues,DefaultValues,FormatValues,NValues}[Unevaluated@x,Sort->False]]],symbolName}]
	];
CheckRAM[] :=
	With[
		{listing=myByteCount/@Names[]},
		Labeled[Grid[Reverse@Take[Sort[listing],-100],
		Frame->True,Alignment->Left],
		Column[{Style["ByteCount for symbols without attributes Protected and ReadProtected in all contexts",16,FontFamily->"Times"],
		Style[Row@{"Total: ",Total[listing[[All,1]]]," bytes for ",Length[listing]," symbols"},Bold]},Center,1.5],Top]
	];
PD[]:=Print["---------------------------------------------------------------------"];

AppendTo[$Messages,OpenWrite[FileNameJoin[{"logs",logfilename}]]];


(*INIT PRINTOUT*)

PD[];

If[$VersionNumber!=12.1,LogPrint["Version error: Script was developed on and for version 12.1.","ERROR"];Quit[],LogPrint["### KAPITZA POTENTIAL SCRIPT ###","RUNNING"];LogPrint["### ACQUIRING PARAMETERS FROM PARAM FILE ###","RUNNING"]];

PD[];

(*\[Rho]inp = Input["Input value for rho (float): "];\[Omega]inp = Input["Input value for om_0 (float): "];*)
highestcoeff = If[IntegerQ[highestcoeffresp],highestcoeffresp,LogPrint["Wrong input in params.txt file!","ERROR"];Quit[]];
(*\[Rho] = If[NumberQ[\[Rho]inp],\[Rho]inp,LogPrint["Wrong input!","ERROR"];Quit[]];
\[Omega] = If[NumberQ[\[Omega]inp],\[Omega]inp,LogPrint["Wrong input!","ERROR"];Quit[]];*)



PD[];
userinfo = "Script has been executed on "<>ToString[$OperatingSystem]<>" by the user " <>ToString[$UserName]<>".";
versioninfo = "Running on Mathematica "<>ToString[$VersionNumber]<>".";
processorinfo = ToString[$ProcessorCount]<>" processors of type "<>ToString[$ProcessorType]<>" available.";
memoryinfo = ToString[N[MemoryAvailable[]/10^6]]<>" MB of memory available for allocation.";
paraminfo = "Calculating up to non-trivial correction "<>ToString[highestcoeff]<>"."(*<>" with params (rho,om_0) = ("<>ToString[\[Rho]]<>","<>ToString[\[Omega]]<>")."*);
esttimeinfo = "In order to estimate the time and memory that will be needed on your machine from available logfiles, run the script paramest.wls.";
exportinfo = "Results will be dumpsaved into the subfolder exports and can be accessed using Results.nb.";
loginfo = "Log for this script execution is saved under exports/logs/"<>logfilename<>".";
LogPrint["### SYSTEM INFORMATION ###","RUNNING"];
LogPrint[userinfo,"RUN:SYS"];LogPrint[versioninfo,"RUN:SYS"];LogPrint[processorinfo,"RUN:SYS"];LogPrint[memoryinfo,"RUN:SYS"];PD[];
LogPrint[paraminfo,"RUN:PRM"];LogPrint[esttimeinfo,"INFORM"];LogPrint[loginfo,"INFORM"];LogPrint[exportinfo,"INFORM"];PD[];
LogPrint["### NOW RUNNING SCRIPT ###","RUNNING"];PD[];



(*SCRIPT*)

LogPrint["Initiating.","PROGRESS"];
\[Epsilon]max = 2 * highestcoeff + 3; (*Max collected \[Epsilon]*)
highestpowercorr=Round[(\[Epsilon]max-2.1)/2]; (*Highest power correction*)
e[n_,t]:=Exp[I \[Omega] n t];
\[CapitalPhi][n_]= \[Epsilon]^(Abs[n]) *Sum[\[Phi][n,m][t]\[Epsilon]^m,{m,0,\[Epsilon]max}];
\[Delta]\[Theta]= Sum[If[n==0,0,e[n,t]\[CapitalPhi][n]],{n,-\[Epsilon]max,\[Epsilon]max}]; (*\[Delta]\[Theta] Ansatz*)
\[Theta]= \[CapitalPhi][0] +\[Delta]\[Theta];
oddphi0repl =Flatten[Table[D[\[Phi][0,2i+1][t],{t,k}]->0,{i,0,highestcoeff+2},{k,0,\[Epsilon]max}]];
newnot = Flatten[Table[D[\[Phi][0,2j][t],{t,k}]->D[\[CurlyPhi][2j][t],{t,k}],{j,0,highestcoeff+2},{k,0,\[Epsilon]max}],2];
enrepl = Table[Exp[I \[Omega] n t]->x^n,{n,-\[Epsilon]max,\[Epsilon]max}];
PD[];

algtime =
Timing[

\[Theta]sercoeff = Series[\[Theta]/.oddphi0repl/.newnot,{\[Epsilon],0,\[Epsilon]max}][[3]];
\[Omega]sercoeff = D[\[Theta]sercoeff,t]/.enrepl;
\[Theta]sercoeff = \[Theta]sercoeff/.enrepl;
dercoeff = Series[\[Epsilon]^2 D[\[CapitalPhi][0],{t,2}] + Sum[If[n==0,0,e[n,t](\[Epsilon]^2 D[\[CapitalPhi][n],{t,2}] + 2 I \[Omega] n \[Epsilon] D[\[CapitalPhi][n],{t,1}]-\[Omega]^2 n^2 \[CapitalPhi][n])],{n,-\[Epsilon]max,\[Epsilon]max}]/.oddphi0repl /.newnot/.enrepl,{\[Epsilon],0,\[Epsilon]max}][[3]];
sinser = Series[Sin[\[CurlyPhi][0][t]+\[Delta]],{\[Delta],0,\[Epsilon]max-1}][[3]];
bellargument = Table[Table[Factorial[i]\[Theta]sercoeff[[i+1]],{i,1,n}],{n,1,\[Epsilon]max}];
sinexpfunc[epsord_] := If[epsord<0,0,If[epsord==0,Sin[\[CurlyPhi][0][t]],1/Factorial[epsord] Sum[Factorial[n]sinser[[n+1]]BellY[epsord,n,bellargument[[epsord-n+1]]],{n,1,epsord}]]];


DumpSave["ini.mx",{(*\[Rho],\[Omega],*)LogPrint,logfilename,\[Theta],\[Delta]\[Theta],\[CapitalPhi],\[Epsilon]max,highestcoeff,myByteCount,CheckRAM,PD,newnot,oddphi0repl}]; (*Dumping basic definitions for later usage.*)


(*Sin treatment*)

LogPrint["Acquiring sin(theta) expansion.","PRG:INIT"];
sinexplist = Table[sinexpfunc[ee],{ee,0,\[Epsilon]max-1}]; (*Need only up to \[Epsilon]max-1, since one \[Epsilon] already on lhs.*)
sinexp[\[Epsilon]ord_] := If[\[Epsilon]ord<0,0,sinexplist[[\[Epsilon]ord+1]]];
dgl[\[Epsilon]ord_,nord_] := (Coefficient[Expand[dercoeff[[\[Epsilon]ord]]],x,nord]==Coefficient[Expand[\[Omega]^2 (sinexp[\[Epsilon]ord-2]-\[Rho]  * 1/(2 I) * (x-1/x)sinexp[\[Epsilon]ord-1])],x,nord]);
Share[]; (*In hope of saving memory.*)



(*Declare \[Phi][n\[NotEqual]0,m,t] solutions*)

PD[];LogPrint["Acquiring algebraic solutions up to eps^"<>ToString[\[Epsilon]max]<>".","PROGRESS"];


For[i=1,i<\[Epsilon]max+1,i++,
	PD[];LogPrint["Finding relevant terms at order eps^"<>ToString[i]<>".","PRG:ALG"];
	ex = Expand[\[Omega]^2 (sinexp[i-2]-\[Rho]  * 1/(2 I) * (x-x^-1)sinexp[i-1]),Trig->False];
	Monitor[
		For[j=0,j<i,j++,
			If[i==1&&j==0,Clear[\[Phi]],<<phi.mx];
			cef = Coefficient[ex,x,i-j];
			If[j==0,
				solex=-1/((i-j)^2 \[Omega]^2) (cef),
				If[j==1,
				solex=-1/((i-j)^2 \[Omega]^2) (cef - 2 I \[Omega] (i-j)D[\[Phi][i-j,j-1][t],t]),
				solex=-1/((i-j)^2 \[Omega]^2) (cef- 2 I \[Omega] (i-j)D[\[Phi][i-j,j-1][t],t]-D[\[Phi][i-j,j-2][t],{t,2}])
				];
			];
			LogPrint["Declaring algebraic solution of phi[+-"<>ToString[i-j]<>","<>ToString[j]<>"][t].","PRG:ALG"];
			\[Phi][i-j,j][t_] = ComplexExpand[solex];
			\[Phi][-(i-j),j][t_] = ComplexExpand[solex\[Conjugate]];
			DumpSave["phi.mx",\[Phi]];
			Clear[\[Phi]];
		];
	,ProgressIndicator[j,{0,i-1}]];
];


ClearAll[ex,cef,solex];
<<phi.mx;

][[1]];

PD[];LogPrint["All phi[n!=0,m,t] have been solved and declared. This took "<>ToString[algtime]<>" s.","INFORM"];PD[];

(*Effective potential treatment*)

DumpSave["ser.mx",{\[Theta]sercoeff,\[Omega]sercoeff}];


vefftime = Timing[


LogPrint["Acquiring acceleration solutions with inserted algebraic solutions.","PROGRESS"];
acc0[ord_] := Coefficient[Expand[\[Omega]^2 (sinexp[ord-2]-\[Rho]  * 1/(2 I) * (x-1/x)sinexp[ord-1]),Trig->False],x,0];
acc0list = Table[acc0[2 ee+2],{ee,0,highestpowercorr}];

For[m=0,m<highestcoeff+1,m++,
	For[n=0,n<\[Epsilon]max-1,n++,
		\!\(\*SuperscriptBox[\(\[CurlyPhi][2  m]\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", "+", "n"}], ")"}],
Derivative],
MultilineFunction->None]\)[t]=D[a[2m][t],{t,n}];
	];
];

PD[];LogPrint["Defining and expanding effective potential.","PRG:ACC"];

Veff[z_] = Sum[\[Epsilon]^(2 m) V[2m][z],{m,0,(\[Epsilon]max-3)/2}];
veffser = Series[D[Veff[z],z]/.z->\[CapitalPhi][0]/.oddphi0repl/.newnot/.\[CurlyPhi][0][t]->x,{\[Epsilon],0,\[Epsilon]max-3}];

PD[];LogPrint["Acquiring derivative expressions from energy conservation.","PRG:ACC"];

(*E0 = Sum[En[2m]\[Epsilon]^(2 m),{m,0,highestcoeff}];*)
derrepleq = Table[Coefficient[D[\[CapitalPhi][0],t]^2,\[Epsilon],2 \[Epsilon]ord]==Coefficient[2(E0-Series[Veff[\[CapitalPhi][0]],{\[Epsilon],0,\[Epsilon]max}]),\[Epsilon],2 \[Epsilon]ord],{\[Epsilon]ord,0,highestcoeff}]/.oddphi0repl/.newnot; 
dertosolve = Table[\[CurlyPhi][2j]'[t],{j,0,highestcoeff}];
derrepl = Simplify[Solve[derrepleq,dertosolve][[1]]];


PD[];LogPrint["Declaring varphi_2m''[t] and derivatives.","PRG:ACC"];

Monitor[
	For[m=0,m<highestcoeff,m++,
		a[2m][t_]=acc0list[[m+1]];
	],
	ProgressIndicator[m,{0,highestcoeff-1}]
];

DumpSave["a.mx",a];
ClearAll[a];
<<a.mx;

PD[];LogPrint["Declaring varphi_2m'[t].","PRG:ACC"];

Monitor[
	For[m=0,m<highestcoeff+1,m++,
		\[CurlyPhi][2m]'[t]=\[CurlyPhi][2m]'[t]/.derrepl[[m+1]]; 
	],
	ProgressIndicator[m,{0,highestcoeff+1}]
];


PD[];LogPrint["Acquiring effective potential expressions.","PROGRESS"];
LogPrint["Expanding acceleration expressions.","PRG:POT"];
acc0list = acc0list/.\[CurlyPhi][0][t]->x;
DumpSave["acclist.mx",acc0list];
ClearAll[acc0list];
<<acclist.mx;

LogPrint["Acquiring expressions for V_2m'[x].","PRG:POT"];

irrelevantrepl = Table[\[CurlyPhi][2m][t]->0,{m,0,highestcoeff}];
vdashtable = Table[V[2m]'[x]==-acc0list[[m+1]]/.irrelevantrepl,{m,0,highestcoeff}];

LogPrint["Applying Simplifications.","PRG:POT"];
DumpSave["vdash.mx",vdashtable];
ClearAll[vdashtable];
<<vdash.mx;

LogPrint["Solving for V_2m[x].","PRG:POT"];

vsol = DSolve[Flatten[{vdashtable,Table[V[2m][0]==0,{m,0,highestcoeff}]}],Table[V[2m],{m,0,highestcoeff}],x][[1]];


For[n=0,n<highestcoeff+1,n++,
	If[n!=0,<<V.mx];
	V[2n][x_]=TrigReduce[(V[2n]/.vsol[[n+1]])[x]];
	DumpSave["V.mx",V];
	ClearAll[V];
];


PD[];

(*CLEANUP*)

<<V.mx;

Veffsimple[x_] = Veff[x]/.E0->En;

LogPrint["Dumpsaving effective potential.","PRG:POT"];

DumpSave["Veff.mx",{Veff,Veffsimple}];



][[1]];

PD[];LogPrint["Effective potential corrections acquired up to non-trivial power correction "<>ToString[highestpowercorr]<>". This took "<>ToString[vefftime]<>" s.","INFORM"];


ClearAll["Global`*"];

<<ini.mx;<<ser.mx;

<<V.mx;<<Veff.mx;
Vsimplist = Table[Coefficient[Expand[Veffsimple[x]],\[Epsilon],2n],{n,0,highestcoeff}];
DumpSave["Vres.mx",Vsimplist];

][[1]];

PD[];PD[];LogPrint["Potential calculation completed. It took "<>ToString[total]<>" s.","INFORM"];
LogPrint["Access solutions using Results.nb. Match initial conditions by running incond.wls","INFORM"];
PD[];
LogPrint["### POTENTIAL CALCULATION DONE ###","RUNNING"];
PD[];

