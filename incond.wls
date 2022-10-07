#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearSystemCache[];
Off[FrontEndObject::notavail];



ClearAll["Global`*"];

If[Not[DirectoryQ["exports"]],
  Print["Error! Need to run potential.wls successfully before matching initial conditions."];Quit[];,SetDirectory["exports"]];




<<ini.mx;<<ser.mx;

logfilename = "inc_"<>DateString["ISODateTime"]<>".log";

If[$VersionNumber!=12.1,LogPrint["Version error: Script was developed on and for version 12.1.","ERROR"];Quit[],PD[];LogPrint["### INITIAL CONDITION SCRIPT ###","RUNNING"];PD[];]


userinfo = "Script has been executed on "<>ToString[$OperatingSystem]<>" by the user " <>ToString[$UserName]<>".";
versioninfo = "Running on Mathematica "<>ToString[$VersionNumber]<>".";
processorinfo = ToString[$ProcessorCount]<>" processors of type "<>ToString[$ProcessorType]<>" available.";
memoryinfo = ToString[N[MemoryAvailable[]/10^6]]<>" MB of memory available for allocation.";
paraminfo = "Performing matching process for initial conditions up to non-trivial correction "<>ToString[highestcoeff]<>" .";
exportinfo = "Results will be dumpsaved into the subfolder exports and can be accessed using Results.nb.";
loginfo = "Log for this script execution is saved under exports/logs/"<>logfilename<>".";

LogPrint["### SYSTEM INFORMATION ###","RUNNING"];
LogPrint[userinfo,"RUN:SYS"];LogPrint[versioninfo,"RUN:SYS"];LogPrint[processorinfo,"RUN:SYS"];LogPrint[memoryinfo,"RUN:SYS"];PD[];
LogPrint[paraminfo,"RUN:PRM"];LogPrint[loginfo,"INFORM"];LogPrint[exportinfo,"INFORM"];PD[];
LogPrint["### NOW RUNNING SCRIPT ###","RUNNING"];PD[];

incontime = 
 	Timing[
   
   
   	For[m = 0, m < highestcoeff + 1, m++,
    		For[n = 0, n < \[Epsilon]max - 1, n++,	
				\!\(\*SuperscriptBox[\(\[CurlyPhi][2  m]\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", "+", "n"}], ")"}],
Derivative],
MultilineFunction->None]\)[0] = D[a[2 m][0], {t, n}];
     		];
    	];
   
   	LogPrint["Acquiring equations.", "PRG:INC"];
   
   	<< V.mx;
   
   	\[Theta]inbef = 
    Table[If[
      m == 0, \[Theta]i == \[Theta]sercoeff[[2 m + 1]] /. x -> 1 /. 
       t -> 0, 0 == \[Theta]sercoeff[[2 m + 1]] /. x -> 1 /. 
       t -> 0], {m, 0, highestcoeff}];
   	\[Omega]inbef = 
    Table[If[
      m == 0, \[Omega]i == \[Omega]sercoeff[[2 m + 1]] /. x -> 1 /. 
       t -> 0, 0 == \[Omega]sercoeff[[2 m + 1]] /. x -> 1 /. 
       t -> 0], {m, 0, highestcoeff}];
   	Einser = 
    Series[V[0][
        Sum[\[Epsilon]^(2 m) \[CurlyPhi][2 m][0], {m, 0, 
          highestcoeff}]] + 
       1/2 (Sum[\[Epsilon]^(2 m) \[CurlyPhi][2 m]'[0], {m, 0, 
            highestcoeff}])^2, {\[Epsilon], 0, 2 highestcoeff}][[3]];
   
   
   	<< phi.mx; << a.mx;
   
   	LogPrint["Solving for varphi_i^(2m).", "PRG:INC"];
   	For[m = 0, m < highestcoeff + 1, m++,
    		\[CurlyPhi]in[2 m] = \[CurlyPhi][2 m][0] = 
       Evaluate[\[CurlyPhi][2 m][0] /. 
         Solve[\[Theta]inbef[[m + 1]], \[CurlyPhi][2 m][0]][[1]]];
    	];
   
   	LogPrint["Solving for omega_i^(2m).", "PRG:INC"];
   	For[m = 0, m < highestcoeff + 1, m++,
    		\[CurlyPhi]in[2 m] = Evaluate[\[CurlyPhi]in[2 m]];
    		\[Omega]in[2 m] = \[CurlyPhi][2 m]'[0] = 
      Evaluate[\[CurlyPhi][2 m]'[0] /. 
        Solve[\[Omega]inbef[[m + 1]], \[CurlyPhi][2 m]'[0]][[1]]];
    	];
   
   	LogPrint["Declaring E0.", "PRG:INC"];
   	E0 = Sum[
     Factor[TrigReduce[Einser[[2 m + 1]]]] \[Epsilon]^(2 m), {m, 0, 
      highestcoeff}];
   
   	DumpSave["incond.mx", {\[CurlyPhi]in, \[Omega]in, E0}];
   
   	LogPrint["Inserting E0 into V_eff.", "PRG:INC"];
   	<< Veff.mx; << incond.mx; << V.mx;
   	LogPrint["Expanding.", "PRG:INC"];
   	veffexp = Series[Veff[x], {\[Epsilon], 0, 2 highestcoeff}];
   	Vreslist = 
    Table[Coefficient[veffexp, \[Epsilon], 2 n], {n, 0, highestcoeff}];
   	Vsimplist = 
    Table[Coefficient[Expand[Veffsimple[x]], \[Epsilon], 2 n], {n, 0, 
      highestcoeff}];
   	DumpSave["Vres.mx", {Vreslist, Vsimplist}];
   
   
   	][[1]];

	PD[]; LogPrint[
 "Initial condition matching completed. It took " <> 
  ToString[incontime] <> " s.", "INFORM"];

 LogPrint["Access solutions using Results.nb.","INFORM"];
 LogPrint["Warning: IF YOU WISH TO KEEP THE RESULTS, COPY THE 'exports' DIRECTORY AND 'Results.nb' SOMEWHERE ELSE NOW", "INFORM"];
PD[];
LogPrint["### INITIAL CONDITION MATCHING DONE ###","RUNNING"];
PD[];
