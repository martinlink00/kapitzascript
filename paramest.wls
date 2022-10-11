Off[FrontEndObject::notavail];

SetDirectory["exports"]; SetDirectory["logs"];

potloglist = FileNames["pot_*.log"];
incloglist = FileNames["inc_*.log"];

extracttime[file_] := Module[ {},
   found = Find[file, "completed. It took"]; Close[file];
   ToExpression[
    StringCases[found, 
      RegularExpression[
       "([+-]?(?=\\.\\d|\\d)(?:\\d+)?(?:\\.?\\d*))(?:[eE]([+-]?\\d+))?\
"]][[-1]]]
   ];

extractorder[file_] := Module[ {},
   found = Find[file, "non-trivial correction"]; Close[file];
   ToExpression[Characters[found][[-2]]]
   ];

extractmaxram[file_] := Module[ {},
   filetab = Import[file, "Table"];
   ramtab = filetab[[All, 5]];
   Max[ramtab]
   ];

Print["Found " <> ToString[Length@potloglist] <> 
   " potential calculation logfiles and " <> 
   ToString[Length@incloglist] <> 
   " initial condition matching logfiles to extract data from."];
Print["Extracting data."];

potdata = {};
For[i = 1, i <= Length@potloglist, i++,
  file = ToString[potloglist[[i]]];
  AppendTo[
   potdata, {extractorder[file], extracttime[file], 
    extractmaxram[file]}];
  ];
incdata = {};
For[i = 1, i <= Length@incloglist, i++,
  file = ToString[incloglist[[i]]];
  AppendTo[
   incdata, {extractorder[file], extracttime[file], 
    extractmaxram[file]}];
  ];

Print["Fitting data."]

timefitpot[
   z_] = (a*Exp[b*z] /. 
    FindFit[potdata[[All, 1 ;; 2]], a*Exp[b*x], {a, b}, x]);
timefitinc[
   z_] = (a*Exp[b*z] /. 
    FindFit[incdata[[All, 1 ;; 2]], a*Exp[b*x], {a, b}, x]);
timeplotfitf = 
  LogPlot[{timefitpot[x], timefitinc[x]}, {x, 0, 15}, 
   AxesLabel -> {"Non-trivial correction", "Computation time (s)"}, 
   GridLines -> Automatic, 
   PlotLegends -> {"Potential", "Initial conditions"}];
timeplotdat = 
  ListLogPlot[{potdata[[All, 1 ;; 2]], incdata[[All, 1 ;; 2]]}, 
   PlotRange -> {{0, 15}, All}, 
   AxesLabel -> {"Non-trivial correction", "Computation time (s)"}];
exptimeplot = Show[timeplotfitf, timeplotdat];

ramfitpot[z_] = 
  a*Exp[b*z] /. 
   FindFit[potdata[[All, 1 ;; 3 ;; 2]], a*Exp[b*x], {a, b}, x];
ramfitinc[z_] = 
  a*Exp[b*z] /. 
   FindFit[incdata[[All, 1 ;; 3 ;; 2]], a*Exp[b*x], {a, b}, x];
ramplotfitf = 
  LogPlot[{ramfitpot[x], ramfitinc[x]}, {x, 0, 15}, 
   AxesLabel -> {"Non-trivial correction", "Max used RAM (GB)"}, 
   GridLines -> Automatic, 
   PlotLegends -> {"Potential", "Initial conditions"}];
ramplotdat = 
  ListLogPlot[{potdata[[All, 1 ;; 3 ;; 2]], 
    incdata[[All, 1 ;; 3 ;; 2]]}, PlotRange -> {{0, 15}, All}, 
   AxesLabel -> {"Non-trivial correction", "Max used RAM (GB)"}];
expramplot = Show[ramplotfitf, ramplotdat];

SetDirectory[ParentDirectory[]];
SetDirectory[ParentDirectory[]];

hc = ReadList["params.txt"][[1]];

Export["comptime.pdf", exptimeplot]; Export["ram.pdf", expramplot];


Print["Computation time and RAM usage of all available logfiles has \
been fitted and plotted into comptime.pdf and ram.pdf."];
Print["For the current inserted highest correction (" <> 
   ToString[hc] <> ") the following parameters are to be expected:"];
Print["Running potential.wls is expected to take " <> 
   ToString[N[timefitpot[hc], 2]] <> " s and use up " <> 
   ToString[N[ramfitpot[hc], 2]] <> " GB of RAM."];
Print["Running incond.wls is expected to take " <> 
   ToString[N[timefitinc[hc], 2]] <> " s and use up " <> 
   ToString[N[ramfitinc[hc], 2]] <> " GB of RAM."];
