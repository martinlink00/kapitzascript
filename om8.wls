SetDirectory["exports"]; << vdash.mx; << ini.mx;
exp8 = Expand[Part[List @@ (vdashtable[[8 + 1]]), 2]];
V[n_][x_] = 0;
Expand[D[exp8, x] /. x -> 0] >> "om8.txt";