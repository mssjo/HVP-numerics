 (* This Mathematica script should be invoked through series.sh. *)

 (* Settings and definitions *)

Get["series.m"];

EnsureSolution[sol_] := If[
    Length[sol] == 1,
    First[sol],
    Print["Failed to solve differential equation!"]; Exit[1]];

InfiniteQ[v_] := Or[v == Infinity, v == -Infinity];
Options[SolveBySeries] = {point -> 0, logord -> 0, minord -> 0, verbose->False};
SolveBySeries[diffeq_, init_, ord_, offs_, var_, opts:OptionsPattern[]] := (

    VPrint = If[OptionValue[verbose], Print, Null&];

    coeffs = Flatten @ Table[coeff[aa,bb], {aa, OptionValue[minord], ord+1}, {bb, 0, OptionValue[logord]}];
    VPrint["> solving for: ", coeffs];
    ansatz = Sum[var^If[InfiniteQ[OptionValue[point]], -aa, aa] * Log[var]^bb * coeff[aa,bb],
        {aa, OptionValue[minord], ord+1}, {bb, 0, OptionValue[logord]}] /. init;
    VPrint["> ansatz: ", ansatz];

    Clear[solution]; solution[0] = init;
    Do [
        eqn = Normal @ Series[diffeq[ansatz], {var, OptionValue[point], iter+offs}];
        VPrint["> equation iter ", iter, ": 0 = ", eqn];
        solution[n] = EnsureSolution[SolveAlways[eqn == 0, var]];
        VPrint["> solution iter ", iter, ": ", solution[n]];
        ansatz = ansatz /. solution[n];
        ,
        {iter, 1, ord}];
    Series[ansatz, {var, OptionValue[point], ord}]);

tSeries[expr_, var_:t] := Series[expr, {var, 0, order}];
bSeries[expr_, var_:b] := Series[expr, {var, Infinity, 2*order}];

SetAttributes[PrintSeries, HoldAll];
PrintSeries[name_, series_, var_, ref_:None, tag_:None] := (
    fail = False;
    vtag = If[tag =!= None, " v. " <> tag, ""];
    If[And[verify, ref =!= None],
        If[Quiet @ Check[PossibleZeroQ @ Normal @ Simplify[series-ref], 0, {Series::ztest}],
            Print[name, " (\\033[32mPASS", vtag, "\\033[0m)", If[verbose, ": ", ""]];
            ,
            tmp1 = verbose; tmp2 = vtag; verbose = True;
            PrintSeries[name <> " (\\033[33mDIFF" <> vtag <> "\\033[0m)", series-ref, var];
            verbose = tmp1; vtag = tmp2;
            Print[name, " (\\033[31mFAIL", vtag, "\\033[0m): "];
            fail=True;
            ,
            tmp = verbose; verbose = True;
            PrintSeries[name <> " (\\033[34mREFERENCE" <> vtag <> "\\033[0m)", ref, var];
            verbose = tmp;
            Print[name, " (\\033[34mUNDECIDED", vtag, "\\033[0m): "];
            fail=True],
        If[verbose, Print[name, ": "]]];
    If[Or[verbose, fail],
        ord = 0;
        While[
            SeriesCoefficient[series, ord] =!= Indeterminate,
            If[SeriesCoefficient[series, ord] != 0,
                Print[" + ", var, "^", ord, "*",
                    Together @ Simplify @ SeriesCoefficient[series, ord],
                    " ~ ", N @ SeriesCoefficient[series, ord]]
                    ,
                    ""];
            ord++]]);

tExport[name_, expr_:None] := (
    file = "series/" <> name <> ".txt";
    Export[file, Simplify /@ CoefficientList[Normal @ If[expr =!= None, expr, Symbol[name]], t]];
    Print["\\033[90mWrote to ", file, "\\033[0m"]);
bExport[name_, expr_:None] := (
    file = "series/" <> name <> ".txt";
    Export[file, Simplify /@ CoefficientList[Normal @ If[expr =!= None, expr, Symbol[name]] /. {b->1/bi}, bi]];
    Print["\\033[90mWrote to ", file, "\\033[0m"]);

(* Jbub series *)

verify = True (*False*); (* Verification is slow here because it relies on Mathematica's expansions of the logs *)

bJbub[n_, b_] := (
    bp = (b + 1)/(2*b);
    bm = (b - 1)/(2*b);
    f2[z_] := PolyLog[2, z] + Log[z]^2/2 + 2 Log[z];
    f3[z_] := PolyLog[3, z] + Log[z]^3/12 - (Pi^2 * Log[z])/12 + f2[z];
    Switch[n,
        1, b*(Log[bp] - Log[bm]) - 2,
        2, 8 - 2*b*(f2[bp] - f2[bm]),
        3, -48 + 12*b*(f3[bp] - f3[bm]) + 3*b*Log[bp]*Log[bm]*(Log[bp] - Log[bm])]);

bseriesJbub[0] = 1;
bdiffeqJbub[f_, n_] := (1 - b^2)*(b*D[f, b] - f) - 2*n*bseriesJbub[n-1];
Do[
    If[debug, Print["n = ", n]];
    bseriesJbub[n] = SolveBySeries[
        bdiffeqJbub[#,n]&, {coeff[0,0] -> 0},
        2*order, -1, b, point->Infinity,
        verbose->debug];
    PrintSeries["Jbub" <> ToString[n], bseriesJbub[n], "1/b",
        Simplify @ bSeries[bJbub[n,b]], "logs"];
    bExport["bseriesJbub" <> ToString[n], bseriesJbub[n]];
    ,
    {n, 1,3}];

tJbub[n_, b_] := bJbub[n,b] /. b -> Sqrt[1 - 4/t];
tseriesJbub[0] = 1;
tdiffeqJbub[f_, n_] := t*(t - 4)*D[f, t] - 2*f - n*t*tseriesJbub[n-1];
Do[
    If[debug, Print["n = ", n]];
    tseriesJbub[n] = SolveBySeries[
        tdiffeqJbub[#,n]&, {coeff[0,0] -> 0},
        order, 0, t,
        verbose->debug];
    PrintSeries["Jbub" <> ToString[n], tseriesJbub[n], "t",
        tSeries[Normal[bseriesJbub[n]] /. b -> Sqrt[1-4/t]], "bseries"];
    tExport["tseriesJbub" <> ToString[n], tseriesJbub[n]];
    ,
    {n, 1,3}];


(* E(2d) series *)

verify = True;

diffeqE1[f_] := (
    + D[f, {t,3}] * (t-16)*(t-4)*t^2
    + D[f, {t,2}] * 6*t*(t^2 - 15*t + 32)
    + D[f, {t,1}] * (7*t^2 - 68*t + 64)
    + D[f, {t,0}] * (t - 4)
    - 24);
tseriesE1 = SolveBySeries[
    diffeqE1, {coeff[0,0] -> -7*Zeta[3](*, coeff[1,0] -> 3/8 - (7*Zeta[3])/16*)},
    order, 0, t,
    verbose->debug];

PrintSeries["E1", tseriesE1, "t",
    Series[
        + t^0 * (-7*Zeta[3])
        + t^1 * (3/8 - (7*Zeta[3])/16)
        + t^2 * (27/512 - (49*Zeta[3])/1024)
        + t^3 * (37/4608 - (7*Zeta[3])/1024)
        + t^4 * (25555/18874368 - (4753*Zeta[3])/4194304)
        + t^5 * (9304913/37748736000 - (13783*Zeta[3])/67108864)
        ,
        {t,0,5}], "old"];

tExport["tseriesE1"];

bseriesE1 = bSeries[Normal[tseriesE1] /. t -> 4/(1-b^2)];
PrintSeries["E1", bseriesE1, "1/b"];
bExport["bseriesE1"];


tseriesE2 = tSeries[-t/4 * D[tseriesE1, t]  - 1/4 * tseriesE1];
PrintSeries["E2", tseriesE2, "t",
    Series[
        + t^0 * (7*Zeta[3])/4
        + t^1 * (-3/16 + (7*Zeta[3])/32)
        + t^2 * (-81/2048 + (147*Zeta[3])/4096)
        + t^3 * (-37/4608 + (7*Zeta[3])/1024)
        + t^4 * (-127775/75497472 + (23765*Zeta[3])/16777216)
        + t^5 * (-9304913/25165824000 + (41349*Zeta[3])/134217728)
        ,
        {t, 0, 5}], "old"];
tExport["tseriesE2"];

tseriesE3 = (t/2 * D[tseriesE1, {t,2}] + (4+t)/8*D[tseriesE1, t] + 1/8*tseriesE1);
PrintSeries["E3", tseriesE3, "t",
    Series[
        + t^0 * (3/16 - (35*Zeta[3])/32)
        + t^1 * (51/256 - (105*Zeta[3])/512)
        + t^2 * (229/4096 - (399*Zeta[3])/8192)
        + t^3 * (35027/2359296 - (6545*Zeta[3])/524288)
        + t^4 * (3953471/1006632960 - (439635*Zeta[3])/134217728)
        + t^5 * (20893409/20132659200 - (463785*Zeta[3])/536870912)
        ,
        {t, 0, 5}], "old"];
tExport["tseriesE3"];

(* Ebar series *)

tseriesE1bar = tSeries[
    + tseriesE1 * (t^3 + 24*t^2 - 600*t + 896)/288
    - tseriesE2 * (t^4 + 12*t^3 - 792*t^2 + 3968*t + 1536)/288
    - tseriesE3 * (t-16)*(t-4)*(t^2 + 40*t + 64)/144
    + (71*t^2 + 18*t + 6102)/216
    + (23 - t)*Pi^2/12
    - 2*Zeta[3]
    ];
PrintSeries["E1bar", tseriesE1bar, "t",
    Series[
        + t^0 * (275 + 23*Pi^2 - 24*Zeta[3])/12
        + t^1 * (-61 - Pi^2 + 42*Zeta[3])/12
        + t^2 * (17/864 - (7*Zeta[3])/64)
        + t^3 * (106 - 63*Zeta[3])/36864
        + t^4 * (134 - 105*Zeta[3])/983040
        + t^5 * (21086 - 17325*Zeta[3])/1887436800
        ,
        {t, 0, 5}], "old"];
tExport["tseriesE1bar"];

tseriesE2bar = (
    + tseriesE1 * (5*t^2 - 80*t + 96)/96
    - tseriesE2 * (5*t^3 - 116*t^2 + 376*t + 896)/96
    - tseriesE3 * (t-16)*(t-4)*(5*t + 16)/48
    + (95*t + 112 + (28 - t)*Pi^2)/48
    - Zeta[3]
    );
PrintSeries["E2bar", tseriesE2bar, "t",
    Series[
        + t^0 * (-5/3 + (7*Pi^2)/12 - Zeta[3])
        + t^1 * (-7 - Pi^2 + 42*Zeta[3])/48
        + t^2 * -1/48
        + t^3 * (-106 + 63*Zeta[3])/147456
        + t^4 * (-134 + 105*Zeta[3])/1966080
        + t^5 * (-21086 + 17325*Zeta[3])/2516582400
        ,
        {t, 0, 5}], "old"];
tExport["tseriesE2bar"];

tseriesE3bar = (
    + tseriesE1 * (t^2 - 12*t - 8)/96
    - tseriesE2 * (t^3 - 24*t^2 + 88*t + 160)/96
    - tseriesE3 * (t-16)*(t-4)*(t+4)/48
    + (39*t - 232 - 10*Pi^2)/48
    );
PrintSeries["E3bar", tseriesE3bar, "t",
    Series[
        + t^0 * (-140 - 5*Pi^2 + 84*Zeta[3])/24
        + t^1 * (-3*(-2 + 7*Zeta[3]))/64
        + t^2 * (106 - 63*Zeta[3])/6144
        + t^3 * (134 - 105*Zeta[3])/98304
        + t^4 * (21086 - 17325*Zeta[3])/125829120
        + t^5 * (-3*(-273174 + 226625*Zeta[3]))/33554432000
        ,
        {t, 0, 5}], "old"];
tExport["tseriesE3bar"];

tseriesE4bar = (
    + tseriesE1 * (t^2 - 20*t + 160)/96
    - tseriesE2 * (t^2 - 28*t + 120)*(t-16)/96
    - tseriesE3 * (t-4)*(t-16)^2/48
    - tseriesJbub[1] * (59 + 3*Pi^2)/8
    + tseriesJbub[2] * 17/8
    - tseriesJbub[3] / 4
    + (25*t - 136 + 14*Pi^2)/24
    - Zeta[3]
    );
PrintSeries["E4bar", tseriesE4bar, "t",
    Series[
        + t^0 * (-5/3 + (7*Pi^2)/12 - Zeta[3])
        + t^1 * (91 + 3*Pi^2 - 42*Zeta[3])/48
        + t^2 * (230 + 12*Pi^2 - 105*Zeta[3])/1920
        + t^3 * (5358 + 512*Pi^2 - 4165*Zeta[3])/573440
        + t^4 * (40382 + 9216*Pi^2 - 72765*Zeta[3])/61931520
        + t^5 * (-4017442 + 9437184*Pi^2 - 73419885*Zeta[3])/348798320640
        ,
        {t, 0, 5}], "old"];
tExport["tseriesE4bar"];

Srat[b_] := (
            + (68*b^4 - 44*b^2 + 24 -  b^2*(b^2 + 5) Pi^2)/(6*(b^2 - 1)^2)
            + (2*b^2*Zeta[3])/(3*(b^2 - 1)));
SJ[b_] :=   (
            + (b^2*bseriesJbub[3])/(3*(b^2 - 1))
            + (2*bseriesJbub[1]*bseriesJbub[2])/(b^2 - 1)
            - ((b^4 + 5) bseriesJbub[2])/(b^2 - 1)^2
            - (12*bseriesJbub[1]^2)/(b^2 - 1)^2
            + ((Pi^2 - 4)  b^4 + (8 - Pi^2)  b^2 - 52)/(2*(b^2 - 1)^2)*bseriesJbub[1]);

Log1[b_] := Log[(b+1)/(b-1)];
sn[n_, b_] := Switch[n,
    0, 2/3*(b^2*(2*b^2 - 1))/(b^2 - 1)^2,
    1, (2*b)/3*(3*b^4 - 7*b^2 + 3)/(b^2 - 1),
    2, (b^2*(4*b^2 - 3))/2];
gn[n_, b_] := Switch[n,
    1, b^2 - 1,
    2, (b^2 - 1)/4 * Log1[b] - b/2];

bS5 = Simplify @ bSeries[
        + Srat[b]
        + SJ[b]
        + Sum[
            sn[k, b] * D[Normal @ bseriesE1, {b,k}],
            {k, 0, 2}]];

tS5 = tSeries[Normal[bS5] /. {b -> Sqrt[1 - 4/t]}];
tdiffeqE5bar[f_] := (
    + D[f, {t,2}] * 12*(t - 4)^2*t^2
    + D[f, {t,1}] * 24*(t - 3)*(t - 4)*t
    + D[f, {t,0}] * 24*(t - 4)
    - 48*tS5);

(* debug = True; *)
If[debug, PrintSeries["S5", tS5, "t"]];

tseriesE5bar = SolveBySeries[
    tdiffeqE5bar, {coeff[0,0] -> 1/3 + Pi^2/12 - 8*Zeta[3]/3},
    order, 0, t,
    verbose->debug];

PrintSeries["E5bar", tseriesE5bar, "t",
    Series[
        + t^0 * (1/3 + Pi^2/12 - 8*Zeta[3]/3)
        + t^1 * (34 + 8*Pi^2 - 105*Zeta[3])/192
        + t^2 * (2770 + 192*Pi^2 - 3675*Zeta[3])/46080
        + t^3 * (23450 - 464*Pi^2 - 7415*Zeta[3])/1228800
        + t^4 * (4450 - 334*Pi^2 + 1485*Zeta[3])/768000
        ,
        {t, 0, 4}], "old"];
tExport["tseriesE5bar"];

homE5[f_] := (
    + D[f, {b,2}] * b^2
    - D[f, {b,0}] * 2*b^2/(b^2 - 1));
bdiffeqE5bar[f_] := (homE5[f] - bS5);

(*PrintSeries["g1", homE5[gn[1,b]], "1/b",
    0, "diffeq"];
PrintSeries["g2", homE5[gn[2,b]], "1/b",
    0, "diffeq"];*)

If[debug, PrintSeries["S5", bS5, "t"]];

bseriesE5bar = SolveBySeries[
    bdiffeqE5bar, {coeff[0,0] -> 1/3 + Pi^2/12 - 8*Zeta[3]/3, coeff[1,0] -> 0},
    2*order, 1, b, point->Infinity,
    verbose->debug];

PrintSeries["E5bar", bseriesE5bar, "1/b",
    bSeries[Normal[tseriesE5bar] /. t -> 4/(1-b^2)], "tseries"];
bExport["bseriesE5bar"];

tseriesE6bar = (
    - (D[tseriesE5bar, t] + 1/t*tseriesE5bar)
    + (
        + (t-16)*(t-4)/(12*t^2) * t*D[t*D[tseriesE1, t], t]
        + (t-10)/(12*t) * t*D[tseriesE1, t]
        + 1/24 * tseriesE1
        )
    + (-2*tseriesJbub[3] + 6*tseriesJbub[2] + 3*(4-Pi^2)*tseriesJbub[1])/(24*t)
    - (Zeta[3]+5)/(3*t) + Pi^2/(12*t));
PrintSeries["E6bar", tseriesE6bar, "t",
    Series[
        + t^0 * (-1/4-1*Pi^2/16+7*Zeta[3]/8)
        + t^1 * (-Pi^2/96+7*Zeta[3]/32-5/32)
        + t^2 * (-Pi^2/480+427*Zeta[3]/8192-13733/184320)
        + t^3 * (-Pi^2/2240+3269*Zeta[3]/262144-137899/5898240)
        ,
        {t, 0, 3}], "old"];
tExport["tseriesE6bar"];

(* Auxiliaries to E5bar *)

hn[n_, b_] := bSeries @ Sum[
        (-1)^k * D[sn[k, b] * gn[n,b]/b^2, {b, k}],
        {k, 0, 2}];
If[debug,
    PrintSeries["h1", hn[1,b], "1/b",
        bSeries[
            (6 - 11*b^2 + 57*b^4 - 54*b^6)/(3*b^2 - 3*b^4)]];
    PrintSeries["h2", hn[2,b], "1/b",
        bSeries[
            (-2*b (6 + 31*b^2 - 93*b^4 + 54*b^6) + (6 - 17*b^2 + 68*b^4 - 111*b^6 + 54*b^8) Log1[b])/(12*b^2*(-1 + b^2)^2)]]];

Irat[n_, b_] := Integrate[ bSeries[Srat[b] * gn[n, b]/b^2], b];
If[debug,
    PrintSeries["I(1)rat", Irat[1,b], "1/b",
        bSeries[
            (Pi^2-8)/2 * Log1[b] + (24 + (68 + 4*Zeta[3] - Pi^2)*b^2)/(6*b)]];
    PrintSeries["I(2)rat", Irat[2,b], "1/b",
        bSeries[
            ((Pi^2-8)/16 * Log1[b] + (24 + (68 + 4*Zeta[3] - Pi^2)*b^2)/(24*b)) * Log1[b] + 2*Log[(b^2-1)/b^2] + (8 - Pi^2)/(4*(b^2 - 1)) - (68 - Pi^2 + 4*Zeta[3])/12]]];
IJ[n_, b_] := Integrate[ SJ[b] * gn[n, b]/b^2, b];



Hdiv[x_] := bSeries[hn[1,x] * (SeriesCoefficient[bseriesE1, 0] + SeriesCoefficient[bseriesE1, 2]/x^2), x];
Hreg[x_] := bSeries[hn[1,x]*ReplaceAll[Normal @ bseriesE1, {b -> x}] - Hdiv[x], x];
If[debug,
    PrintSeries["Hdiv", Hdiv[b], "1/b"];
    PrintSeries["Hreg", Hreg[b], "1/b"]];

IE[n_, b_] := bSeries @ Switch[n,
    1, Integrate[Hreg[b], b],
    2, Integrate[ bseriesE1 * hn[2,b], b]];
Idiv[b_] := bSeries @ Integrate[Hdiv[b], b];

PrintSeries["I(1)E", IE[1,b], "1/b"];
PrintSeries["I(1)E", IE[2,b], "1/b"];
PrintSeries["Idiv", Idiv[b], "1/b"];

iseriesE5bar = bSeries[
    + sn[2,b]/b^2 * bseriesE1
    - Sum[
        (-1)^n * gn[3-n,b] * (Irat[n,b] + IJ[n,b] + IE[n,b] + If[n == 1, Idiv[b], 0]),
        {n, 1,2}]];

PrintSeries["(4.5)", iseriesE5bar, "1/b",
    bseriesE5bar, "bseries"];

(*Gn[n_,b_] := -(b^2-1)^2/(8*b) * Switch[n,
    1, 2*b,
    2, Normal @ bseriesJbub[1]/2] - gn[n,b] / (4/(1-b^2));

Print @ Normal @ bSeries @ Simplify[(
        + (t-16)*(t-4)/(12*t^2) * t*D[t*D[Normal[tseriesE1], t], t]
        + ((t-10)/(12*t) - (t-16)/(2*t**2)) * t*D[Normal[tseriesE1], t]
        + (1/24 - (t**2 - 28*t + 48)/(3*t**2*(t-4)) - sn[2,b]/b^2/t) * Normal[tseriesE1]
        ) /. t -> 4/(1-b^2)];
Quit[];
iseriesE6bar = bSeries[
    + bSeries @ Simplify[ Normal[
        + (t-16)*(t-4)/(12*t^2) * t*D[t*D[tseriesE1, t], t]
        + ((t-10)/(12*t) - (t-16)/(2*t**2)) * t*D[tseriesE1, t]
        + (1/24 - (t**2 - 28*t + 48)/(3*t**2*(t-4)) - sn[2,b]/b^2/t) * tseriesE1
        ] /. t -> 4/(1-b^2)]
    + (
        -  4 * bseriesJbub[3]
        + 12 * bseriesJbub[2]
        + 6*(4-Pi^2) * bseriesJbub[1]
        ) / (48*t)
    + (Pi^2 - 4*Zeta[3] - 20) / (12*t)
    - Sum[
        (-1)^n * Gn[3-n,b] * (Irat[n,b] + IJ[n,b] + IE[n,b] + If[n == 1, Idiv[b], 0]),
        {n, 1,2}]];


PrintSeries["(4.9-10)", iseriesE6bar, "1/b",
    bseriesE6bar, "bseries"];*)

PrintSeries["(4.20a)", bSeries[gn[2,b]*Irat[1,b] - gn[1,b]*Irat[2,b]], "1/b",
    Series[
        + (Pi^2 - 4*Zeta[3] - 68)/12
        + (4 - Pi^2)/(4*b^2),
        {b, Infinity, 2}], "paper"];
PrintSeries["(4.20b)", bSeries[gn[2,b]*IJ[1,b] - gn[1,b]*IJ[2,b]], "1/b",
    Series[
        - (4 - Pi^2)/(12*b^2),
        {b, Infinity, 2}], "paper"];
PrintSeries["(4.20c)", bSeries[gn[2,b]*IE[1,b] - gn[1,b]*IE[2,b]], "1/b",
    Series[
        - 63*Zeta[3]/10
        + (15863*Zeta[3] - 6642)/(1440*b^2),
        {b, Infinity, 2}], "paper"];
PrintSeries["(4.20d)", bSeries[sn[2,b]/b^2 * bseriesE1 + gn[2,b]*Idiv[b]], "1/b",
    Series[
        + (180 + 119*Zeta[3])/30
        + (4662-12713*Zeta[3])/(1440*b^2),
        {b, Infinity, 2}], "paper"];


bseriesSJg1 = bSeries[SJ[b] * gn[1,b]/b^2];
PrintSeries["SJg1", bseriesSJg1, "1/b",
    Series[
        + 1/b^2  * (Pi^2/3 - 4/3)
        + 1/b^4  * Pi^2/5
        + 1/b^6  * (Pi^2/7 - 352/21)
        + 1/b^8  * (Pi^2/9 - 4624/135)
        + 1/b^10 * (Pi^2/11 - 876368/17325)
        + 1/b^12 * (Pi^2/13 - 12085792/184275)
        + 1/b^14 * (Pi^2/15 - 86635328/1091475)
        + 1/b^16 * (Pi^2/17 - 7406063248/80405325)
        + 1/b^18 * (Pi^2/19 - 84060807056/808782975)
        + 1/b^20 * (Pi^2/21 - 138669411904/1206079875)
        ,
        {b, Infinity, 20}], "old"];
bExport["bseriesSJg1"];

bseriesSJg2 = bSeries[SJ[b] * gn[2,b]/b^2];
PrintSeries["SJg2", bseriesSJg2, "1/b",
    Series[
        + 1/b^5  * (4 - Pi^2)/9
        + 1/b^7  * (8 - 3*Pi^2)/15
        + 1/b^9  * (388/63 - 142*Pi^2/525)
        + 1/b^11 * 2*(-26528 + 465*Pi^2)/2835
        + 1/b^13 * (1550372 - 15215*Pi^2)/40425
        + 1/b^15 * (1971483496/30405375 - (2689*Pi^2)/6435)
        + 1/b^17 * 4*(-1735258237 + 8054725*Pi^2)/70945875
        + 1/b^19 * -4*(-13766949808 + 48874455*Pi^2)/402026625
        ,
        {b, Infinity, 19}], "old"];
bExport["bseriesSJg2"];

bseriesE1h2 = bSeries[bseriesE1 * hn[2,b]];
PrintSeries["E1h2", bseriesE1h2, "1/b",
    Series[
        + 1/b^3  * -63*Zeta[3]/5
        + 1/b^5  * -(486 + 857*Zeta[3])/180
        + 1/b^7  * -(19330 + 28413*Zeta[3])/6720
        + 1/b^9  * -(1021/336 + 75947*Zeta[3]/15840)
        + 1/b^11 * (-22304071946 - 37545403665*Zeta[3])/6642155520
        + 1/b^13 * -(7*(86369562598 + 148521738375*Zeta[3]))/158146560000
        ,
        {b, Infinity, 13}], "old"];
bExport["bseriesE1h2"];

bseriesHreg = bSeries[Hreg[b]];
PrintSeries["Hreg", bseriesHreg, "1/b",
    Series[
        + 1/b^2 * 189*(3*Zeta[3] - 2)/32
        + 1/b^4 * (693*Zeta[3] - 334)/64
        + 1/b^6 * (776223*Zeta[3] - 322346)/73728
        + 1/b^8 * 7*(133565625*Zeta[3] - 39203686)/110592000
        ,
        {b, Infinity, 8}], "old"];
bExport["bseriesHreg"];
