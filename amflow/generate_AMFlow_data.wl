(*load the package*)
current = If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName;

Get["/home/mssjo/Code/amflow/AMFlow.m"];

(*set ibp reducer, could be "Blade", "FiniteFlow+LiteRed", "FIRE+LiteRed" or "Kira"*)
SetReductionOptions["IBPReducer" -> "FIRE+LiteRed"];

(*read results from cached files whenever possible, False by default*)
SetAMFOptions["UseCache" -> True];

(*configuration of the integral family*)
AMFlowInfo["Family"] = HVP3loop;
AMFlowInfo["Loop"] = {l1,l2,l3};
AMFlowInfo["Leg"] = {p1,p2};
AMFlowInfo["Conservation"] = {p2 -> -p1};
AMFlowInfo["Replacement"] = {p1^2 -> m2*t};
AMFlowInfo["Propagator"] = {
    l1^2-m2, l2^2-m2, l3^2-m2,
    (l1-p1)^2-m2, (l2-p1)^2-m2, (l3-p1)^2-m2,
    (l1+l3)^2-m2, (l2+l3)^2-m2, (l1+l2)^2-m2};
AMFlowInfo["NThread"] = 8;

(*SolveIntegrals: computes given integrals with given precision goal up to given eps order*)
masters2d = {j[HVP3loop, 1,0,0,0,1,0,1,1,0],
             j[HVP3loop, 2,0,0,0,1,0,1,1,0],
             j[HVP3loop, 3,0,0,0,1,0,1,1,0]};
masters4d = {j[HVP3loop, 1,1,0,1,0,0,1,1,0],
             j[HVP3loop, 1,1,0,1,1,0,1,1,0],
             j[HVP3loop, 2,1,0,1,1,0,1,1,0]};

precision = 20;
epsorder = 6;

reim[x_] := SeriesCoefficient[Series[x, {eps,0,0}], 0] // {Re[#], Im[#]}&

data = Table[
    AMFlowInfo["Numeric"] = {t -> tt, m2 -> 1};
    data2d = AbsoluteTiming[
            SetAMFOptions["D0" -> 2];
            sol = SolveIntegrals[masters2d, precision, epsorder];
            masters2d /. sol];
    data4d = AbsoluteTiming[
            SetAMFOptions["D0" -> 4];
            sol = SolveIntegrals[masters4d, precision, epsorder];
            masters4d /. sol];
    Flatten @ Join[{tt}, reim /@ data2d[[2]], reim /@ data4d[[2]], {data2d[[1]], data4d[[1]]}]
    ,
    {tt, {-6,-2,2,6,10,14,18,22}}];

Export["AMFlow_data.dat", data];
Quit[];
