(* Import AMFlow.m *)
<< AMFlow`;

(* Set IBP reducer, options are "Blade", "FiniteFlow+LiteRed", "FIRE+LiteRed" or "Kira"*)
SetReductionOptions["IBPReducer" -> "FIRE+LiteRed"];

(* Read results from cached files whenever possible *)
SetAMFOptions["UseCache" -> True];

(* Define integral family *)
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

(* Function to compute integral *)
Compute[integral_, tt_, prec_] := (
    AMFlowInfo["Numeric"] = {t -> tt, m2 -> 1};
    sol = SolveIntegrals[{integral}, prec, 6];
    Export["integral.dat", integral /. sol]);

