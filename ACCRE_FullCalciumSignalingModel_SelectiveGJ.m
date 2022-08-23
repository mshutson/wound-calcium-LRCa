(* ::Package:: *)

(* ::Input::Initialization:: *)
(*Parameters that will be varied*)
extParams={r\[Mu]T,\[Tau]Heal,rPMCA,kPMCA,rSOC,\[Eta]NCX,nNCX,kNCX,rlkPM,\[Alpha],\[Alpha]0,kl,nl,Kc,\[Eta]IPR,\[Eta]SERCA,kSERCA1,kSERCA2,\[Eta]lkER,Bx,\[Eta]GJIP3,\[Eta]GJc,connectAblated};
(*Parameters that will not be varied*)
intParams=
{
kdeg->1.25 (*s^-1*),
(*\[Delta]\[Rule]1.234*10^-3 (*dimensionless; one of 3 options given by Lemon (also rr)*),*)
\[Nu]->4.3*10^-9*15*10^-5 (*dm^3;Assumes a circular cell surface of 4.3*10^-9 dm^2 cell height of 15 \[Mu]m)*),
a0->4.3*10^-9(*dm^2*),

\[Tau]d->.001 (*s; really fast*),
nPMCA->2. (*dimensionless; chosen to match nSERCA*),
(*kPMCA\[Rule].45 (*\[Mu]M; matching Han et al. 2017, although various values in various published models*),*)
cExt->1000. (*\[Mu]M; 1 mM*),
nSOC->3.8 (*dimensionless; from "Recent developments in models of calcium signalling"*),
kSOC->187. (*\[Mu]M (per ER Volume); from "Recent developments in models of calcium signalling"*),
(*\[Eta]lkPM\[Rule].0000019(*s^-1. Obtained from constant flux from Han et al 2017 where I have assumed 1mM external calcium concentration*),*)
(*nNCX\[Rule]1. (*Just something to use for now*),*)

\[Epsilon]->1./.185 (*Dimensionless; ratio between cyt volume and ER volume; various values used in published models*),
nSERCA->2. (*dimensionless*),
(*kSERCA\[Rule].1 (*\[Mu]M various values in various published models*),*)
d1->.13 (*\[Mu]M*),
d2->1.05 (*\[Mu]M*),
d3->0.943 (*\[Mu]M*),
d5->0.0823 (*\[Mu]M*),
a2->0.2(*\[Mu]M^-1s^-1*),

Be->150.(*\[Mu]M*),
Ke->10.(*\[Mu]M*),
Kx->.167(*\[Mu]M. Taken from Chen et al.2014*),
nx->2.96(*dimensionless; Taken from Chen et al. 2014*)
};

woundSizeScale=1.;
gbpScale=1.;
\[Rho]r=lr;
rh=\[Alpha]*(c[t]/(Kc+c[t]))*(\[Rho]r^nl/(kl^nl+\[Rho]r^nl)+\[Alpha]0);


dip3=ip3'[t]==rh-kdeg*ip3[t]+jGJIP3;

(*Time dependence of \[Mu]T damage*)
fHeal=Exp[-(t/\[Tau]Heal)]*(1-Exp[-(t/\[Tau]d)])*UnitStep[t];

(*PMCA Flux*)
jPMCA=rPMCA*c[t]^nPMCA/(kPMCA^nPMCA+c[t]^nPMCA);

(*NCX Flux*)
jNCX=c[t]^nNCX/(kNCX^nNCX+c[t]^nNCX);

(*Leak current across the PM.*)
jlkPM=rlkPM*(1-c[t]/cExt);

(*microtear flux for damaged cells. Add in wound size scale?*)
j\[Mu]T=woundSizeScale*r\[Mu]T*(1-c[t]/cExt)*fHeal;

(*SOC flux*)
jSOC=rSOC*kSOC^nSOC/(kSOC^nSOC+cER[t]^nSOC);

(*Full PM Flux*)
jPM=\[Eta]NCX*(j\[Mu]T+jlkPM+jSOC-jPMCA-jNCX);
(*IPR flux terms*)
\[Zeta]=d2*(ip3[t]+d1)/(ip3[t]+d3);

\[Tau]h=1/(a2*(\[Zeta]+c[t]));

h\[Infinity]=\[Zeta]/(\[Zeta]+c[t]);

m\[Infinity]=(ip3[t]/(d1+ip3[t]))*(c[t]/(d5+c[t]));

dh=h'[t]==(h\[Infinity]-h[t])/\[Tau]h;
(*IPR flux*)
jIPR=\[Eta]IPR*m\[Infinity]^3*h[t]^3*(cER[t]-c[t]);

(*leak flux*)
jlkER=\[Eta]lkER*(cER[t]-c[t]);

(*SERCA flux*)
jSERCA=\[Eta]SERCA*(c[t]^nSERCA-kSERCA1^nSERCA*kSERCA2^nSERCA*cER[t]^nSERCA)/(kSERCA1^nSERCA+c[t]^nSERCA);

(*Total ER Flux*)
jER=jIPR+jlkER-jSERCA;
(*Rapid Buffer Terms*)
\[Beta]=(1+(Ke*Be)/(Ke+c[t])^2+(nx*Kx^nx*Bx*c[t]^(nx-1))/(Kx^nx+c[t]^nx)^2)^-1; (*Modified to account for Hill kinetics for general n. Derived by me with cooperative binding and rapid buffer assumption
*)

dc=c'[t]==\[Beta]*(jER+jPM+jGJc);

dcER=cER'[t]==\[Epsilon]*-jER;
(*Just values to start with. When initial values are actually important we will just run the unwounded model for a long time.*)
cEQ=.1;
cEREQ=200.;
ip3EQ=.01;
hEQ=.67;

initc=c[0]==cEQ;
initcER=cER[0]==cEREQ;
initIP3=ip3[0]==ip3EQ;
inith=h[0]==hEQ;
deqns={dc,dcER,dip3,dh,initc,initcER,initIP3,inith};
c0=.1;
minInt=0.09971653246990764`;
maxInt=0.14249861653280668`;
(*This function does not work for fluorescence values larger than maxInt*)
FtoC[f_]:=(((maxInt-minInt)*c0^nx+(f-minInt)*Kx^nx)/(maxInt-f))^(1/nx)/.intParams
CtoF[c_]:=((c^nx-c0^nx)*maxInt+(c0^nx+Kx^nx)*minInt)/(c^nx+Kx^nx)/.intParams
tLong=1000.;
newInits=
ParametricNDSolveValue[
deqns/.intParams/.{r\[Mu]T->0.,lr->0.,jGJc->0.,jGJIP3->0.},
Table[x[tLong],{x,{c,cER,ip3,h}}],
{t,0,tLong},
extParams];
(*produces a set of rules to replace the default initial conditions.*)
getNewInits[params_?AssociationQ]:=Table[({c,cER,ip3,h}[[i]][0]==val_)->({c,cER,ip3,h}[[i]][0]==#[[i]]),{i,Length@#}]&@(newInits@@(params/@extParams))
(*video=Import["C:\\Users\\Greg Stevens\\Box\\Research Backup\\Current Work\\NewCalciumModel\\Flares\\Correlations\\Videos\\1-16-15_GCaMP6m_speed2_gain90_zoom1_40x_movie.tif"];*)
(*Avg cell diamater in \[Mu]m*)
\[CapitalDelta]x=7.3;
(*Manipulate[
{Show[
ImageAdjust[video[[i]],{c,b}],
Graphics[{Yellow,Circle[{246,258},rad],Red,Point[{246,258}]}]
],

Grid[{{"Pixels","\[Mu]m","cells (rounded)","time (s)"},Round[#,.01]&/@{rad,rad/1.6,Round[rad/(1.6*\[CapitalDelta]x)]+1,2.14*(i-67)}},Frame\[Rule]All]
}

,{i,1,360,1}
,{c,0,1}
,{b,0,1}
,{{rad,50},5,256}
]*)
(*Ring info pulled from Manipuate above*)
distanceToRing[x_]:=Round[x/\[CapitalDelta]x]+1
ringToDistance[r_]:=\[CapitalDelta]x*(r-1)

cavRadMicrons=51.25*woundSizeScale;
ringsDamaged=distanceToRing[cavRadMicrons];
ablatedRadius=24.5*woundSizeScale;
ringsAblated=distanceToRing[ablatedRadius];
numRingsFirstExp=10;
numRingsMax=21;
maxRad2ndExp=111.;
maxRadFlares=160.;

(*Pull from O'Connor et al Zones of Damage paper. Assume "high damage" region is nuclear envelope breakdown?*)
zoneEqn=rPM==1.1*rNM+14;
radiusHighDamage=rNM/.First@Solve[zoneEqn/.rPM->cavRadMicrons,rNM];
ringsHighDamage=Round[radiusHighDamage/\[CapitalDelta]x]+1;
\[Mu]TDistData=Import["microtearDistribution_Andrew_dataPoints.xlsx"]//First;
\[Mu]TDistDataScaled=MapAt[Rescale[#,{Min@#,Max@#}&@\[Mu]TDistData[[;;,2]]]&,\[Mu]TDistData,{;;,2}];
(*Try fitting the part after the initial rise to a sigmoid*)
(*Manipulate[ListLinePlot[\[Mu]TDistDataScaled[[i;;]],PlotRange\[Rule]{0,Automatic}],{i,1,Length@\[Mu]TDistData,1}]*)
formFrame=10;
\[Mu]TFit=NonlinearModelFit[\[Mu]TDistDataScaled[[formFrame;;]],{a*Erf[-(x-x0)/\[Chi]]+b,\[Chi]>0,a>0},{{\[Chi],10},{x0,46},{a,.5},{b,.5}},x];
\[Mu]TDistParams=\[Mu]TFit["BestFitParameters"];
(*Find fitting function value at first 0 point in Andrew's data*)
cavRadPixels=(#[[FirstPosition[#[[formFrame;;]],{x_,y_}/;y<.01]+formFrame-1]]&@\[Mu]TDistDataScaled)[[1,1]];
(*Right now just do a horizontal dilation / conraction to "convert" between pixel values and microns. Use the cavitation radius from Andrew's data (pixel distance where the function falls below some arbitrary threshold) and the video (radius of the initial calcium influx)*)
\[Mu]TDist[x_]:=If[x<=cavRadMicrons,(a*Erf[-(x*(cavRadPixels/cavRadMicrons)-x0)/\[Chi]]+b)/.\[Mu]TDistParams,0.];
numRings=numRingsMax;

(*Denote ring position as the middle of the ring. Note that the first ring's width is only half a cell, since we are defining the first ring to be a single cell*)
ringPositions=Range[0.,(numRings-1)*\[CapitalDelta]x,\[CapitalDelta]x];

(*Get \[Sigma]\[Mu]T for each cell. If ring is not damaged, set to 0. Ablated cells will be taken care of when equations are made*)
\[Sigma]\[Mu]TList=ReplacePart[\[Mu]TDist/@ringPositions,Table[i->0.,{i,ringsDamaged+1,numRings}]];
(*This is so that we can relate the \[Mu]T parameter from the single-cell DKO model to the actual maximum parameter value*)
Interpolation[\[Mu]TDistData/.{t_,f_}->{t,f/Max@\[Mu]TDistData[[;;,2]]}];
\[Mu]TConversion=%[cavRadMicrons*.96];
rd=10.;(* \[Mu]m. spatial scale.*)
t0=47; (*s. median time for 2nd expansion signal to start from Erica's paper. *)
DL=260; (*\[Mu]m^2/s. Estimated free diffusion constant of GBP*)
LRth =.5; (*Dimensionless. Ratio of [L.R] Subscript[to [R], T] needed for signal to be "on"*)
dr=.001; (*"small" r since we cannot use the origin in polar coordinates, \[Mu]m*)
\[Rho]Far=1000/rd (*r that is "infinity" for numeric solving. scaled*);

(*Differential Equations*)
s=\[Gamma]d/\[Gamma]w^2*Exp[-(r-dr)^2/(2*\[Gamma]w^2)-\[Gamma]d*t];
dx=D[x[r,t],t]==(1+(\[Gamma]pL*yT[r,t])/(1+\[Gamma]P*x[r,t])^2)^-1*(\[Gamma]DP*Laplacian[x[r,t],{r,\[Theta]},"Polar"]+s+\[Gamma]c*\[Gamma]pL/\[Gamma]P*((\[Gamma]P*x[r,t])/(1+\[Gamma]P*x[r,t]))^2*yT[r,t]);
dyT=D[yT[r,t],t]==(-\[Gamma]c*\[Gamma]P*x[r,t]*yT[r,t])/(1+\[Gamma]P*x[r,t]);
dz=D[z[r,t],t]==(1+\[Gamma]R/(1+\[Gamma]L*z[r,t])^2)^-1*(\[Gamma]DL*Laplacian[z[r,t],{r,\[Theta]},"Polar"]+(\[Gamma]c*\[Gamma]P*x[r,t]*yT[r,t])/(1+\[Gamma]P*x[r,t]));

(*Initial Conditions and Boundary Conditions*)
initx=x[r,0]==0.; (*Initially no active protease*)
bc1x=(D[x[r,t],r]==0)/.r->dr; (*Radially symmetric --> no flux through origin*)
bc2x=x[\[Rho]Far,t]==0.; (*Concentration to 0 far away*)

inityT=yT[r,0]==1.; (*proligand initially exists at all points in space*)

initz=z[r,0]==0; (*Initially no ligand*)
bc1z=(D[z[r,t],r]==0)/.r->dr; (*Radially symmetric --> no flux through origin*)
bc2z=z[\[Rho]Far,t]==0.; (*Concentration to 0 far away*)

(*Put it all together*)
deqnsXY={dx,dyT,initx,bc1x,bc2x,inityT};
deqnsZ={dz,initz,bc1z,bc2z};
deqnsLR=Join[deqnsXY,deqnsZ];
(*lrParams=Import["C:\\Users\\Greg Stevens\\Dropbox\\Research_Backup\\CurrentWork\\NewCalciumModel\\Second Expansion\\LRParams.m"];*)
lrParams=Import["LRParams_lowthresh.m"];
lrParamsChosen=lrParams[[1,1]];
lrParamsChosenAll1=lrParams[[1]];
lrEndTime=10.*60./t0;
tEnd=2*lrEndTime*t0;
ligwoundSizeScale=.9;
ligandSol=NDSolveValue[deqnsLR/.{\[Gamma]w->\[Gamma]w*ligwoundSizeScale,\[Gamma]P->\[Gamma]P*ligwoundSizeScale^2,\[Gamma]pL->\[Gamma]pL*gbpScale,\[Gamma]L->\[Gamma]L*gbpScale}/.lrParamsChosen,z,{r,dr,\[Rho]Far},{t,0,tEnd},Method->{"IndexReduction"->Automatic,"EquationSimplification"->"Residual","PDEDiscretization"->{"MethodOfLines","SpatialDiscretization"->{"TensorProductGrid","MinPoints"->200,"MaxPoints"->500,"DifferenceOrder"->2}}}];
(*LR OUTPUT TO GO INTO CALCIUM MODEL*)
lrSol[r_?NumericQ,t_?NumericQ]:=((\[Gamma]L*gbpScale*ligandSol[r/rd,t/t0])/(1+\[Gamma]L*gbpScale*ligandSol[r/rd,t/t0]))/.lrParamsChosen
caRadData=Import["controlCaRadData.m"];
\[CapitalDelta]t=caRadData[[2,1]]-caRadData[[1,1]];
maxFirstExpRad=Max@caRadData[[;;20,2]];
(*connectAblated=0.015; (*parameter to connect the ablated region to the rest of the tissue. The range can be 1 (full connection) to 0 (no connection)*)*)
(*connectAblated=0.*)

finalParamAssociation=<|r\[Mu]T->41.14497773461118`,\[Tau]Heal->5.25`,rPMCA->0.12`,kPMCA->0.1`,rSOC->0.15000000000000002`,\[Eta]NCX->10.`,nNCX->2.5`,kNCX->1.6`,rlkPM->0.00025`,\[Alpha]->0.9`,\[Alpha]0->0.1`,kl->0.05`,nl->5.`,Kc->0.4`,\[Eta]IPR->4.`,\[Eta]SERCA->5.`,kSERCA1->0.1`,kSERCA2->0.`,\[Eta]lkER->0.01`,Bx->5.`,\[Eta]GJIP3->4.`,\[Eta]GJc->2.`,connectAblated->0.025`|>;
(*List of strings to specify model components to knockout by 70% in one half of the tissue. Current options are GJ, PLC, SERCA, SOC, PMCA. Note that GJs are handled separately when GJ terms are made, since the GJ parameters are specific to cell-cell boundaries rather than single cells.*)

knockoutSpec={"GJIP3","GJc"}[[ToExpression[Environment["SLURM_ARRAY_TASK_ID"]]]]

knockoutVal=0.
knockoutList={knockoutSpec};

(*knockout commands that do not require user input*)
knockoutAssociations=
<|"GJ"->Nothing,
"PLC"->({\[Alpha]->#*\[Alpha](*,\[Alpha]0\[Rule]#*\[Alpha]0*)}),
"SERCA"->(\[Eta]SERCA->#*\[Eta]SERCA),
"SOC"->(\[Eta]SOC->#*\[Eta]SOC),
"PMCA"->(\[Eta]PMCA->#*\[Eta]PMCA),
"GJIP3"->(\[Eta]GJIP3->#*\[Eta]GJIP3),
"GJc"->(\[Eta]GJc->#*\[Eta]CJc)|>&@knockoutVal;

knockoutRules=Flatten[knockoutAssociations/@knockoutList];


(*List of strings to specify model components to overexpress by 100% in one half of the tissue. Current options are PLC.*)
overExpressVal=2.(*overexpression factor*);
overExpressList={};

(*knockout commands that do not require user input*)
overExpressAssociations=<|"PLC"->({\[Alpha]->#*\[Alpha](*,\[Alpha]0\[Rule]#*\[Alpha]0*)})|>&@overExpressVal;

overExpressRules=Flatten[overExpressAssociations/@overExpressList];

gjRand=True; (*True = give each cell a random GJ parameter the tissue model. Value for a border determined by smallest value of the two cells that share that border*)
\[Sigma]GJ=.4;(*\[Sigma] parameter for GJ randomness*)
gjPolarize=False; (*True = multiply GJ params by sin\[Theta], where \[Theta] is the angle between the radial vector to a cell edge and the direction of the cell edge*)

gcampRand=True; (*True = randomize the amount of GCaMP in each cell. The random scale will also scale the fluorescence values of the model output*)
\[Sigma]gcamp=.1;(*\[Sigma] parameter for gcamp randomizer*)

plcRand=False; (*True = randomize PLC (parameter \[Alpha]) between cells*)
varPLC=.225; (*standard deviation for PLC randomizer*)

(*connectAblated=0.015; (*parameter to connect the ablated region to the rest of the tissue. The range can be 1 (full connection) to 0 (no connection)*)*)
connectAblated=.;

mesh=Import["mesh_4826cells_450x450.m"]
meshRatio=1.;

(*Get information about each cell*)
range=450.;(*From mesh file name (side length of grid in \[Mu]m)*)

(*Get distances from center of mesh to each cell centroid*)
distances=Sqrt[(#[[1]]-range/2)^2+(#[[2]]-range/2)^2]&/@PropertyValue[{mesh,2},MeshCellCentroid]*meshRatio;
positions=PropertyValue[{mesh,2},MeshCellCentroid];

(*Get number of cells*)
numCells2D=Length@distances;
interiorCells=MeshCellIndex[mesh,{2,"Interior"}][[;;,2]];
areas=Area[#]&/@MeshPrimitives[mesh,2]*(10^-5)^2*meshRatio^2; (*Convert from \[Mu]m^2 to dm^2*)
avgArea=Mean[areas[[interiorCells]]];
\[CapitalDelta]x=N[Sqrt[avgArea*(10^5)^2/\[Pi]]]*2 ;(*Should be around 7.4*)
(*Determine cells "labels"*)
ablatedPos=Position[distances,x_/;x<ablatedRadius]//Flatten;
woundPos=Position[distances,x_/;ablatedRadius<=x<cavRadMicrons]//Flatten;
highDamagePos=Position[distances,x_/;ablatedRadius<=x<radiusHighDamage]//Flatten;

(*Produce adjacency matrix with shared edge lengths as elements*)
Block[{conn,lengths,directions,positionVectors,sines,edges,adj},

(*Get connectivity matrix, which is a matrix whose rows represent a cell and column represents an edge associated with the cell*)
conn=mesh["ConnectivityMatrix"[2,1]]//Unitize//Normal;

(*Lengths of each edge*)
lengths=EuclideanDistance[#[[1]],#[[2]]]&@@@MeshPrimitives[mesh,1];

(*direction vectors of each edge*)
directions=#[[2]]-#[[1]]&@@@MeshPrimitives[mesh,1];

(*position vectors of the midpoints of each edge*)
positionVectors=((#[[1]]+#[[2]])/2-{range/2,range/2})&@@@MeshPrimitives[mesh,1];

(*List of Sin of angle between position vector of edge midpoint (relative to the wound center) and the direciton of the shared edge*)
sines=Table[Norm@Cross[Append[directions[[i]],0],Append[positionVectors[[i]],0]]/(Norm[directions[[i]]]*Norm[positionVectors[[i]]]),{i,Length@directions}];

(*List where each element is an edge and tells you which cell(s) have that edge*)
edges=(Position[#,1]&/@Transpose[conn])/.{e_}->e;

(*Initialize matrix*)
adj=ConstantArray[0.,{numCells2D,numCells2D}];

(*Function that replaces corresponding elements in the initial adjMat with an edge length if the cells in the row&column share that edge*)
alter[pos_,edge_]:=If[Length[pos]==2,adj[[pos[[1]],pos[[2]]]]=adj[[pos[[2]],pos[[1]]]]=lengths[[edge]]*If[gjPolarize,sines[[edge]],1.]];

MapIndexed[alter[#1,First[#2]]&,edges];

adjMat=SparseArray[adj]*meshRatio;
]
(*(*Check to make sure adjacencies are correct*)
Manipulate[HighlightMesh[HighlightMesh[mesh,Join[Thread[{2,Position[Normal[adjMat[[i]]],x_/;x\[NotEqual]0]//Flatten}]]],Style[{2,i},Red]],{i,1,numCells2D,1}]*)
(*random GJ values. Make sure none of them are negative*)
(*randGJmults=If[gjRand,Abs@RandomVariate[NormalDistribution[1.,\[Sigma]GJ],numCells2D],ConstantArray[1.,numCells2D]];*)
randGJmults=If[gjRand,RandomVariate[LogNormalDistribution[-\[Sigma]GJ^2/2,\[Sigma]GJ],numCells2D],ConstantArray[1.,numCells2D]];

(*\[Eta] values are based on a single "ideal" cell of diameter \[CapitalDelta]x, with flux through its entire perimeter. Therefore, in order to properly scale based on cell volume and shared perimeter between two cells, the GJ fluxes have to be scaled by (v0/Subscript[v, i])*(aij/a0), where the 0 subscript denotes the "ideal" cell, i is the cell in question, and j denotes any cell adjacent to i. Since all cells have the same height, volume ratios are reduced to ratios of PM area, and shared areas are reduced to lengths. So v0 becomes the PM area of the ideal cell, and a0 becomes the circumference of the ideal cell.

a0 and areas are in dm^2, adjMat lengths and \[CapitalDelta]x are in \[Mu]m. So the ratios work out to be dimensionless*)
fluxGJCaTerms2D=Table[\[Eta]GJc*Total[If[MemberQ[knockoutList,"GJ"]&&(positions[[i,1]]<range/2.||positions[[#,1]]<range/2.),knockoutVal,1.]*If[MemberQ[ablatedPos,i]||MemberQ[ablatedPos,#],connectAblated,1.]*(Min[randGJmults[[i]],randGJmults[[#]]]*(a0/areas[[i]])*(adjMat[[i,#]]/(\[Pi]*\[CapitalDelta]x))*(Subscript[c, #][t]-Subscript[c, i][t]))&/@(Flatten@adjMat[[i]]["NonzeroPositions"])],{i,numCells2D}];
fluxIP3Terms2D=Table[\[Eta]GJIP3*Total[If[MemberQ[knockoutList,"GJ"]&&(positions[[i,1]]<range/2.||positions[[#,1]]<range/2.),knockoutVal,1.]*If[MemberQ[ablatedPos,i]||MemberQ[ablatedPos,#],connectAblated,1.]*(Min[randGJmults[[i]],randGJmults[[#]]]*(a0/areas[[i]])*(adjMat[[i,#]]/(\[Pi]*\[CapitalDelta]x))(Subscript[ip3, #][t]-Subscript[ip3, i][t]))&/@(Flatten@adjMat[[i]]["NonzeroPositions"])],{i,numCells2D}];

(*Unwounded terms; only difference is the part with ablation connection*)
fluxGJCaTerms2DUW=Table[\[Eta]GJc*Total[If[MemberQ[knockoutList,"GJ"]&&(positions[[i,1]]<range/2.||positions[[#,1]]<range/2.),knockoutVal,1.]*(Min[randGJmults[[i]],randGJmults[[#]]]*(a0/areas[[i]])*(adjMat[[i,#]]/(\[Pi]*\[CapitalDelta]x))*(Subscript[c, #][t]-Subscript[c, i][t]))&/@(Flatten@adjMat[[i]]["NonzeroPositions"])],{i,numCells2D}];
fluxIP3Terms2DUW=Table[\[Eta]GJIP3*Total[If[MemberQ[knockoutList,"GJ"]&&(positions[[i,1]]<range/2.||positions[[#,1]]<range/2.),knockoutVal,1.]*(Min[randGJmults[[i]],randGJmults[[#]]]*(a0/areas[[i]])*(adjMat[[i,#]]/(\[Pi]*\[CapitalDelta]x))(Subscript[ip3, #][t]-Subscript[ip3, i][t]))&/@(Flatten@adjMat[[i]]["NonzeroPositions"])],{i,numCells2D}];
getRandomCells[numCells_,var_]:=RandomVariate[LogNormalDistribution[-.5*Log[var+1],Sqrt[Log[var+1]]],numCells];

randExtParams={\[Alpha]};

randVals=Association[Table[randExtParams[[param]]->getRandomCells[numCells2D,.2],{param,Length@randExtParams}]];
randRules=Table[(#->#*randVals[[Key@#,cell]])&@randExtParams[[param]],{cell,numCells2D},{param,Length@randExtParams}];

gcampRandVals=RandomVariate[LogNormalDistribution[-\[Sigma]gcamp^2/2,\[Sigma]gcamp],numCells2D];
(*plcRandVals=RandomVariate[LogNormalDistribution[-.5*Log[varPLC+1],Sqrt[Log[varPLC+1]]],numCells2D];*)

maxGCAMPVal=Max@gcampRandVals; (*Needed to determine how to scale video output later*)

deqsTissueEQ=deqns/.{r\[Mu]T->0.,lr->0.};
(*Variable replacements for multiple cells*)
depVars2D=Table[Subscript[spec, i][t],{i,numCells2D},{spec,{ip3,h,c,cER}}];
variableReplacements2D=Flatten/@Table[{spec[t]->Subscript[spec, i][t],spec'[t]->Subscript[spec, i]'[t]},{i,numCells2D},{spec,{ip3,h,c,cER}}];
initReplacement2D=Flatten/@Table[{spec[0]->Subscript[spec, i][0]},{i,numCells2D},{spec,{ip3,h,c,cER}}]; 
multiCellReplacements2D=
Table[

Join[
variableReplacements2D[[i]],
initReplacement2D[[i]],
{Bx->Bx*If[gcampRand,gcampRandVals[[i]],1.]},
(*{\[Alpha]\[Rule]\[Alpha]*If[plcRand,plcRandVals[[i]],1.]},*)
(*{\[Alpha]0\[Rule]\[Alpha]0*If[plcRand,plcRandVals[[i]],1.]},*)
randRules[[i]],
{jGJIP3->fluxIP3Terms2DUW[[i]],jGJc->fluxGJCaTerms2DUW[[i]]},
{r->distances[[i]]}
],

{i,numCells2D}];
tissueCellEquationsEQ=deqsTissueEQ/.multiCellReplacements2D;

tissueModelEquationsEQ=Flatten[Table[tissueCellEquationsEQ[[i]]/.If[positions[[i,1]]<range/2,knockoutRules,{}]/.If[positions[[i,1]]<range/2,overExpressRules,{}],{i,numCells2D}]]/.intParams;
tLong=300.; (*seconds*)
framesBeforeWounding=20;
tissueModelSolEQ=
ParametricNDSolveValue[
tissueModelEquationsEQ,
depVars2D/.Table[{t->\[Tau]},{\[Tau],tLong-(framesBeforeWounding-1)*\[CapitalDelta]t,tLong,\[CapitalDelta]t}],
{t,0,tLong},
extParams,
DependentVariables->Flatten@depVars2D,
"Method"->{"EquationSimplification"->{Automatic,"TimeConstraint"->100}}];

newEQSoln=tissueModelSolEQ@@(finalParamAssociation/@extParams);
newEQs=Table[{ip30->#[[i,1]],ca0->#[[i,3]],cER0->#[[i,4]],h0->#[[i,2]]},{i,numCells2D}]&@(Last@newEQSoln);
deqsTissue=deqns/.{lr->lrSol[r,t],ip3[0]==x_->ip3[0]==ip30,c[0]==x_->c[0]==ca0,cER[0]==x_->cER[0]==cER0,h[0]==x_->h[0]==h0};

(*Making ablated cells have constant calcium. Use Dispatch to make things faster. Might need to modify things other than Calcium*)
ablatedChanges=Dispatch@Flatten@Table[If[MemberQ[ablatedPos,i],{Subscript[c, i]'[t]==x_->Subscript[c, i]'[t]==0,Subscript[c, i][0]==x_->Subscript[c, i][0]==cExt},{Subscript[c, i]'[t]==x_->Subscript[c, i]'[t]==x,Subscript[c, i][0]==x_->Subscript[c, i][0]==x}],{i,numCells2D}];
(*Comapred to the equilibrium tissue, we add in the new resting levels from the previous section, and we take out the initial condition scaling (since we already have the initial conditions for each cell)*)
multiCellReplacements2D=
Table[

Join[
newEQs[[i]],
variableReplacements2D[[i]],
initReplacement2D[[i]],
{r\[Mu]T->r\[Mu]T*\[Mu]TDist[distances[[i]]]},
{Bx->Bx*If[gcampRand,gcampRandVals[[i]],1.]},
(*{\[Alpha]\[Rule]\[Alpha]*If[plcRand,plcRandVals[[i]],1.]},*)
(*{\[Alpha]0\[Rule]\[Alpha]0*If[plcRand,plcRandVals[[i]],1.]},*)
randRules[[i]],
{jGJIP3->fluxIP3Terms2D[[i]],jGJc->fluxGJCaTerms2D[[i]]},
{r->distances[[i]]}
],

{i,numCells2D}];
(*Break into separate steps to speed up ablation replacements*)
tissueCellEquations=deqsTissue/.multiCellReplacements2D;

tissueModelEquations=Flatten[Table[tissueCellEquations[[i]]/.If[MemberQ[ablatedPos,i],{Subscript[c, i]'[t]==x_->Subscript[c, i]'[t]==0,Subscript[c, i][0]==x_->Subscript[c, i][0]==cExt},{}]/.If[positions[[i,1]]<range/2,knockoutRules,{}]/.If[positions[[i,1]]<range/2,overExpressRules,{}],{i,numCells2D}]]/.intParams;
tissueModelSol=
ParametricNDSolveValue[
tissueModelEquations,
depVars2D[[;;,3]],
{t,0,400.},
extParams,
DependentVariables->Flatten@depVars2D,
"Method"->{"EquationSimplification"->{Automatic,"TimeConstraint"->100}}];

test=tissueModelSol@@(finalParamAssociation/@extParams);
testG=CtoF/@test;
(*Set frames before wounding to be resting calcium levels*)
animate2D[soln_,stamps_?BooleanQ]:=
Block[{timePoints,minF,maxF,styles,animList,frames},
timePoints=Join[CtoF@(newEQSoln[[;;,;;,3]]),Table[soln[[i]]/.t->\[Tau],{\[Tau],0.1,400.1,\[CapitalDelta]t},{i,numCells2D}]];

minF=.12;
maxF=1.;

styles=
Table[Style[{2,i},If[MemberQ[ablatedPos,i]&&\[Tau]>framesBeforeWounding,Black,GrayLevel[Min[maxF,If[gcampRand,gcampRandVals[[i]],1.]*Rescale[timePoints[[\[Tau],i]],{minInt,maxInt*maxGCAMPVal},{minF,maxF}]]]]],{\[Tau],Length@timePoints},{i,numCells2D}];

animList=HighlightMesh[mesh,#,MeshCellStyle->{1->Directive[Opacity[0],Antialiasing->False]}]&/@styles;

frames=If[stamps,MapIndexed[Show[#1,Graphics[{White,Text[Style["t = "<>ToString[Round[\[CapitalDelta]t*(First[#2]-(framesBeforeWounding+1)),.1]]<>" s",White,18],{75,35}],White,Thickness@.025,Line[{{.75*range,.1*range},{.75*range+50,.1*range}}]}]]&,animList],

animList];

Return[frames]

(*ListAnimate[frames,7]*)
]
frames2D=animate2D[testG,False];

Export["FullModel_param6_"<>knockoutSpec<>"_KDN_GCaMP_PLC_Varied.tif",frames2D,"FrameRate"->10];
