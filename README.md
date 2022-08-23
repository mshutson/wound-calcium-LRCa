# wound-calcium-LRCa

<B>Mathematica notebooks</B>

•	FullCalciumSignalingModel_ParameterBuildUp_11.m: Notebook to test parameter values within various sub-models as well as set up implementation onto the ACCRE cluster

•	ACCRE_FullCalciumSignalingModel_Control.m: Notebook that results in a control model response. In this and subsequent notebooks, the list “randExtParams” determines which parameters will be varied across the entire tissue, but GCaMP variation is handled separately with the “gcampRand” flag.

•	ACCRE_FullCalciumSignalingModel_Knockdowns.m: Notebook specifically for the GJ and PLC knockdowns, although knockdown of any parameter on half of the tissue can be obtained by specifying parameters in the list “knockoutSpec” and adding a corresponding line to the “knockoutAssociations” association. Notebook was set up to run each knockdown simulation in parallel jobs on ACCRE.

•	ACCRE_FullCalciumSignalingModel_JTG.m: Notebook for “Jump the Gap” model output. Notebook was set up to test various wound to pnr border distances in parallel jobs on ACCRE.

•	ACCRE_FullCalciumSignalingModel_freeCalcium.m: Notebook that outputs a free calcium heat map video along with its corresponding GCaMP video.

•	ACCRE_FullCalciumSignalingModel_GCaMPKnockdown.m: Notebook specifically for the GCaMP knockdowns. It also outputs calcium vs time for pre-selected cells (cell numbers in file cellPairs.m)

•	ACCRE_FullCalciumSignalingModel_GJ Fluxes.m: Notebook that outputs gap junction fluxes across each cell edge in a control wound simulation. Output files “randGJmults.m” and “gcampRandVals.m” are needed for flux analysis.

•	ACCRE_FullCalciumSignalingModel_SelectiveGJ.m: Notebook for selective gap junctions (knockdown of calcium or IP3 transfer). Notebook is set up to run each knockdown simulation in parallel jobs on ACCRE

•	ACCRE_FullCalciumSignalingModel_ParameterVariation.m: Notebook for runs where a single parameter is varied across the tissue. Notebook is set up to run each simulation in parallel jobs on ACCRE

•	fig*.m: Notebook for analysis and figure creation of Figure * in the main text.
 

<B>Other Files</B>

•	microtearDistribution_Andrew_dataPoints.xlsx: Micro-tear distribution data used to determine the microtear distribution function

•	controlCaRadData.m: Calcium signal radius vs time experimental data of a control wound. Used to initialize some values in the notebooks

•	LRParams_lowthresh.m: Parameters for the Ligand-Receptor reaction-diffusion model. These result from fits to a distal calcium response with a signaling threshold of 5% receptor activation. The first parameter set is used in notebooks above.

•	mesh_4826cells_450x450.m: Voronoi mesh of the tissue

•	gcampRandVals.m: Random GCaMP multipliers (needed for GJ flux analysis)

•	randGJmults.m: Random GJ multipliers (needed for GJ flux analysis)

•	cellPairs.m: Pairs of cell labels for cells along the horizontal mid-line of the tissue mesh that are roughly the same distance from the wound. Used for free cytosolic calcium concentration analysis. 
