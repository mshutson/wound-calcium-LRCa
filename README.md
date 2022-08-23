# wound-calcium-LRCa

The files in this repository are Mathematica notebooks (.nb) and package files (.m). These notebooks were written to model wound-induced calcuim responses in epithelia that follow from both physical damage to the cells' plasma membranes and protease cleavage of a pro-ligand to release active ligand that then binds to a G-protein coupled receptor. The models are described in full in the publication: </BR></BR>
  "A mathematical model of calcium signals around laser-induced epithelial wounds"</BR>
  Aaron C. Stevens, James T. O’Connor, Andrew D. Pumford, Andrea Page-McCaw, M. Shane Hutson</BR>  
<BR/>

<B>Mathematica notebooks</B>
<UL>
<LI>FullCalciumSignalingModel_ParameterBuildUp_11.m: Notebook to test parameter values within various sub-models as well as set up implementation onto the ACCRE cluster

<LI>ACCRE_FullCalciumSignalingModel_Control.m: Notebook that results in a control model response. In this and subsequent notebooks, the list “randExtParams” determines which parameters will be varied across the entire tissue, but GCaMP variation is handled separately with the “gcampRand” flag.

<LI>ACCRE_FullCalciumSignalingModel_Knockdowns.m: Notebook specifically for the GJ and PLC knockdowns, although knockdown of any parameter on half of the tissue can be obtained by specifying parameters in the list “knockoutSpec” and adding a corresponding line to the “knockoutAssociations” association. Notebook was set up to run each knockdown simulation in parallel jobs on ACCRE.

<LI>ACCRE_FullCalciumSignalingModel_JTG.m: Notebook for “Jump the Gap” model output. Notebook was set up to test various wound to pnr border distances in parallel jobs on ACCRE.

<LI>ACCRE_FullCalciumSignalingModel_freeCalcium.m: Notebook that outputs a free calcium heat map video along with its corresponding GCaMP video.

<LI>ACCRE_FullCalciumSignalingModel_GCaMPKnockdown.m: Notebook specifically for the GCaMP knockdowns. It also outputs calcium vs time for pre-selected cells (cell numbers in file cellPairs.m)

<LI>ACCRE_FullCalciumSignalingModel_GJ Fluxes.m: Notebook that outputs gap junction fluxes across each cell edge in a control wound simulation. Output files “randGJmults.m” and “gcampRandVals.m” are needed for flux analysis.

<LI>ACCRE_FullCalciumSignalingModel_SelectiveGJ.m: Notebook for selective gap junctions (knockdown of calcium or IP3 transfer). Notebook is set up to run each knockdown simulation in parallel jobs on ACCRE

<LI>ACCRE_FullCalciumSignalingModel_ParameterVariation.m: Notebook for runs where a single parameter is varied across the tissue. Notebook is set up to run each simulation in parallel jobs on ACCRE

<LI>fig*.m: Notebook for analysis and figure creation of Figure * in the main text.
</UL>


<B>Other Files</B>
<UL>
<LI>microtearDistribution_Andrew_dataPoints.xlsx: Micro-tear distribution data used to determine the microtear distribution function

<LI>controlCaRadData.m: Calcium signal radius vs time experimental data of a control wound. Used to initialize some values in the notebooks

<LI>LRParams_lowthresh.m: Parameters for the Ligand-Receptor reaction-diffusion model. These result from fits to a distal calcium response with a signaling threshold of 5% receptor activation. The first parameter set is used in notebooks above.

<LI>mesh_4826cells_450x450.m: Voronoi mesh of the tissue

<LI>gcampRandVals.m: Random GCaMP multipliers (needed for GJ flux analysis)

<LI>randGJmults.m: Random GJ multipliers (needed for GJ flux analysis)

<LI>cellPairs.m: Pairs of cell labels for cells along the horizontal mid-line of the tissue mesh that are roughly the same distance from the wound. Used for free cytosolic calcium concentration analysis. 
  </UL>
