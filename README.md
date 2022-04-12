# WSMEmodel
WSME Model Codes - SSA + DSA + DSAw/L Approximations

Cite:  Thermodynamics and folding landscapes of large proteins from a statistical mechanical model. Gopi S, Aranganathan A, Naganathan AN. Curr Res Struct Biol. 2019 Oct 23;1:6-12. doi: 10.1016/j.crstbi.2019.10.002

Step 1: Identify a PDB file and ensure that there are no missing atoms/residues or breaks (preferably between 50 - 700 residues)

Step 2: Generate the STRIDE output as per the format attached and save as 'struct.txt' (this needs to be generated for every protein of interest)

Step 3: Modify the input parameters in cmapCalcElecBlock.m and exceute it. Multiple outputs will be generated. This will take only a few seconds.

Step 4: Execute FesCalc_Block.m to generate free-energy profiles and surfaces of the protein under consideration. Longer the protein, more will be the time taken. For CI2, a 65 residue protein, it takes less than 3 minutes to generate the conformational landscape.

Caution 1 - Large block sizes (>6) have not been tested. Larger the block size, more will be the smoothening of the free-energy profile/surface - but the outputs can be generated faster.

Caution 2 - For a more quantiative comparison with experiments (for example, stability), the parameters need to be iteratively modified.

Caution 3 - The current approximation underestimates the configurational entropy of the unfolded state and hence should be used with care.

Caution 4 - The version of the code here can be directly used to simulate thermodynamic unfolding behaviors of water-soluble proteins with a dielectric constant of 29 (DOI: 10.1021/ct300676w). For membrane proteins, the dielectric constant needs to be changed to 4 (instead of 29) in the files FesCalc_Block.m (Line 87) and in cmapCalcElecBlock.m (Line 123).
