# WSMEmodel
WSME Model Codes - SSA + DSA + DSAw/L Approximations

Step 1: Identify a PDB file and ensure there are no missing atoms/residues or breaks (preferably between 50 - 700 residues)

Step 2: Generate the STRIDE output as per the format attached and save as 'struct.txt'

Step 3: Modify the input parameters in cmapCalcElecBlock.m and exceute it. Multiple outputs will be generated.

Step 4: Execute FesCalc_Block.m to generate free-energy profiles and surfaces of the protein under consideration. 

Caution 1 - Large block sizes have not been tested. Larger the block size, more will be the smoothening of the free-energy profile/surface - but the outputs can be generated faster.

Cauton 2 - For a more quantiative comparison with experiments (for example, stability), the parameters need to be iterative modified.
