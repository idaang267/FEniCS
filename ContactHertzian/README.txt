Contact: Ida Ang (ia267@cornell.edu)
Affiliation:
  * Former graduate student from the Bouklas Laboratory at Cornell University
  * For more details on previous research an alternate contact is
    Nikolaos Bouklas (nb589@cornell.edu)

Details:
  * Related to the Numerical tours of continuum mechanics using FEniCS series
    "Hertzian contact with a rigid indenter using a penalty approach"
    (see citation below)
      * This was a fully 3 dimensional formulation
  * Paraview recommended for visualization of results

Specific FEniCS Details and Changes Made:
  * These code files demonstrate the penalty approach for contact under a linear
    elastic assumption for the constitutive relationships
      * Simplified to 2D in 2DContactHertzPenalty.py
      * Original with adjustments for visualization and post-processing in
        3DContactHertzPenalty.py
      * PostGapCases.py is for Postprocessing for visualization
  * SNES code specifically shows projection of boundaries to the mesh for
    visualization of boundaries where constraints can be applied
  * Basic for loop to solve problem at multiple iterations
  * Multiple fields projected to the same mesh for visualization
  * Demonstrates how to save multiple fields at multiple time steps to .txt
    files for future post-processing with PostGapCases.py

@manual{bleyer2018numericaltours,
title={Numerical Tours of Computational Mechanics with {FE}ni{CS}},
DOI={10.5281/zenodo.1287832},
howpublished = {https://comet-fenics.readthedocs.io},
publisher={Zenodo},
author={Jeremy Bleyer},
year={2018}}
