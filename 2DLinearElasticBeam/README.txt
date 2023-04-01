Contact: Ida Ang (ia267@cornell.edu)
Affiliation:
  * Former graduate student from the Bouklas Laboratory at Cornell University
  * For more details on previous research an alternate contact is
    Nikolaos Bouklas (nb589@cornell.edu)

Details:
  * Related to the Numerical tours of continuum mechanics using FEniCS:
    2D Linear Elasticity
    (see citation below)
  * Paraview recommended for visualization of results

Specific FEniCS Details and Changes Made:
  * These code files demonstrates a 2D linear elastic formulation
      * plane strain and plane stress assumptions
      * PostProc.py is for post-processing for visualization
  * Adds functionality for user parameters to be parsed from command line
      * python3 Demo2DLE.py --VariableName VariableValue
      * python3 Demo2DLE.py --model plane_strain
  * Basic for loop to solve problem at multiple iterations
  * Demonstrates how to save multiple fields at multiple time steps to .txt
    files for future post-processing with PostProc.py

@manual{bleyer2018numericaltours,
title={Numerical Tours of Computational Mechanics with {FE}ni{CS}},
DOI={10.5281/zenodo.1287832},
howpublished = {https://comet-fenics.readthedocs.io},
publisher={Zenodo},
author={Jeremy Bleyer},
year={2018}}
