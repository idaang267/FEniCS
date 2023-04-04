Contact: Ida Ang (ia267@cornell.edu)
Affiliation:
  * Former graduate student from the Bouklas Laboratory at Cornell University
  * For more details on previous research an alternate contact is
    Nikolaos Bouklas (nb589@cornell.edu)

Details:
  * Related to the Programmer’s reference for DOLFIN (Python)
      * Poisson's Equation
      * Poisson's Equation with periodic boundary conditions
      * Mixed Formulation for Poisson's Equation
      * NOTE: I am not taking credit for this work, all original demos can be
        found in the links provided in the pdf document and the original code
        with all credits and information can still be accessed and downloaded.
  * Paraview recommended for visualization of results

Specific FEniCS Details and Changes Made:
  * The code files in order are
      * DemoPoisson.py
      * DemoPeriodic.py
      * DemoMixedPoisson.py
  * Projection of boundaries to the mesh for visualization of boundaries
    where boundary conditions are applied
  * Multiple fields projected to the same mesh for visualization
  * Demonstrates some postprocessing in DemoPeriodic.py
  * NOTE: for a list of important changes to the code, see the pdf document.
    DemoPeriodic.py was working in the past but with a warning message, which
    has now been fixed. DemoMixedPoisson.py had the most changes to
    implementation. 
