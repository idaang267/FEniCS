Read me for .py codes

Plane Strain Formulations:
  2D-planestrain-discrete.py
    RectangleMesh or 2DShearTestHalf.xml
    (Displacement + Pressure)
  2D-planestrain-stabilized.py
    RectangleMesh or 2DShearTestRef.xml
    (Displacement + Pressure) + (Damage) + Stabilization Formulation

Plane Stress Formulations:
  2D-planestress-discrete-shear-test.py
    RectangleMesh or 2DShearTestHalf.xml
    (Displacement + Pressure + F33)

  2D-planestress-discrete-stabilized.py
    RectangleMesh or 2DShearTestHalf.xml (Requires a regular mesh)
    (Displacement + Pressure + F33)
    Stated in a weak form sense? not sure what was happening here

  2D-planestress-stabilized.py
    RectangleMesh (Requires a regular mesh)
    (Displacement + Pressure + F33) + (Damage) + Stabilization Formulation

  2D-planestress-TH.py
    RectangleMesh or 2DShearTestRef.xml
    (Displacement + Pressure + F33) + (Damage) Formulation

  2D-planestress-TH-asymptotic.py
    RectangleMesh or 2DSquare3Ref.xml
    (Displacement + Pressure + F33) + (Damage) Formulation
    Above with asymptotic boundary conditions applied

  2D-planestress-u.py
    RectangleMesh or 2DShearTestRef.xml
    (Displacement + F33) + (Damage)

Original 3D Formulations
  3D-TH-Con.py
    BoxMesh or ShearTestRefStruct.xml
    (Displacement + Pressure) + (Damage) Formulation
    Special conditional statement
  3D-TH-Orig.py
    BoxMesh or ShearTestRefStruct.xml
    (Displacement + Pressure) + (Damage) Formulation
  3D-traction-stabilized.py
    BoxMesh
    (Displacement + Pressure) + (Damage) + Stabilization Formulation
    Original 3D code with original formulation
  3D-decomposition-stabilized.py
    BoxMesh
    (Displacement + Pressure) + (Damage) + Stabilization Formulation
    Considers the decomposition of energy into active and passive terms
  3D-Taylor-Hood-Hybrid.py
    BoxMesh
    (Displacement + Pressure) + (Damage) Formulation
    Bin's version of the code above
      Not strictly accurate because of second term of elastic energy

Miscellaneous Code
  isotropic-traction.py
    My version
  isotropic-traction_j.py
    Jason's version
  isotropic-traction_orig.py
    Bin's original version
  weak-kinking.py and anisotropic-gradient-damage.py
    from Bin's papers
  convert.py
    converts xml for viewing in FEniCS
