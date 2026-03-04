**DissolveAI**

**Machine Learning-Based Solubility Prediction System**

DissolveAI is a supervised machine learning system developed to predict the solubility of chemical compounds in various solvents at room temperature. The goal is to assist early-stage research by providing fast computational estimates of solubility using cheminformatics and regression modeling.

**Problem Formulation**

Solubility prediction is treated as a regression problem.

**Inputs**

  Molecular structure (SMILES)
  
  Solvent identity

**Output**

  Continuous solubility value

**Feature Engineering**

Molecules are converted to numerical features using RDKit:

**Physicochemical Descriptors (8 features)**

  Molecular Weight
  
  LogP
  
  TPSA
  
  H-bond donors/acceptors
  
  Rotatable bonds
  
  Aromatic rings
  
  Fraction sp³

**Morgan Fingerprints**  

Radius = 2

1024-bit vector

**Solvent Encoding**

  One-hot encoding (scikit-learn)

Final feature matrix combines descriptors + fingerprints + solvent encoding.

**Model**

  Random Forest Regressor
  
  500 trees
  
  80/20 train-test split
  
  Evaluated using RMSE and 5-fold cross-validation

The model captures nonlinear relationships between structure, solvent, and solubility.

**System Architecture**

  Offline training pipeline
  
  Flask-based prediction backend
  
  Simple web interface prototype
  
  Suggests alternative solvents if predicted solubility is low

**Limitations**

  Performance depends on dataset coverage (BigSolDB)
  
  Limited extrapolation beyond training chemical space
  
  Purely data-driven (no explicit thermodynamic modeling)

**Technologies**

Python | RDKit | Scikit-learn | Flask | NumPy | Next.js

**Future Work**

Temperature-dependent modeling

Integration with graph neural networks

Expansion to larger solubility datasets
