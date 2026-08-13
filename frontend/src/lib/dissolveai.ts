/**
 * Facts in this file are taken from the DissolveAI repository
 * (github.com/Iram768/dissolveai) — nothing here is invented.
 *
 * Sources:
 *  - solvent_project/model_training_BigSolDB/app.py            (Flask API contract)
 *  - solvent_project/model_training_BigSolDB/models/metadata.json
 *  - solvent_project/data_processed/BigSolDB_clean.csv          (solvent vocabulary)
 *  - solvent_project/model_training_BigSolDB/train_full_pipeline_BigSol.py
 *  - README.md
 */

/** Exact 70 solvent categories present in BigSolDB_clean.csv / the one-hot encoder. */
export const SUPPORTED_SOLVENTS = [
  "water",
  "ethanol",
  "methanol",
  "isopropanol",
  "n-propanol",
  "n-butanol",
  "isobutanol",
  "sec-butanol",
  "tert-butanol",
  "n-pentanol",
  "isopentanol",
  "n-hexanol",
  "n-heptanol",
  "n-octanol",
  "n-nonanol",
  "n-decanol",
  "ethylene glycol",
  "propylene glycol",
  "acetone",
  "2-butanone",
  "cyclohexanone",
  "cyclopentanone",
  "MIBK",
  "acetonitrile",
  "DMSO",
  "DMF",
  "DMAc",
  "NMP",
  "formamide",
  "THF",
  "1,4-dioxane",
  "diethyl ether",
  "MTBE",
  "anisole",
  "2-methoxyethanol",
  "2-ethoxyethanol",
  "2-propoxyethanol",
  "2-butoxyethanol",
  "transcutol",
  "ethyl acetate",
  "methyl acetate",
  "n-propyl acetate",
  "isopropyl acetate",
  "n-butyl acetate",
  "isobutyl acetate",
  "n-pentyl acetate",
  "ethyl formate",
  "dimethyl carbonate",
  "gamma-butyrolactone",
  "acetic acid",
  "formic acid",
  "propionic acid",
  "n-butyric acid",
  "benzene",
  "toluene",
  "ethylbenzene",
  "chlorobenzene",
  "o-xylene",
  "m-xylene",
  "p-xylene",
  "chloroform",
  "dichloromethane",
  "1,2-dichloroethane",
  "tetrachloromethane",
  "cyclohexane",
  "n-pentane",
  "n-hexane",
  "n-heptane",
  "n-octane",
  "n-dodecane",
] as const;

export type Solvent = (typeof SUPPORTED_SOLVENTS)[number];

/** Groupings are presentational only; every name above is a real encoder category. */
export const SOLVENT_GROUPS: { label: string; solvents: string[] }[] = [
  {
    label: "Aqueous & alcohols",
    solvents: [
      "water",
      "ethanol",
      "methanol",
      "isopropanol",
      "n-propanol",
      "n-butanol",
      "isobutanol",
      "sec-butanol",
      "tert-butanol",
      "n-pentanol",
      "isopentanol",
      "n-hexanol",
      "n-heptanol",
      "n-octanol",
      "n-nonanol",
      "n-decanol",
      "ethylene glycol",
      "propylene glycol",
    ],
  },
  {
    label: "Polar aprotic",
    solvents: [
      "acetone",
      "2-butanone",
      "cyclohexanone",
      "cyclopentanone",
      "MIBK",
      "acetonitrile",
      "DMSO",
      "DMF",
      "DMAc",
      "NMP",
      "formamide",
      "gamma-butyrolactone",
    ],
  },
  {
    label: "Ethers & glycol ethers",
    solvents: [
      "THF",
      "1,4-dioxane",
      "diethyl ether",
      "MTBE",
      "anisole",
      "2-methoxyethanol",
      "2-ethoxyethanol",
      "2-propoxyethanol",
      "2-butoxyethanol",
      "transcutol",
    ],
  },
  {
    label: "Esters & carbonates",
    solvents: [
      "ethyl acetate",
      "methyl acetate",
      "n-propyl acetate",
      "isopropyl acetate",
      "n-butyl acetate",
      "isobutyl acetate",
      "n-pentyl acetate",
      "ethyl formate",
      "dimethyl carbonate",
    ],
  },
  {
    label: "Acids",
    solvents: [
      "acetic acid",
      "formic acid",
      "propionic acid",
      "n-butyric acid",
    ],
  },
  {
    label: "Aromatics",
    solvents: [
      "benzene",
      "toluene",
      "ethylbenzene",
      "chlorobenzene",
      "o-xylene",
      "m-xylene",
      "p-xylene",
    ],
  },
  {
    label: "Halogenated",
    solvents: [
      "chloroform",
      "dichloromethane",
      "1,2-dichloroethane",
      "tetrachloromethane",
    ],
  },
  {
    label: "Alkanes",
    solvents: [
      "cyclohexane",
      "n-pentane",
      "n-hexane",
      "n-heptane",
      "n-octane",
      "n-dodecane",
    ],
  },
];

/** From models/metadata.json in the repository. */
export const MODEL_METADATA = {
  temperatureMinK: 243.15,
  temperatureMaxK: 425.77,
  solubilityLowerBound: -12.0,
  solubilityUpperBound: 2.0,
  solubilityThreshold: -2.0,
  totalSolvents: 70,
  defaultTemperatureK: 298.0,
} as const;

/** Rows in data_processed/BigSolDB_clean.csv (measured solubility values, mol/L). */
export const DATASET_ROWS = 100983;

export const EXAMPLE_MOLECULES: { smiles: string; label: string }[] = [
  { smiles: "CCO", label: "ethanol" },
  { smiles: "c1ccccc1", label: "benzene" },
  { smiles: "O", label: "water" },
  { smiles: "CC(=O)Oc1ccccc1C(=O)O", label: "aspirin" },
];

export const REPO_URL = "https://github.com/Iram768/dissolveai";
export const PROFILE_URL = "https://github.com/Iram768";
