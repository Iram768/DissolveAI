import { r as __toESM } from "../_runtime.mjs";
import { n as AnimatePresence } from "../_libs/framer-motion+[...].mjs";
import {
  n as require_jsx_runtime,
  r as require_react,
} from "../_libs/react+tanstack__react-query.mjs";
import { h as Link } from "../_libs/@tanstack/react-router+[...].mjs";
import { t as motion } from "../_libs/motion.mjs";
//#region node_modules/.nitro/vite/services/ssr/assets/routes-BOwAXmoD.js
var import_react = /* @__PURE__ */ __toESM(require_react());
var import_jsx_runtime = require_jsx_runtime();
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
var SUPPORTED_SOLVENTS = [
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
];
/** Groupings are presentational only; every name above is a real encoder category. */
var SOLVENT_GROUPS = [
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
var MODEL_METADATA = {
  temperatureMinK: 243.15,
  temperatureMaxK: 425.77,
  solubilityLowerBound: -12,
  solubilityUpperBound: 2,
  solubilityThreshold: -2,
  totalSolvents: 70,
  defaultTemperatureK: 298,
};
/** Rows in data_processed/BigSolDB_clean.csv (measured solubility values, mol/L). */
var DATASET_ROWS = 100983;
var EXAMPLE_MOLECULES = [
  {
    smiles: "CCO",
    label: "ethanol",
  },
  {
    smiles: "c1ccccc1",
    label: "benzene",
  },
  {
    smiles: "O",
    label: "water",
  },
  {
    smiles: "CC(=O)Oc1ccccc1C(=O)O",
    label: "aspirin",
  },
];
var REPO_URL = "https://github.com/Iram768/dissolveai";
var PROFILE_URL = "https://github.com/Iram768";
var NAV$1 = [
  {
    to: "#home",
    label: "Home",
  },
  {
    to: "#predictor",
    label: "Predictor",
  },
  {
    to: "#how-it-works",
    label: "How It Works",
  },
  {
    to: "#model",
    label: "Model",
  },
  {
    to: "#about",
    label: "About",
  },
];
function SiteHeader() {
  const [scrolled, setScrolled] = (0, import_react.useState)(false);
  const [open, setOpen] = (0, import_react.useState)(false);
  (0, import_react.useEffect)(() => {
    const onScroll = () => setScrolled(window.scrollY > 12);
    onScroll();
    window.addEventListener("scroll", onScroll, { passive: true });
    return () => window.removeEventListener("scroll", onScroll);
  }, []);
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("header", {
    className: `sticky top-0 z-50 transition-colors duration-300 ${scrolled ? "border-b border-border bg-background/85 backdrop-blur-md" : "border-b border-transparent"}`,
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
        className:
          "mx-auto flex max-w-6xl items-center gap-6 px-5 py-4 sm:px-8",
        children: [
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(Link, {
            to: "/",
            hash: "home",
            className: "group flex items-center gap-3",
            "aria-label": "DissolveAI home",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(FlaskMark, {}),
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", {
                className: "leading-tight",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                    className:
                      "block font-display text-lg tracking-tight text-foreground",
                    children: "DissolveAI",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                    className: "label-caps hidden sm:block",
                    children: "Molecular Solubility Prediction",
                  }),
                ],
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("nav", {
            className: "ml-auto hidden items-center gap-7 md:flex",
            "aria-label": "Main",
            children: [
              NAV$1.map((item) =>
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                  Link,
                  {
                    to: "/",
                    hash: item.to.replace("#", ""),
                    className:
                      "text-sm text-muted-foreground transition-colors hover:text-foreground",
                    children: item.label,
                  },
                  item.to,
                ),
              ),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("a", {
                href: REPO_URL,
                target: "_blank",
                rel: "noreferrer noopener",
                className:
                  "border border-input px-3 py-1.5 font-mono text-xs uppercase tracking-[0.14em] text-foreground transition-colors hover:border-primary hover:text-primary",
                children: "GitHub",
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
            type: "button",
            onClick: () => setOpen((v) => !v),
            "aria-expanded": open,
            "aria-controls": "mobile-nav",
            className:
              "ml-auto flex h-9 w-9 flex-col items-center justify-center gap-1.5 border border-input md:hidden",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                className: "sr-only",
                children: open ? "Close menu" : "Open menu",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                className: "block h-px w-4 bg-foreground",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                className: "block h-px w-4 bg-foreground",
              }),
            ],
          }),
        ],
      }),
      open &&
        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("nav", {
          id: "mobile-nav",
          "aria-label": "Main, mobile",
          className: "border-t border-border bg-background md:hidden",
          children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("ul", {
            className: "mx-auto max-w-6xl px-5 py-2 sm:px-8",
            children: [
              NAV$1.map((item) =>
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                  "li",
                  {
                    className: "border-b border-border last:border-b-0",
                    children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                      Link,
                      {
                        to: "/",
                        hash: item.to.replace("#", ""),
                        onClick: () => setOpen(false),
                        className: "block py-3 text-sm text-foreground",
                        children: item.label,
                      },
                    ),
                  },
                  item.to,
                ),
              ),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("li", {
                className: "py-3",
                children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("a", {
                  href: "https://github.com/Iram768/dissolveai",
                  target: "_blank",
                  rel: "noreferrer noopener",
                  className:
                    "font-mono text-xs uppercase tracking-[0.14em] text-primary",
                  children: "GitHub repository",
                }),
              }),
            ],
          }),
        }),
    ],
  });
}
function FlaskMark() {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("svg", {
    width: "26",
    height: "26",
    viewBox: "0 0 24 24",
    "aria-hidden": "true",
    className: "text-primary",
    fill: "none",
    stroke: "currentColor",
    strokeWidth: "1.1",
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
        d: "M9 3h6M10 3v5.2L4.9 18.4A2 2 0 0 0 6.7 21.4h10.6a2 2 0 0 0 1.8-3L14 8.2V3",
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
        d: "M7.4 15.4h9.2",
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("circle", {
        cx: "10.6",
        cy: "17.9",
        r: "0.9",
        fill: "currentColor",
        stroke: "none",
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("circle", {
        cx: "13.9",
        cy: "18.6",
        r: "0.6",
        fill: "currentColor",
        stroke: "none",
      }),
    ],
  });
}
var NAV = [
  {
    to: "#home",
    label: "Home",
  },
  {
    to: "#predictor",
    label: "Predictor",
  },
  {
    to: "#how-it-works",
    label: "How It Works",
  },
  {
    to: "#model",
    label: "Model",
  },
  {
    to: "#about",
    label: "About",
  },
];
var CONTACT = [
  {
    label: "Email",
    href: "mailto:iramjaved751@gmail.com",
    text: "iramjaved751@gmail.com",
  },
  {
    label: "LinkedIn",
    href: "https://www.linkedin.com/in/iram-javed/",
    text: "linkedin.com/in/iram-javed",
  },
  {
    label: "Portfolio",
    href: "https://iram-javed-portfolio.vercel.app/",
    text: "iram-javed-portfolio.vercel.app",
  },
  {
    label: "GitHub",
    href: PROFILE_URL,
    text: "github.com/Iram768",
  },
  {
    label: "Repository",
    href: REPO_URL,
    text: "github.com/Iram768/dissolveai",
  },
];
function SiteFooter() {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("footer", {
    className: "mt-24 border-t border-border bg-secondary/60",
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
        className:
          "mx-auto grid max-w-6xl gap-10 px-5 py-14 sm:px-8 md:grid-cols-3",
        children: [
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                className: "font-display text-xl text-foreground",
                children: "DissolveAI",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                className: "label-caps mt-1",
                children: "Molecular Solubility Prediction",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                className:
                  "mt-4 max-w-xs text-sm leading-relaxed text-muted-foreground",
                children:
                  "A supervised machine-learning system for estimating compound solubility across solvents using RDKit-derived features and a Random Forest regressor.",
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("nav", {
            "aria-label": "Footer",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", {
                className: "label-caps",
                children: "Navigation",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("ul", {
                className: "mt-4 space-y-2",
                children: NAV.map((item) =>
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                    "li",
                    {
                      children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                        Link,
                        {
                          to: "/",
                          hash: item.to.replace("#", ""),
                          className:
                            "text-sm text-muted-foreground transition-colors hover:text-primary",
                          children: item.label,
                        },
                      ),
                    },
                    item.to,
                  ),
                ),
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", {
                className: "label-caps",
                children: "Contact",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("ul", {
                className: "mt-4 space-y-2",
                children: CONTACT.map((item) =>
                  /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                    "li",
                    {
                      className: "text-sm",
                      children: [
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                          className:
                            "mr-2 font-mono text-xs uppercase tracking-[0.12em] text-muted-foreground",
                          children: item.label,
                        }),
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("a", {
                          href: item.href,
                          target: item.href.startsWith("mailto:")
                            ? void 0
                            : "_blank",
                          rel: "noreferrer noopener",
                          className:
                            "break-all text-foreground underline decoration-border underline-offset-4 transition-colors hover:text-primary hover:decoration-primary",
                          children: item.text,
                        }),
                      ],
                    },
                    item.label,
                  ),
                ),
              }),
            ],
          }),
        ],
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
        className: "border-t border-border",
        children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
          className:
            "mx-auto max-w-6xl px-5 py-5 font-mono text-xs text-muted-foreground sm:px-8",
          children:
            "DissolveAI — research prototype. Predictions are computational estimates and are not a substitute for experimental measurement.",
        }),
      }),
    ],
  });
}
/**
 * Hero visual: a restrained scientific illustration — a slowly rotating
 * molecular skeleton inside measurement rules, with a few drifting solvent
 * particles. Pure SVG/CSS, no WebGL, ~20 animated nodes.
 */
function HeroLab() {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
    className: "relative",
    children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
      className: "rule-grid border border-border bg-card p-4 sm:p-6",
      children: [
        /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
          className:
            "flex items-center justify-between border-b border-border pb-3",
          children: [
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
              className: "label-caps",
              children: "Fig. 1 — solute in solvent field",
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
              className: "font-mono text-[0.65rem] text-muted-foreground",
              children: "298 K",
            }),
          ],
        }),
        /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("svg", {
          viewBox: "0 0 320 300",
          role: "img",
          "aria-label":
            "Illustration of a molecular skeleton surrounded by solvent particles and measurement guides",
          className: "mt-4 h-auto w-full",
          children: [
            /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("g", {
              stroke: "var(--rule-strong)",
              strokeWidth: "1",
              fill: "none",
              children: [
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                  d: "M24 34 H296",
                  strokeDasharray: "2 6",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                  d: "M24 266 H296",
                  strokeDasharray: "2 6",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                  d: "M24 34 V44 M296 34 V44",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                  d: "M24 256 V266 M296 256 V266",
                }),
              ],
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("text", {
              x: "160",
              y: "26",
              textAnchor: "middle",
              className: "fill-muted-foreground",
              style: {
                font: "9px var(--font-mono)",
                letterSpacing: "0.14em",
              },
              children: "SOLVENT SHELL",
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("text", {
              x: "160",
              y: "286",
              textAnchor: "middle",
              className: "fill-muted-foreground",
              style: {
                font: "9px var(--font-mono)",
                letterSpacing: "0.14em",
              },
              children: "log mol/L ESTIMATE",
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("g", {
              className: "animate-drift",
              style: { transformOrigin: "160px 150px" },
              children: [
                [52, 78],
                [268, 96],
                [70, 214],
                [258, 210],
                [160, 58],
                [104, 128],
                [220, 176],
                [186, 240],
              ].map(([cx, cy], i) =>
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                  "g",
                  {
                    children: [
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("circle", {
                        cx,
                        cy,
                        r: i % 3 === 0 ? 4 : 2.6,
                        fill: "var(--teal-soft)",
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("circle", {
                        cx,
                        cy,
                        r: i % 3 === 0 ? 4 : 2.6,
                        fill: "none",
                        stroke: "var(--teal)",
                        strokeWidth: "0.6",
                      }),
                    ],
                  },
                  `${cx}-${cy}`,
                ),
              ),
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("g", {
              style: {
                transformOrigin: "160px 150px",
                animation: "dissolve-spin 64s linear infinite",
              },
              children: [
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("g", {
                  stroke: "var(--ink-soft)",
                  strokeWidth: "1.4",
                  fill: "none",
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                      d: "M110 150 L136 106 L188 106 L214 150 L188 194 L136 194 Z",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                      d: "M140 112 L184 112 M140 188 L184 188",
                      strokeWidth: "1",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                      d: "M214 150 L252 150",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                      d: "M252 150 L272 118",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                      d: "M110 150 L78 150",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("path", {
                      d: "M136 106 L120 70",
                    }),
                  ],
                }),
                [
                  [110, 150],
                  [136, 106],
                  [188, 106],
                  [214, 150],
                  [188, 194],
                  [136, 194],
                  [252, 150],
                ].map(([cx, cy]) =>
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                    "circle",
                    {
                      cx,
                      cy,
                      r: "3.4",
                      fill: "var(--ink)",
                    },
                    `${cx}-${cy}`,
                  ),
                ),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("circle", {
                  cx: "272",
                  cy: "118",
                  r: "5",
                  fill: "var(--paper)",
                  stroke: "var(--teal)",
                  strokeWidth: "1.4",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("text", {
                  x: "280",
                  y: "112",
                  className: "fill-primary",
                  style: { font: "9px var(--font-mono)" },
                  children: "OH",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("circle", {
                  cx: "78",
                  cy: "150",
                  r: "5",
                  fill: "var(--paper)",
                  stroke: "var(--ochre)",
                  strokeWidth: "1.4",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("circle", {
                  cx: "120",
                  cy: "70",
                  r: "4",
                  fill: "var(--paper)",
                  stroke: "var(--ink-soft)",
                  strokeWidth: "1.2",
                }),
              ],
            }),
          ],
        }),
        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
          className: "mt-4 grid grid-cols-3 gap-3 border-t border-border pt-3",
          children: [
            ["RDKit", "descriptors"],
            ["Morgan", "radius 2 · 1024 bit"],
            ["Random Forest", "500 trees"],
          ].map(([k, v]) =>
            /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
              "div",
              {
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className: "font-mono text-[0.68rem] text-foreground",
                    children: k,
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className: "font-mono text-[0.6rem] text-muted-foreground",
                    children: v,
                  }),
                ],
              },
              k,
            ),
          ),
        }),
      ],
    }),
  });
}
/**
 * 2D structure rendered client-side from the SMILES the user typed, using
 * smiles-drawer. Nothing decorative: if the string cannot be parsed we say so.
 */
function MoleculeStructure({ smiles, height = 260, className }) {
  const svgRef = (0, import_react.useRef)(null);
  const [parseError, setParseError] = (0, import_react.useState)(false);
  (0, import_react.useEffect)(() => {
    let cancelled = false;
    const svg = svgRef.current;
    if (!svg) return;
    if (!smiles.trim()) {
      svg.innerHTML = "";
      setParseError(false);
      return;
    }
    (async () => {
      const mod = await import("../_libs/smiles-drawer.mjs").then((n) => n.t);
      if (cancelled || !svgRef.current) return;
      const lib = mod.default ?? mod;
      const drawer = new lib.SvgDrawer({
        width: 480,
        height,
        bondThickness: 1.1,
        bondLength: 18,
        atomVisualization: "default",
        fontSizeLarge: 6,
        fontSizeSmall: 4,
        padding: 18,
        themes: {
          light: {
            C: "#3a3f47",
            O: "#0f6b6b",
            N: "#1f4f7a",
            F: "#2f7a4f",
            CL: "#2f7a4f",
            BR: "#8a5a1f",
            I: "#6b3f7a",
            P: "#a06a1f",
            S: "#a08a1f",
            B: "#7a5a3f",
            SI: "#7a7a7a",
            H: "#5b6068",
            BACKGROUND: "transparent",
          },
        },
      });
      lib.parse(
        smiles.trim(),
        (tree) => {
          if (cancelled || !svgRef.current) return;
          try {
            svgRef.current.innerHTML = "";
            drawer.draw(tree, svgRef.current, "light");
            setParseError(false);
          } catch {
            setParseError(true);
          }
        },
        () => {
          if (cancelled) return;
          if (svgRef.current) svgRef.current.innerHTML = "";
          setParseError(true);
        },
      );
    })();
    return () => {
      cancelled = true;
    };
  }, [smiles, height]);
  return /* @__PURE__ */ (0, import_jsx_runtime.jsx)("figure", {
    className,
    children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
      className:
        "relative overflow-hidden border border-border bg-card rule-grid",
      children: [
        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("svg", {
          ref: svgRef,
          role: "img",
          "aria-label": smiles.trim()
            ? `Two-dimensional structure drawn from SMILES ${smiles}`
            : "No molecular structure entered yet",
          viewBox: "0 0 480 260",
          className: "h-auto w-full",
          style: { maxHeight: height },
        }),
        !smiles.trim() &&
          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("figcaption", {
            className:
              "absolute inset-0 flex items-center justify-center px-6 text-center text-sm text-muted-foreground",
            children: "The structure appears here as you type a SMILES string.",
          }),
        parseError &&
          smiles.trim() &&
          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("figcaption", {
            className:
              "absolute inset-0 flex items-center justify-center px-6 text-center text-sm text-destructive",
            children:
              "This SMILES string could not be parsed into a structure.",
          }),
      ],
    }),
  });
}
/**
 * API client for the existing DissolveAI Flask backend
 * (solvent_project/model_training_BigSolDB/app.py).
 *
 * Contract, verbatim from that file:
 *   POST {base}/predict
 *   body: { solute_smiles: string, solvent_name: string, temperature?: number }
 *   200: { solute_name, solvent_name, temperature_used_K, predicted_solubility,
 *          message?, note?, suggested_solvents?: [{solvent, predicted_solubility}] }
 *   200 (unknown solvent): { solute_name, message }
 *   400: { error }        500: { error, detail }
 *
 * Units: predicted_solubility is reported by the backend in log mol/L and is
 * clipped to [solubility_lower_bound, solubility_upper_bound] from metadata.json.
 */
var DEFAULT_BASE_URL = "http://127.0.0.1:5000";
var STORAGE_KEY = "dissolveai.backend-url";
function getBackendUrl() {
  if (typeof window !== "undefined") {
    const stored = window.localStorage.getItem(STORAGE_KEY);
    if (stored && stored.trim()) return stored.trim().replace(/\/+$/, "");
  }
  return (
    {
      BASE_URL: "/",
      DEV: false,
      MODE: "production",
      PROD: true,
      SSR: true,
      TSS_DEV_SERVER: "false",
      TSS_DEV_SSR_STYLES_BASEPATH: "/",
      TSS_DEV_SSR_STYLES_ENABLED: "true",
      TSS_DISABLE_CSRF_MIDDLEWARE_WARNING: "false",
      TSS_INLINE_CSS_ENABLED: "false",
      TSS_ROUTER_BASEPATH: "",
      TSS_SERVER_FN_BASE: "/_serverFn/",
    }["VITE_DISSOLVEAI_API_URL"]?.trim() || DEFAULT_BASE_URL
  ).replace(/\/+$/, "");
}
function setBackendUrl(url) {
  if (typeof window === "undefined") return;
  if (url.trim()) window.localStorage.setItem(STORAGE_KEY, url.trim());
  else window.localStorage.removeItem(STORAGE_KEY);
}
var DEFAULT_BACKEND_URL = DEFAULT_BASE_URL;
var PredictionError = class extends Error {
  title;
  detail;
  constructor(title, message, detail) {
    super(message);
    this.name = "PredictionError";
    this.title = title;
    this.detail = detail;
  }
};
async function predictSolubility(input) {
  const url = `${getBackendUrl()}/predict`;
  let response;
  try {
    response = await fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        solute_smiles: input.smiles,
        solvent_name: input.solvent,
        temperature: input.temperatureK ?? null,
      }),
    });
  } catch {
    throw new PredictionError(
      "Prediction service unreachable",
      `No response from the DissolveAI backend at ${getBackendUrl()}. Start the Flask service (app.py) and confirm the endpoint below.`,
    );
  }
  let payload;
  try {
    payload = await response.json();
  } catch {
    throw new PredictionError(
      "Unreadable response",
      "The backend replied with a payload that is not valid JSON.",
    );
  }
  if (!response.ok || payload.error) {
    const error =
      payload.error ?? `Request failed with status ${response.status}`;
    if (/invalid smiles/i.test(error))
      throw new PredictionError(
        "Unable to analyze molecule",
        "The provided SMILES could not be parsed by RDKit. Please check the molecular structure and try again.",
      );
    if (/no solute smiles/i.test(error))
      throw new PredictionError(
        "No molecular structure supplied",
        "Enter a SMILES string before running a prediction.",
      );
    throw new PredictionError("Prediction failed", error, payload.detail);
  }
  if (payload.predicted_solubility === void 0)
    throw new PredictionError(
      "Solvent not covered by the model",
      payload.message ??
        "The selected solvent is not present in the trained solvent encoder.",
    );
  return payload;
}
function Predictor() {
  const [smiles, setSmiles] = (0, import_react.useState)("");
  const [solvent, setSolvent] = (0, import_react.useState)("water");
  const [temperature, setTemperature] = (0, import_react.useState)("");
  const [status, setStatus] = (0, import_react.useState)("idle");
  const [result, setResult] = (0, import_react.useState)(null);
  const [error, setError] = (0, import_react.useState)(null);
  const [smilesError, setSmilesError] = (0, import_react.useState)(null);
  const [showEndpoint, setShowEndpoint] = (0, import_react.useState)(false);
  const [endpoint, setEndpoint] = (0, import_react.useState)(() =>
    typeof window === "undefined" ? DEFAULT_BACKEND_URL : getBackendUrl(),
  );
  const tempNumber = temperature.trim() === "" ? null : Number(temperature);
  const tempOutOfRange =
    tempNumber !== null &&
    (Number.isNaN(tempNumber) ||
      tempNumber < MODEL_METADATA.temperatureMinK ||
      tempNumber > MODEL_METADATA.temperatureMaxK);
  async function run(event) {
    event.preventDefault();
    setError(null);
    if (!smiles.trim()) {
      setSmilesError("Enter a SMILES string to describe the molecule.");
      return;
    }
    setSmilesError(null);
    if (!SUPPORTED_SOLVENTS.includes(solvent)) {
      setError({
        title: "Solvent not supported",
        message:
          "Select one of the solvents present in the trained solvent encoder.",
      });
      return;
    }
    setStatus("running");
    setResult(null);
    try {
      const payload = await predictSolubility({
        smiles: smiles.trim(),
        solvent,
        temperatureK:
          tempNumber !== null && !Number.isNaN(tempNumber) ? tempNumber : null,
      });
      setResult(payload);
      setStatus("done");
    } catch (err) {
      const pe =
        err instanceof PredictionError
          ? err
          : new PredictionError(
              "Prediction failed",
              "An unexpected error occurred.",
            );
      setError({
        title: pe.title,
        message: pe.message,
      });
      setStatus("error");
    }
  }
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
    className:
      "grid gap-10 lg:grid-cols-[minmax(0,1fr)_minmax(0,0.95fr)] lg:gap-14",
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("form", {
        onSubmit: run,
        className: "space-y-9",
        noValidate: true,
        children: [
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("fieldset", {
            className: "space-y-3",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("legend", {
                className: "label-caps",
                children: "Molecular structure",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("label", {
                htmlFor: "smiles",
                className: "block text-sm text-muted-foreground",
                children: "SMILES notation for the solute",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
                id: "smiles",
                name: "smiles",
                value: smiles,
                onChange: (e) => {
                  setSmiles(e.target.value);
                  setSmilesError(null);
                },
                spellCheck: false,
                autoComplete: "off",
                placeholder: "CCO",
                "aria-invalid": smilesError ? true : void 0,
                "aria-describedby": smilesError
                  ? "smiles-error"
                  : "smiles-help",
                className:
                  "w-full border-b-2 border-input bg-transparent pb-2 font-mono text-2xl text-foreground transition-colors placeholder:text-muted-foreground/60 focus:border-primary focus:outline-none aria-[invalid=true]:border-destructive sm:text-3xl",
              }),
              smilesError
                ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("p", {
                    id: "smiles-error",
                    className: "text-sm text-destructive",
                    children: ["✕ ", smilesError],
                  })
                : /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("p", {
                    id: "smiles-help",
                    className: "flex flex-wrap items-center gap-2 text-sm",
                    children: [
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                        className: "text-muted-foreground",
                        children: "Example molecules:",
                      }),
                      EXAMPLE_MOLECULES.map((m) =>
                        /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                          "button",
                          {
                            type: "button",
                            onClick: () => setSmiles(m.smiles),
                            className:
                              "border border-border px-2 py-0.5 font-mono text-xs text-foreground transition-colors hover:border-primary hover:text-primary",
                            children: [
                              m.smiles,
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                "span",
                                {
                                  className:
                                    "ml-1.5 font-sans text-[0.65rem] text-muted-foreground",
                                  children: m.label,
                                },
                              ),
                            ],
                          },
                          m.smiles,
                        ),
                      ),
                    ],
                  }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("fieldset", {
            className: "space-y-4",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("legend", {
                className: "label-caps",
                children: [
                  "Solvent · ",
                  MODEL_METADATA.totalSolvents,
                  " supported",
                ],
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
                className:
                  "max-h-72 space-y-5 overflow-y-auto border border-border bg-card p-4",
                children: SOLVENT_GROUPS.map((group) =>
                  /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                    "div",
                    {
                      children: [
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                          className:
                            "font-mono text-[0.65rem] uppercase tracking-[0.14em] text-muted-foreground",
                          children: group.label,
                        }),
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
                          className: "mt-2 flex flex-wrap gap-1.5",
                          children: group.solvents.map((s) => {
                            const active = s === solvent;
                            return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                              "button",
                              {
                                type: "button",
                                onClick: () => setSolvent(s),
                                "aria-pressed": active,
                                className: active
                                  ? "border border-primary bg-primary px-2.5 py-1 text-xs text-primary-foreground"
                                  : "border border-border px-2.5 py-1 text-xs text-foreground transition-colors hover:border-primary/60 hover:text-primary",
                                children: [
                                  active &&
                                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                      "span",
                                      {
                                        "aria-hidden": "true",
                                        children: "✓ ",
                                      },
                                    ),
                                  s,
                                ],
                              },
                              s,
                            );
                          }),
                        }),
                      ],
                    },
                    group.label,
                  ),
                ),
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("fieldset", {
            className: "space-y-2",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("legend", {
                className: "label-caps",
                children: "Temperature — optional",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                className: "flex flex-wrap items-baseline gap-3",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
                    id: "temperature",
                    type: "number",
                    inputMode: "decimal",
                    step: "0.01",
                    value: temperature,
                    onChange: (e) => setTemperature(e.target.value),
                    placeholder: String(MODEL_METADATA.defaultTemperatureK),
                    "aria-describedby": "temp-help",
                    className:
                      "w-36 border-b border-input bg-transparent pb-1 font-mono text-lg focus:border-primary focus:outline-none",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                    className: "font-mono text-sm text-muted-foreground",
                    children: "K",
                  }),
                ],
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("p", {
                id: "temp-help",
                className: tempOutOfRange
                  ? "text-sm text-accent-foreground"
                  : "text-sm text-muted-foreground",
                children: [
                  "Training data spans ",
                  MODEL_METADATA.temperatureMinK,
                  "–",
                  MODEL_METADATA.temperatureMaxK,
                  " K. Left blank, the backend uses ",
                  MODEL_METADATA.defaultTemperatureK,
                  " K; values outside the range are clamped by the backend.",
                ],
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
            className:
              "flex flex-wrap items-center gap-4 border-t border-border pt-6",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
                type: "submit",
                disabled: status === "running",
                className:
                  "border border-primary bg-primary px-6 py-3 font-mono text-xs uppercase tracking-[0.18em] text-primary-foreground transition-opacity hover:opacity-90 disabled:opacity-60",
                children: status === "running" ? "Running…" : "Run prediction",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
                type: "button",
                onClick: () => setShowEndpoint((v) => !v),
                "aria-expanded": showEndpoint,
                className:
                  "font-mono text-xs uppercase tracking-[0.14em] text-muted-foreground underline decoration-border underline-offset-4 hover:text-primary",
                children: "Prediction service",
              }),
            ],
          }),
          showEndpoint &&
            /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
              className: "border border-border bg-secondary/50 p-4",
              children: [
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("label", {
                  htmlFor: "endpoint",
                  className: "label-caps",
                  children: "Flask backend base URL",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
                  id: "endpoint",
                  value: endpoint,
                  onChange: (e) => setEndpoint(e.target.value),
                  onBlur: () => setBackendUrl(endpoint),
                  spellCheck: false,
                  className:
                    "mt-2 w-full border-b border-input bg-transparent pb-1 font-mono text-sm focus:border-primary focus:outline-none",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("p", {
                  className:
                    "mt-3 text-sm leading-relaxed text-muted-foreground",
                  children: [
                    "Predictions are computed by this project's own Flask service (",
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", {
                      className: "font-mono",
                      children:
                        "solvent_project/model_training_BigSolDB/app.py",
                    }),
                    "), which loads the trained Random Forest model and solvent encoder. Run it locally and this interface posts to ",
                    /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("code", {
                      className: "font-mono",
                      children: [endpoint, "/predict"],
                    }),
                    ". No values are generated in the browser.",
                  ],
                }),
              ],
            }),
        ],
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
        className: "space-y-6",
        children: [
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                className: "label-caps mb-3",
                children: "Structure preview — rendered from your SMILES",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                className: "relative",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                    MoleculeStructure,
                    { smiles },
                  ),
                  status === "running" &&
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
                      className:
                        "pointer-events-none absolute inset-0 overflow-hidden",
                      children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                        "div",
                        { className: "animate-scan h-px w-full bg-primary/70" },
                      ),
                    }),
                ],
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(AnimatePresence, {
            mode: "wait",
            children: [
              status === "running" &&
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                  motion.div,
                  {
                    initial: {
                      opacity: 0,
                      y: 8,
                    },
                    animate: {
                      opacity: 1,
                      y: 0,
                    },
                    exit: { opacity: 0 },
                    className: "paper-panel p-6",
                    children: [
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                        className: "label-caps",
                        children: "Request in flight",
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("ul", {
                        className:
                          "mt-4 space-y-2 text-sm text-muted-foreground",
                        children: [
                          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("li", {
                            className: "flex items-center gap-3",
                            children: [
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                Dot,
                                {},
                              ),
                              " RDKit parses the structure and computes descriptors",
                            ],
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("li", {
                            className: "flex items-center gap-3",
                            children: [
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                Dot,
                                {},
                              ),
                              " Solvent one-hot encoding is appended to the feature row",
                            ],
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("li", {
                            className: "flex items-center gap-3",
                            children: [
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                Dot,
                                {},
                              ),
                              " Random Forest regressor estimates solubility",
                            ],
                          }),
                        ],
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                        className:
                          "mt-4 font-mono text-xs text-muted-foreground",
                        children:
                          "Awaiting the backend response — stages are the pipeline in app.py, not a progress simulation.",
                      }),
                    ],
                  },
                  "running",
                ),
              status === "error" &&
                error &&
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                  motion.div,
                  {
                    initial: {
                      opacity: 0,
                      y: 8,
                    },
                    animate: {
                      opacity: 1,
                      y: 0,
                    },
                    exit: { opacity: 0 },
                    role: "alert",
                    className:
                      "border border-destructive/40 bg-destructive/5 p-6",
                    children: [
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                        className: "font-display text-xl text-foreground",
                        children: error.title,
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                        className:
                          "mt-2 text-sm leading-relaxed text-muted-foreground",
                        children: error.message,
                      }),
                    ],
                  },
                  "error",
                ),
              status === "done" &&
                result &&
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                  motion.div,
                  {
                    initial: {
                      opacity: 0,
                      y: 10,
                    },
                    animate: {
                      opacity: 1,
                      y: 0,
                    },
                    transition: {
                      duration: 0.4,
                      ease: [0.22, 1, 0.36, 1],
                    },
                    className: "space-y-6",
                    children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                      ResultPanel,
                      { result },
                    ),
                  },
                  "result",
                ),
            ],
          }),
        ],
      }),
    ],
  });
}
function Dot() {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
    "aria-hidden": "true",
    className: "h-1.5 w-1.5 shrink-0 rounded-full bg-primary",
  });
}
function ResultPanel({ result }) {
  const value = result.predicted_solubility ?? 0;
  const suggestions = result.suggested_solvents ?? [];
  const lower = MODEL_METADATA.solubilityLowerBound;
  const upper = MODEL_METADATA.solubilityUpperBound;
  const position = Math.min(
    100,
    Math.max(0, ((value - lower) / (upper - lower)) * 100),
  );
  const thresholdPos =
    ((MODEL_METADATA.solubilityThreshold - lower) / (upper - lower)) * 100;
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
    import_jsx_runtime.Fragment,
    {
      children: [
        /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
          className: "paper-panel p-6 sm:p-8",
          "aria-live": "polite",
          children: [
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
              className: "label-caps",
              children: "Predicted solubility",
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
              className:
                "mt-2 font-mono text-5xl tracking-tight text-foreground sm:text-6xl",
              children: value.toFixed(3),
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
              className: "mt-1 font-mono text-sm text-muted-foreground",
              children: "log mol/L",
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
              className: "mt-8",
              children: [
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                  className:
                    "relative h-8 border border-border bg-secondary/60",
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
                      className: "absolute inset-y-0 left-0 bg-primary/20",
                      style: { width: `${position}%` },
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
                      className: "absolute inset-y-0 w-px bg-primary",
                      style: { left: `${position}%` },
                      "aria-hidden": "true",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
                      className:
                        "absolute inset-y-0 w-px border-l border-dashed border-accent-foreground/60",
                      style: { left: `${thresholdPos}%` },
                      "aria-hidden": "true",
                    }),
                  ],
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                  className:
                    "mt-2 flex justify-between font-mono text-[0.65rem] text-muted-foreground",
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                      children: lower.toFixed(1),
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", {
                      children: [
                        "threshold ",
                        MODEL_METADATA.solubilityThreshold.toFixed(1),
                      ],
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
                      children: upper.toFixed(1),
                    }),
                  ],
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("p", {
                  className: "mt-2 text-xs text-muted-foreground",
                  children: [
                    "Output range is clipped by the backend to [",
                    lower,
                    ", ",
                    upper,
                    "] log mol/L; the dashed mark is the repository's poor-solubility threshold.",
                  ],
                }),
              ],
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("dl", {
              className:
                "mt-8 grid gap-5 border-t border-border pt-6 sm:grid-cols-2",
              children: [
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", {
                      className: "label-caps",
                      children: "Compound",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dd", {
                      className: "mt-1 break-words text-sm text-foreground",
                      children: result.solute_name,
                    }),
                  ],
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", {
                      className: "label-caps",
                      children: "Solvent",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dd", {
                      className: "mt-1 text-sm text-foreground",
                      children: result.solvent_name,
                    }),
                  ],
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", {
                      className: "label-caps",
                      children: "Temperature used",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("dd", {
                      className: "mt-1 font-mono text-sm text-foreground",
                      children: [result.temperature_used_K, " K"],
                    }),
                  ],
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", {
                      className: "label-caps",
                      children: "Model",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dd", {
                      className: "mt-1 text-sm text-foreground",
                      children: "Random Forest Regressor",
                    }),
                  ],
                }),
              ],
            }),
            result.note &&
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                className:
                  "mt-6 border-l-2 border-accent-foreground/50 bg-accent/40 px-4 py-3 text-sm text-accent-foreground",
                children: result.note,
              }),
            result.message &&
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                className: "mt-4 text-sm leading-relaxed text-muted-foreground",
                children: result.message,
              }),
          ],
        }),
        suggestions.length > 0 &&
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
            className: "paper-panel p-6 sm:p-8",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
                className: "font-display text-2xl text-foreground",
                children: "Looking for a better solvent?",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                className: "mt-2 text-sm leading-relaxed text-muted-foreground",
                children:
                  "The backend re-ran the same molecule through every remaining solvent in the encoder and returned those that clear its solubility threshold.",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("ul", {
                className: "mt-6 space-y-3",
                children: suggestions.map((s) => {
                  const width = Math.min(
                    100,
                    Math.max(
                      2,
                      ((s.predicted_solubility - lower) / (upper - lower)) *
                        100,
                    ),
                  );
                  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                    "li",
                    {
                      children: [
                        /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                          className:
                            "flex items-baseline justify-between gap-4",
                          children: [
                            /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                              "span",
                              {
                                className: "text-sm text-foreground",
                                children: s.solvent,
                              },
                            ),
                            /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                              "span",
                              {
                                className: "font-mono text-sm text-foreground",
                                children: [
                                  s.predicted_solubility.toFixed(3),
                                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                    "span",
                                    {
                                      className:
                                        "ml-1 text-xs text-muted-foreground",
                                      children: "log mol/L",
                                    },
                                  ),
                                ],
                              },
                            ),
                          ],
                        }),
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
                          className: "mt-1.5 h-1.5 bg-secondary",
                          children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                            "div",
                            {
                              className: "h-full bg-primary/70",
                              style: { width: `${width}%` },
                            },
                          ),
                        }),
                      ],
                    },
                    s.solvent,
                  );
                }),
              }),
            ],
          }),
      ],
    },
  );
}
var STAGES = [
  {
    label: "Molecular structure",
    detail:
      "A SMILES string identifies the solute. RDKit parses it into a molecule object.",
  },
  {
    label: "RDKit",
    detail:
      "The cheminformatics layer. If parsing fails, the request is rejected as an invalid structure rather than guessed.",
  },
  {
    label: "Molecular descriptors + Morgan fingerprint",
    detail:
      "Eight physicochemical descriptors describe bulk properties; a 1024-bit Morgan fingerprint (radius 2) encodes local substructures.",
  },
  {
    label: "Solvent encoding",
    detail:
      "The selected solvent is one-hot encoded with the scikit-learn encoder fitted during training, so only trained solvents are accepted.",
  },
  {
    label: "Random Forest",
    detail:
      "A 500-tree regressor trained on the combined feature matrix maps structure and solvent to a solubility value.",
  },
  {
    label: "Solubility prediction",
    detail:
      "A continuous value in log mol/L, clipped to the bounds recorded in the model metadata, plus alternative solvents when solubility is poor.",
  },
];
function Home() {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
    className: "min-h-screen scroll-smooth",
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)(SiteHeader, {}),
      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("main", {
        children: [
          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("section", {
            id: "home",
            className:
              "mx-auto max-w-6xl px-5 pb-16 pt-12 sm:px-8 sm:pt-20 scroll-mt-20",
            children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
              className:
                "grid items-center gap-14 lg:grid-cols-[minmax(0,1.05fr)_minmax(0,1fr)]",
              children: [
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                      className: "label-caps",
                      children: "Computational solubility · research prototype",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h1", {
                      className:
                        "mt-5 text-4xl leading-[1.08] text-foreground sm:text-5xl lg:text-6xl",
                      children: "Understand how molecules dissolve.",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                      className:
                        "mt-6 max-w-xl text-base leading-relaxed text-muted-foreground sm:text-lg",
                      children:
                        "DissolveAI estimates the solubility of a compound in a chosen solvent from its molecular structure alone. A SMILES string is converted into cheminformatics features with RDKit, combined with an encoded solvent identity, and passed to a trained Random Forest regressor.",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                      className: "mt-9 flex flex-wrap items-center gap-4",
                      children: [
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("a", {
                          href: "#predictor",
                          className:
                            "border border-primary bg-primary px-6 py-3 font-mono text-xs uppercase tracking-[0.18em] text-primary-foreground transition-opacity hover:opacity-90",
                          children: "Explore predictor",
                        }),
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("a", {
                          href: "#how-it-works",
                          className:
                            "border border-input px-6 py-3 font-mono text-xs uppercase tracking-[0.18em] text-foreground transition-colors hover:border-primary hover:text-primary",
                          children: "How it works",
                        }),
                      ],
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("dl", {
                      className:
                        "mt-12 grid grid-cols-2 gap-x-8 gap-y-6 border-t border-border pt-8 sm:grid-cols-3",
                      children: [
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Stat, {
                          value: MODEL_METADATA.totalSolvents.toString(),
                          label: "Solvents in the encoder",
                        }),
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Stat, {
                          value: DATASET_ROWS.toLocaleString("en-US"),
                          label: "Cleaned training rows",
                        }),
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Stat, {
                          value: `${MODEL_METADATA.temperatureMinK}–${MODEL_METADATA.temperatureMaxK} K`,
                          label: "Temperature coverage",
                        }),
                      ],
                    }),
                  ],
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)(HeroLab, {}),
              ],
            }),
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("section", {
            className: "border-y border-border bg-secondary/50",
            children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
              className:
                "mx-auto grid max-w-6xl gap-10 px-5 py-16 sm:px-8 md:grid-cols-3",
              children: [
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Card, {
                  index: "01",
                  title: "Structure in",
                  body: "A SMILES string is parsed by RDKit. Physicochemical descriptors and a 1024-bit Morgan fingerprint describe the molecule numerically.",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Card, {
                  index: "02",
                  title: "Solvent context",
                  body: "The chosen solvent is one-hot encoded from the vocabulary the model was trained on, so solvent identity is part of the feature row.",
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Card, {
                  index: "03",
                  title: "Estimate out",
                  body: "A Random Forest regressor returns a continuous solubility value in log mol/L, with alternative solvents suggested when solubility is poor.",
                }),
              ],
            }),
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
            id: "predictor",
            className:
              "mx-auto max-w-6xl px-5 py-20 sm:px-8 scroll-mt-20 border-b border-border",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("header", {
                className: "max-w-2xl mb-12",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className: "label-caps",
                    children: "Experiment bench",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", {
                    className: "mt-4 text-3xl text-foreground sm:text-4xl",
                    children: "Run a solubility prediction",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className:
                      "mt-4 text-base leading-relaxed text-muted-foreground",
                    children:
                      "Enter a molecular structure and select a solvent to estimate its predicted solubility. Every number shown below comes from this project's trained model, served by its Flask prediction service.",
                  }),
                ],
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Predictor, {}),
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                className: "mt-12 border-t border-border pt-8",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
                    className: "text-xl text-foreground mb-3",
                    children: "Molecular features used by the request",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className:
                      "max-w-2xl text-sm leading-relaxed text-muted-foreground",
                    children:
                      "The deployed prediction service computes molecular weight, LogP and TPSA with RDKit, appends the one-hot encoded solvent, and predicts with the trained Random Forest. Those descriptor values are used internally. The full training pipeline additionally derives H-bond donors and acceptors, rotatable bonds, aromatic rings, fraction sp³ carbon, and a 1024-bit Morgan fingerprint.",
                  }),
                ],
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
            id: "how-it-works",
            className:
              "mx-auto max-w-6xl px-5 py-20 sm:px-8 scroll-mt-20 border-b border-border",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("header", {
                className: "max-w-2xl",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className: "label-caps",
                    children: "Pipeline",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", {
                    className: "mt-4 text-3xl text-foreground sm:text-4xl",
                    children: "From molecule to solubility",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className:
                      "mt-4 text-base leading-relaxed text-muted-foreground",
                    children:
                      "Solubility prediction is treated as a regression problem. The inputs are a molecular structure and a solvent identity; the output is a continuous solubility value.",
                  }),
                ],
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("ol", {
                className: "mt-14 space-y-0",
                children: STAGES.map((stage, i) =>
                  /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                    "li",
                    {
                      className:
                        "relative grid gap-4 pb-10 sm:grid-cols-[7rem_1fr]",
                      children: [
                        /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                          className: "flex items-start gap-4",
                          children: [
                            /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                              "span",
                              {
                                className: "font-mono text-xs text-primary",
                                children: String(i + 1).padStart(2, "0"),
                              },
                            ),
                            i < STAGES.length - 1 &&
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                "span",
                                {
                                  "aria-hidden": "true",
                                  className:
                                    "absolute left-[0.4rem] top-6 h-full w-px bg-border",
                                },
                              ),
                          ],
                        }),
                        /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                          className: "border-b border-border pb-6",
                          children: [
                            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
                              className: "text-xl text-foreground",
                              children: stage.label,
                            }),
                            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                              className:
                                "mt-2 max-w-2xl text-sm leading-relaxed text-muted-foreground",
                              children: stage.detail,
                            }),
                          ],
                        }),
                      ],
                    },
                    stage.label,
                  ),
                ),
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                className: "mt-16 pt-12",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
                    className: "text-2xl text-foreground sm:text-3xl",
                    children: "Feature engineering",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className:
                      "mt-3 max-w-2xl text-base leading-relaxed text-muted-foreground",
                    children:
                      "The model does not read a SMILES string directly. The structure is first transformed into numerical features, and those features are what the regressor sees.",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                    className: "mt-10 grid gap-8 md:grid-cols-3",
                    children: [
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Feature, {
                        title: "Molecular descriptors",
                        body: "Molecular weight, LogP, TPSA, hydrogen-bond donors and acceptors, rotatable bonds, aromatic rings, and fraction of sp³ carbon — eight properties describing overall molecular character.",
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Feature, {
                        title: "Morgan fingerprints",
                        body: "A 1024-bit vector at radius 2 recording which circular substructures are present, giving the model structural detail that bulk descriptors miss.",
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Feature, {
                        title: "Solvent encoding",
                        body: `One-hot encoding over the ${MODEL_METADATA.totalSolvents} solvents in the training data, so the same molecule can be evaluated in different chemical environments.`,
                      }),
                    ],
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className:
                      "mt-8 max-w-2xl text-sm leading-relaxed text-muted-foreground",
                    children:
                      "These blocks are concatenated into a single feature matrix and passed to the trained Random Forest model. Temperature is carried through as a numerical feature, which is why the predictor exposes an optional temperature field.",
                  }),
                ],
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
            id: "model",
            className:
              "mx-auto max-w-6xl px-5 py-20 sm:px-8 scroll-mt-20 border-b border-border",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("header", {
                className: "max-w-2xl",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className: "label-caps",
                    children: "Model card",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", {
                    className: "mt-4 text-3xl text-foreground sm:text-4xl",
                    children: "The model behind DissolveAI",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className:
                      "mt-4 text-base leading-relaxed text-muted-foreground",
                    children:
                      "A supervised regression model, trained offline and served by a Flask prediction endpoint. Everything below is taken from the project's training code and model metadata; no performance figures are claimed beyond what the repository records.",
                  }),
                ],
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("dl", {
                className:
                  "mt-12 grid gap-x-10 gap-y-8 border-t border-border pt-10 sm:grid-cols-2 lg:grid-cols-4",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                    label: "Algorithm",
                    value: "Random Forest Regressor",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                    label: "Trees",
                    value: "500",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                    label: "Split",
                    value: "80 / 20 train–test",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                    label: "Validation",
                    value: "5-fold cross-validation",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                    label: "Evaluation metric",
                    value: "RMSE",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                    label: "Target unit",
                    value: "log mol/L",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                    label: "Prediction bounds",
                    value: `${MODEL_METADATA.solubilityLowerBound} to ${MODEL_METADATA.solubilityUpperBound}`,
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                    label: "Poor-solubility threshold",
                    value: MODEL_METADATA.solubilityThreshold.toString(),
                  }),
                ],
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                className:
                  "mt-6 max-w-2xl text-sm leading-relaxed text-muted-foreground",
                children:
                  "Cross-validated RMSE values are produced when the training script is run and are not published as fixed numbers in the repository, so no accuracy figure is shown here.",
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                className: "mt-16 pt-12 border-t border-border",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
                    className: "text-2xl text-foreground sm:text-3xl",
                    children: "Dataset",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                    className:
                      "mt-6 grid gap-10 lg:grid-cols-[minmax(0,1.1fr)_minmax(0,1fr)]",
                    children: [
                      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                        className:
                          "space-y-5 text-base leading-relaxed text-muted-foreground",
                        children: [
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                            children:
                              "The served model is trained on the cleaned BigSolDB solubility dataset: measured solubility values for solute–solvent pairs at recorded temperatures. Cleaning is handled by the project's own data-processing scripts before feature generation.",
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                            children:
                              "An additional AqSolDB pipeline exists in the repository for aqueous solubility work; the prediction service used by this interface is the BigSolDB model, which is what makes multi-solvent prediction possible.",
                          }),
                        ],
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("dl", {
                        className:
                          "grid gap-6 border-t border-rule-strong pt-6 sm:grid-cols-2 lg:border-l lg:border-t-0 lg:pl-8 lg:pt-0",
                        children: [
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                            label: "Source",
                            value: "BigSolDB (cleaned)",
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                            label: "Rows after cleaning",
                            value: DATASET_ROWS.toLocaleString("en-US"),
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                            label: "Solvents",
                            value: MODEL_METADATA.totalSolvents.toString(),
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                            label: "Temperature range",
                            value: `${MODEL_METADATA.temperatureMinK} – ${MODEL_METADATA.temperatureMaxK} K`,
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                            label: "Measured unit",
                            value: "mol/L",
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Spec, {
                            label: "Held-out test share",
                            value: "20%",
                          }),
                        ],
                      }),
                    ],
                  }),
                ],
              }),
            ],
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("section", {
            id: "limitations",
            className:
              "mx-auto max-w-6xl px-5 py-20 sm:px-8 border-b border-border",
            children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
              className: "grid gap-12 md:grid-cols-2",
              children: [
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", {
                      className: "text-2xl text-foreground sm:text-3xl",
                      children: "Why solubility matters",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                      className:
                        "mt-6 space-y-5 text-base leading-relaxed text-muted-foreground",
                      children: [
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                          children:
                            "Solubility determines whether a compound can be processed, formulated, dosed, or even measured. In drug development it governs bioavailability and the choice of a formulation vehicle. In chemical formulation and materials research it decides which solvent system a synthesis or coating can use at all.",
                        }),
                        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                          children:
                            "Measuring solubility experimentally for every candidate compound and solvent pair is slow and material-intensive. A computational estimate narrows the search space before anything is weighed out, which is the purpose DissolveAI is built for: fast screening guidance for early-stage research, not a replacement for measurement.",
                        }),
                      ],
                    }),
                  ],
                }),
                /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                  children: [
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", {
                      className: "text-2xl text-foreground sm:text-3xl",
                      children: "Limitations",
                    }),
                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)("ul", {
                      className:
                        "mt-6 space-y-4 text-sm leading-relaxed text-muted-foreground",
                      children: [
                        "Performance depends on the coverage of the training dataset. Chemistry that is sparsely represented in BigSolDB will be predicted less reliably.",
                        "Extrapolation beyond the training chemical space is limited; the model interpolates patterns it has seen rather than reasoning about new chemistry.",
                        "The approach is purely data-driven, with no explicit thermodynamic modelling of dissolution.",
                        "Only solvents present in the fitted encoder can be evaluated; any other solvent is refused rather than approximated.",
                        `Temperatures outside ${MODEL_METADATA.temperatureMinK}–${MODEL_METADATA.temperatureMaxK} K are clamped to that range by the prediction service.`,
                      ].map((item) =>
                        /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                          "li",
                          {
                            className: "flex gap-3",
                            children: [
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                "span",
                                {
                                  "aria-hidden": "true",
                                  className:
                                    "mt-2 h-1 w-3 shrink-0 bg-accent-foreground/60",
                                },
                              ),
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                "span",
                                { children: item },
                              ),
                            ],
                          },
                          item,
                        ),
                      ),
                    }),
                  ],
                }),
              ],
            }),
          }),
          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
            id: "about",
            className: "mx-auto max-w-6xl px-5 py-20 sm:px-8 scroll-mt-20",
            children: [
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("header", {
                className: "max-w-2xl mb-10",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                    className: "label-caps",
                    children: "Project",
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", {
                    className: "mt-4 text-3xl text-foreground sm:text-4xl",
                    children: "About DissolveAI",
                  }),
                ],
              }),
              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                className:
                  "grid gap-14 lg:grid-cols-[minmax(0,1.15fr)_minmax(0,0.85fr)]",
                children: [
                  /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                    className:
                      "space-y-6 text-base leading-relaxed text-muted-foreground",
                    children: [
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                        children:
                          "DissolveAI is a supervised machine-learning system for predicting the solubility of chemical compounds in a range of solvents at room temperature. It was built to support early-stage research, where the number of candidate compound–solvent combinations is far larger than the number that can reasonably be measured in the lab.",
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                        children:
                          "The system pairs cheminformatics with regression modelling. Molecular structures given as SMILES are converted into physicochemical descriptors and Morgan fingerprints with RDKit; the solvent is one-hot encoded; a Random Forest regressor trained on cleaned BigSolDB measurements maps the combined features to a solubility value. When the predicted solubility in the requested solvent is poor, the service evaluates the same molecule in the remaining trained solvents and reports better options.",
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
                        children:
                          "The architecture is deliberately separated: an offline training pipeline produces the model and encoder artefacts, a Flask service loads them and exposes a single prediction endpoint, and this web interface is a thin client over that endpoint. No prediction logic lives in the browser, so what you see in the predictor is exactly what the model returns.",
                      }),
                    ],
                  }),
                  /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("aside", {
                    className: "space-y-8",
                    children: [
                      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                        className: "border-t border-rule-strong pt-5",
                        children: [
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
                            className: "label-caps",
                            children: "Technologies",
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("ul", {
                            className: "mt-3 flex flex-wrap gap-2",
                            children: [
                              "Python",
                              "RDKit",
                              "scikit-learn",
                              "Flask",
                              "NumPy",
                              "pandas",
                              "React",
                            ].map((t) =>
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                "li",
                                {
                                  className:
                                    "border border-border px-2.5 py-1 font-mono text-xs text-foreground",
                                  children: t,
                                },
                                t,
                              ),
                            ),
                          }),
                        ],
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                        className: "border-t border-rule-strong pt-5",
                        children: [
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
                            className: "label-caps",
                            children: "Problem formulation",
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("dl", {
                            className: "mt-3 space-y-3 text-sm",
                            children: [
                              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                                "div",
                                {
                                  children: [
                                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                      "dt",
                                      {
                                        className: "text-muted-foreground",
                                        children: "Inputs",
                                      },
                                    ),
                                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                      "dd",
                                      {
                                        className: "text-foreground",
                                        children:
                                          "Molecular structure (SMILES), solvent identity",
                                      },
                                    ),
                                  ],
                                },
                              ),
                              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                                "div",
                                {
                                  children: [
                                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                      "dt",
                                      {
                                        className: "text-muted-foreground",
                                        children: "Output",
                                      },
                                    ),
                                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                      "dd",
                                      {
                                        className: "text-foreground",
                                        children:
                                          "Continuous solubility value (log mol/L)",
                                      },
                                    ),
                                  ],
                                },
                              ),
                              /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(
                                "div",
                                {
                                  children: [
                                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                      "dt",
                                      {
                                        className: "text-muted-foreground",
                                        children: "Task",
                                      },
                                    ),
                                    /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                      "dd",
                                      {
                                        className: "text-foreground",
                                        children: "Regression",
                                      },
                                    ),
                                  ],
                                },
                              ),
                            ],
                          }),
                        ],
                      }),
                      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
                        className: "border-t border-rule-strong pt-5",
                        children: [
                          /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
                            className: "label-caps",
                            children: "Source",
                          }),
                          /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("ul", {
                            className: "mt-3 space-y-2 text-sm",
                            children: [
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                "li",
                                {
                                  children: /* @__PURE__ */ (0,
                                  import_jsx_runtime.jsx)("a", {
                                    href: REPO_URL,
                                    target: "_blank",
                                    rel: "noreferrer noopener",
                                    className:
                                      "text-foreground underline decoration-border underline-offset-4 hover:text-primary",
                                    children: "Project repository",
                                  }),
                                },
                              ),
                              /* @__PURE__ */ (0, import_jsx_runtime.jsx)(
                                "li",
                                {
                                  children: /* @__PURE__ */ (0,
                                  import_jsx_runtime.jsx)("a", {
                                    href: PROFILE_URL,
                                    target: "_blank",
                                    rel: "noreferrer noopener",
                                    className:
                                      "text-foreground underline decoration-border underline-offset-4 hover:text-primary",
                                    children: "GitHub profile",
                                  }),
                                },
                              ),
                            ],
                          }),
                        ],
                      }),
                    ],
                  }),
                ],
              }),
            ],
          }),
        ],
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)(SiteFooter, {}),
    ],
  });
}
function Stat({ value, label }) {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", {
        className: "label-caps",
        children: label,
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dd", {
        className: "mt-1 font-mono text-2xl text-foreground",
        children: value,
      }),
    ],
  });
}
function Card({ index, title, body }) {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("article", {
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
        className: "font-mono text-xs text-primary",
        children: index,
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
        className: "mt-3 text-xl text-foreground",
        children: title,
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
        className: "mt-3 text-sm leading-relaxed text-muted-foreground",
        children: body,
      }),
    ],
  });
}
function Feature({ title, body }) {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("article", {
    className: "border-t border-rule-strong pt-4",
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h3", {
        className: "text-lg text-foreground",
        children: title,
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
        className: "mt-2 text-sm leading-relaxed text-muted-foreground",
        children: body,
      }),
    ],
  });
}
function Spec({ label, value }) {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", {
        className: "label-caps",
        children: label,
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dd", {
        className: "mt-1 font-mono text-base text-foreground",
        children: value,
      }),
    ],
  });
}
//#endregion
export { Home as component };
