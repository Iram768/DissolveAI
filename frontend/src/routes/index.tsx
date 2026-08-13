import { createFileRoute } from "@tanstack/react-router";
import { SiteHeader } from "@/components/site-header";
import { SiteFooter } from "@/components/site-footer";
import { HeroLab } from "@/components/hero-lab";
import { Predictor } from "@/components/predictor";
import {
  DATASET_ROWS,
  MODEL_METADATA,
  PROFILE_URL,
  REPO_URL,
} from "@/lib/dissolveai";

export const Route = createFileRoute("/")({
  head: () => ({
    meta: [
      { title: "DissolveAI — Molecular Solubility Prediction" },
      {
        name: "description",
        content:
          "DissolveAI estimates compound solubility across 70 solvents from SMILES structures using RDKit features and a Random Forest regressor.",
      },
      {
        property: "og:title",
        content: "DissolveAI — Molecular Solubility Prediction",
      },
      {
        property: "og:description",
        content:
          "Predict molecular solubility before the experiment: RDKit descriptors, Morgan fingerprints, solvent encoding, Random Forest regression.",
      },
    ],
  }),
  component: Home,
});

const STAGES = [
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
  return (
    <div className="min-h-screen scroll-smooth">
      <SiteHeader />
      <main>
        {/* HERO SECTION */}
        <section
          id="home"
          className="mx-auto max-w-6xl px-5 pb-16 pt-12 sm:px-8 sm:pt-20 scroll-mt-20"
        >
          <div className="grid items-center gap-14 lg:grid-cols-[minmax(0,1.05fr)_minmax(0,1fr)]">
            <div>
              <p className="label-caps">
                Computational solubility · research prototype
              </p>
              <h1 className="mt-5 text-4xl leading-[1.08] text-foreground sm:text-5xl lg:text-6xl">
                Understand how molecules dissolve.
              </h1>
              <p className="mt-6 max-w-xl text-base leading-relaxed text-muted-foreground sm:text-lg">
                DissolveAI estimates the solubility of a compound in a chosen
                solvent from its molecular structure alone. A SMILES string is
                converted into cheminformatics features with RDKit, combined
                with an encoded solvent identity, and passed to a trained Random
                Forest regressor.
              </p>
              <div className="mt-9 flex flex-wrap items-center gap-4">
                <a
                  href="#predictor"
                  className="border border-primary bg-primary px-6 py-3 font-mono text-xs uppercase tracking-[0.18em] text-primary-foreground transition-opacity hover:opacity-90"
                >
                  Explore predictor
                </a>
                <a
                  href="#how-it-works"
                  className="border border-input px-6 py-3 font-mono text-xs uppercase tracking-[0.18em] text-foreground transition-colors hover:border-primary hover:text-primary"
                >
                  How it works
                </a>
              </div>

              <dl className="mt-12 grid grid-cols-2 gap-x-8 gap-y-6 border-t border-border pt-8 sm:grid-cols-3">
                <Stat
                  value={MODEL_METADATA.totalSolvents.toString()}
                  label="Solvents in the encoder"
                />
                <Stat
                  value={DATASET_ROWS.toLocaleString("en-US")}
                  label="Cleaned training rows"
                />
                <Stat
                  value={`${MODEL_METADATA.temperatureMinK}–${MODEL_METADATA.temperatureMaxK} K`}
                  label="Temperature coverage"
                />
              </dl>
            </div>

            <HeroLab />
          </div>
        </section>

        {/* WORKFLOW SUMMARY CARDS */}
        <section className="border-y border-border bg-secondary/50">
          <div className="mx-auto grid max-w-6xl gap-10 px-5 py-16 sm:px-8 md:grid-cols-3">
            <Card
              index="01"
              title="Structure in"
              body="A SMILES string is parsed by RDKit. Physicochemical descriptors and a 1024-bit Morgan fingerprint describe the molecule numerically."
            />
            <Card
              index="02"
              title="Solvent context"
              body="The chosen solvent is one-hot encoded from the vocabulary the model was trained on, so solvent identity is part of the feature row."
            />
            <Card
              index="03"
              title="Estimate out"
              body="A Random Forest regressor returns a continuous solubility value in log mol/L, with alternative solvents suggested when solubility is poor."
            />
          </div>
        </section>

        {/* PREDICTOR SECTION */}
        <section
          id="predictor"
          className="mx-auto max-w-6xl px-5 py-20 sm:px-8 scroll-mt-20 border-b border-border"
        >
          <header className="max-w-2xl mb-12">
            <p className="label-caps">Experiment bench</p>
            <h2 className="mt-4 text-3xl text-foreground sm:text-4xl">
              Run a solubility prediction
            </h2>
            <p className="mt-4 text-base leading-relaxed text-muted-foreground">
              Enter a molecular structure and select a solvent to estimate its
              predicted solubility. Every number shown below comes from this
              project&apos;s trained model, served by its Flask prediction
              service.
            </p>
          </header>

          <Predictor />

          <div className="mt-12 border-t border-border pt-8">
            <h3 className="text-xl text-foreground mb-3">
              Molecular features used by the request
            </h3>
            <p className="max-w-2xl text-sm leading-relaxed text-muted-foreground">
              The deployed prediction service computes molecular weight, LogP
              and TPSA with RDKit, appends the one-hot encoded solvent, and
              predicts with the trained Random Forest. Those descriptor values
              are used internally. The full training pipeline additionally
              derives H-bond donors and acceptors, rotatable bonds, aromatic
              rings, fraction sp³ carbon, and a 1024-bit Morgan fingerprint.
            </p>
          </div>
        </section>

        {/* HOW IT WORKS SECTION */}
        <section
          id="how-it-works"
          className="mx-auto max-w-6xl px-5 py-20 sm:px-8 scroll-mt-20 border-b border-border"
        >
          <header className="max-w-2xl">
            <p className="label-caps">Pipeline</p>
            <h2 className="mt-4 text-3xl text-foreground sm:text-4xl">
              From molecule to solubility
            </h2>
            <p className="mt-4 text-base leading-relaxed text-muted-foreground">
              Solubility prediction is treated as a regression problem. The
              inputs are a molecular structure and a solvent identity; the
              output is a continuous solubility value.
            </p>
          </header>

          <ol className="mt-14 space-y-0">
            {STAGES.map((stage, i) => (
              <li
                key={stage.label}
                className="relative grid gap-4 pb-10 sm:grid-cols-[7rem_1fr]"
              >
                <div className="flex items-start gap-4">
                  <span className="font-mono text-xs text-primary">
                    {String(i + 1).padStart(2, "0")}
                  </span>
                  {i < STAGES.length - 1 && (
                    <span
                      aria-hidden="true"
                      className="absolute left-[0.4rem] top-6 h-full w-px bg-border"
                    />
                  )}
                </div>
                <div className="border-b border-border pb-6">
                  <h3 className="text-xl text-foreground">{stage.label}</h3>
                  <p className="mt-2 max-w-2xl text-sm leading-relaxed text-muted-foreground">
                    {stage.detail}
                  </p>
                </div>
              </li>
            ))}
          </ol>

          <div className="mt-16 pt-12">
            <h3 className="text-2xl text-foreground sm:text-3xl">
              Feature engineering
            </h3>
            <p className="mt-3 max-w-2xl text-base leading-relaxed text-muted-foreground">
              The model does not read a SMILES string directly. The structure is
              first transformed into numerical features, and those features are
              what the regressor sees.
            </p>
            <div className="mt-10 grid gap-8 md:grid-cols-3">
              <Feature
                title="Molecular descriptors"
                body="Molecular weight, LogP, TPSA, hydrogen-bond donors and acceptors, rotatable bonds, aromatic rings, and fraction of sp³ carbon — eight properties describing overall molecular character."
              />
              <Feature
                title="Morgan fingerprints"
                body="A 1024-bit vector at radius 2 recording which circular substructures are present, giving the model structural detail that bulk descriptors miss."
              />
              <Feature
                title="Solvent encoding"
                body={`One-hot encoding over the ${MODEL_METADATA.totalSolvents} solvents in the training data, so the same molecule can be evaluated in different chemical environments.`}
              />
            </div>
            <p className="mt-8 max-w-2xl text-sm leading-relaxed text-muted-foreground">
              These blocks are concatenated into a single feature matrix and
              passed to the trained Random Forest model. Temperature is carried
              through as a numerical feature, which is why the predictor exposes
              an optional temperature field.
            </p>
          </div>
        </section>

        {/* MODEL SECTION */}
        <section
          id="model"
          className="mx-auto max-w-6xl px-5 py-20 sm:px-8 scroll-mt-20 border-b border-border"
        >
          <header className="max-w-2xl">
            <p className="label-caps">Model card</p>
            <h2 className="mt-4 text-3xl text-foreground sm:text-4xl">
              The model behind DissolveAI
            </h2>
            <p className="mt-4 text-base leading-relaxed text-muted-foreground">
              A supervised regression model, trained offline and served by a
              Flask prediction endpoint. Everything below is taken from the
              project&apos;s training code and model metadata; no performance
              figures are claimed beyond what the repository records.
            </p>
          </header>

          <dl className="mt-12 grid gap-x-10 gap-y-8 border-t border-border pt-10 sm:grid-cols-2 lg:grid-cols-4">
            <Spec label="Algorithm" value="Random Forest Regressor" />
            <Spec label="Trees" value="500" />
            <Spec label="Split" value="80 / 20 train–test" />
            <Spec label="Validation" value="5-fold cross-validation" />
            <Spec label="Evaluation metric" value="RMSE" />
            <Spec label="Target unit" value="log mol/L" />
            <Spec
              label="Prediction bounds"
              value={`${MODEL_METADATA.solubilityLowerBound} to ${MODEL_METADATA.solubilityUpperBound}`}
            />
            <Spec
              label="Poor-solubility threshold"
              value={MODEL_METADATA.solubilityThreshold.toString()}
            />
          </dl>
          <p className="mt-6 max-w-2xl text-sm leading-relaxed text-muted-foreground">
            Cross-validated RMSE values are produced when the training script is
            run and are not published as fixed numbers in the repository, so no
            accuracy figure is shown here.
          </p>

          <div className="mt-16 pt-12 border-t border-border">
            <h3 className="text-2xl text-foreground sm:text-3xl">Dataset</h3>
            <div className="mt-6 grid gap-10 lg:grid-cols-[minmax(0,1.1fr)_minmax(0,1fr)]">
              <div className="space-y-5 text-base leading-relaxed text-muted-foreground">
                <p>
                  The served model is trained on the cleaned BigSolDB solubility
                  dataset: measured solubility values for solute–solvent pairs
                  at recorded temperatures. Cleaning is handled by the
                  project&apos;s own data-processing scripts before feature
                  generation.
                </p>
                <p>
                  An additional AqSolDB pipeline exists in the repository for
                  aqueous solubility work; the prediction service used by this
                  interface is the BigSolDB model, which is what makes
                  multi-solvent prediction possible.
                </p>
              </div>
              <dl className="grid gap-6 border-t border-rule-strong pt-6 sm:grid-cols-2 lg:border-l lg:border-t-0 lg:pl-8 lg:pt-0">
                <Spec label="Source" value="BigSolDB (cleaned)" />
                <Spec
                  label="Rows after cleaning"
                  value={DATASET_ROWS.toLocaleString("en-US")}
                />
                <Spec
                  label="Solvents"
                  value={MODEL_METADATA.totalSolvents.toString()}
                />
                <Spec
                  label="Temperature range"
                  value={`${MODEL_METADATA.temperatureMinK} – ${MODEL_METADATA.temperatureMaxK} K`}
                />
                <Spec label="Measured unit" value="mol/L" />
                <Spec label="Held-out test share" value="20%" />
              </dl>
            </div>
          </div>
        </section>

        {/* LIMITATIONS & WHY SOLUBILITY MATTERS */}
        <section
          id="limitations"
          className="mx-auto max-w-6xl px-5 py-20 sm:px-8 border-b border-border"
        >
          <div className="grid gap-12 md:grid-cols-2">
            <div>
              <h2 className="text-2xl text-foreground sm:text-3xl">
                Why solubility matters
              </h2>
              <div className="mt-6 space-y-5 text-base leading-relaxed text-muted-foreground">
                <p>
                  Solubility determines whether a compound can be processed,
                  formulated, dosed, or even measured. In drug development it
                  governs bioavailability and the choice of a formulation
                  vehicle. In chemical formulation and materials research it
                  decides which solvent system a synthesis or coating can use at
                  all.
                </p>
                <p>
                  Measuring solubility experimentally for every candidate
                  compound and solvent pair is slow and material-intensive. A
                  computational estimate narrows the search space before
                  anything is weighed out, which is the purpose DissolveAI is
                  built for: fast screening guidance for early-stage research,
                  not a replacement for measurement.
                </p>
              </div>
            </div>
            <div>
              <h2 className="text-2xl text-foreground sm:text-3xl">
                Limitations
              </h2>
              <ul className="mt-6 space-y-4 text-sm leading-relaxed text-muted-foreground">
                {[
                  "Performance depends on the coverage of the training dataset. Chemistry that is sparsely represented in BigSolDB will be predicted less reliably.",
                  "Extrapolation beyond the training chemical space is limited; the model interpolates patterns it has seen rather than reasoning about new chemistry.",
                  "The approach is purely data-driven, with no explicit thermodynamic modelling of dissolution.",
                  "Only solvents present in the fitted encoder can be evaluated; any other solvent is refused rather than approximated.",
                  `Temperatures outside ${MODEL_METADATA.temperatureMinK}–${MODEL_METADATA.temperatureMaxK} K are clamped to that range by the prediction service.`,
                ].map((item) => (
                  <li key={item} className="flex gap-3">
                    <span
                      aria-hidden="true"
                      className="mt-2 h-1 w-3 shrink-0 bg-accent-foreground/60"
                    />
                    <span>{item}</span>
                  </li>
                ))}
              </ul>
            </div>
          </div>
        </section>

        {/* ABOUT SECTION */}
        <section
          id="about"
          className="mx-auto max-w-6xl px-5 py-20 sm:px-8 scroll-mt-20"
        >
          <header className="max-w-2xl mb-10">
            <p className="label-caps">Project</p>
            <h2 className="mt-4 text-3xl text-foreground sm:text-4xl">
              About DissolveAI
            </h2>
          </header>

          <div className="grid gap-14 lg:grid-cols-[minmax(0,1.15fr)_minmax(0,0.85fr)]">
            <div className="space-y-6 text-base leading-relaxed text-muted-foreground">
              <p>
                DissolveAI is a supervised machine-learning system for
                predicting the solubility of chemical compounds in a range of
                solvents at room temperature. It was built to support
                early-stage research, where the number of candidate
                compound–solvent combinations is far larger than the number that
                can reasonably be measured in the lab.
              </p>
              <p>
                The system pairs cheminformatics with regression modelling.
                Molecular structures given as SMILES are converted into
                physicochemical descriptors and Morgan fingerprints with RDKit;
                the solvent is one-hot encoded; a Random Forest regressor
                trained on cleaned BigSolDB measurements maps the combined
                features to a solubility value. When the predicted solubility in
                the requested solvent is poor, the service evaluates the same
                molecule in the remaining trained solvents and reports better
                options.
              </p>
              <p>
                The architecture is deliberately separated: an offline training
                pipeline produces the model and encoder artefacts, a Flask
                service loads them and exposes a single prediction endpoint, and
                this web interface is a thin client over that endpoint. No
                prediction logic lives in the browser, so what you see in the
                predictor is exactly what the model returns.
              </p>
            </div>

            <aside className="space-y-8">
              <div className="border-t border-rule-strong pt-5">
                <h3 className="label-caps">Technologies</h3>
                <ul className="mt-3 flex flex-wrap gap-2">
                  {[
                    "Python",
                    "RDKit",
                    "scikit-learn",
                    "Flask",
                    "NumPy",
                    "pandas",
                    "React",
                  ].map((t) => (
                    <li
                      key={t}
                      className="border border-border px-2.5 py-1 font-mono text-xs text-foreground"
                    >
                      {t}
                    </li>
                  ))}
                </ul>
              </div>

              <div className="border-t border-rule-strong pt-5">
                <h3 className="label-caps">Problem formulation</h3>
                <dl className="mt-3 space-y-3 text-sm">
                  <div>
                    <dt className="text-muted-foreground">Inputs</dt>
                    <dd className="text-foreground">
                      Molecular structure (SMILES), solvent identity
                    </dd>
                  </div>
                  <div>
                    <dt className="text-muted-foreground">Output</dt>
                    <dd className="text-foreground">
                      Continuous solubility value (log mol/L)
                    </dd>
                  </div>
                  <div>
                    <dt className="text-muted-foreground">Task</dt>
                    <dd className="text-foreground">Regression</dd>
                  </div>
                </dl>
              </div>

              <div className="border-t border-rule-strong pt-5">
                <h3 className="label-caps">Source</h3>
                <ul className="mt-3 space-y-2 text-sm">
                  <li>
                    <a
                      href={REPO_URL}
                      target="_blank"
                      rel="noreferrer noopener"
                      className="text-foreground underline decoration-border underline-offset-4 hover:text-primary"
                    >
                      Project repository
                    </a>
                  </li>
                  <li>
                    <a
                      href={PROFILE_URL}
                      target="_blank"
                      rel="noreferrer noopener"
                      className="text-foreground underline decoration-border underline-offset-4 hover:text-primary"
                    >
                      GitHub profile
                    </a>
                  </li>
                </ul>
              </div>
            </aside>
          </div>
        </section>
      </main>
      <SiteFooter />
    </div>
  );
}

function Stat({ value, label }: { value: string; label: string }) {
  return (
    <div>
      <dt className="label-caps">{label}</dt>
      <dd className="mt-1 font-mono text-2xl text-foreground">{value}</dd>
    </div>
  );
}

function Card({
  index,
  title,
  body,
}: {
  index: string;
  title: string;
  body: string;
}) {
  return (
    <article>
      <p className="font-mono text-xs text-primary">{index}</p>
      <h3 className="mt-3 text-xl text-foreground">{title}</h3>
      <p className="mt-3 text-sm leading-relaxed text-muted-foreground">
        {body}
      </p>
    </article>
  );
}

function Feature({ title, body }: { title: string; body: string }) {
  return (
    <article className="border-t border-rule-strong pt-4">
      <h3 className="text-lg text-foreground">{title}</h3>
      <p className="mt-2 text-sm leading-relaxed text-muted-foreground">
        {body}
      </p>
    </article>
  );
}

function Spec({ label, value }: { label: string; value: string }) {
  return (
    <div>
      <dt className="label-caps">{label}</dt>
      <dd className="mt-1 font-mono text-base text-foreground">{value}</dd>
    </div>
  );
}
