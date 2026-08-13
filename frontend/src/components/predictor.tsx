import { useState } from "react";
import { AnimatePresence, motion } from "motion/react";
import { MoleculeStructure } from "@/components/molecule-structure";
import {
  DEFAULT_BACKEND_URL,
  PredictionError,
  getBackendUrl,
  predictSolubility,
  setBackendUrl,
  type PredictionResponse,
} from "@/lib/predict-api";
import {
  EXAMPLE_MOLECULES,
  MODEL_METADATA,
  SOLVENT_GROUPS,
  SUPPORTED_SOLVENTS,
} from "@/lib/dissolveai";

type Status = "idle" | "running" | "done" | "error";

export function Predictor() {
  const [smiles, setSmiles] = useState("");
  const [solvent, setSolvent] = useState<string>("water");
  const [temperature, setTemperature] = useState("");
  const [status, setStatus] = useState<Status>("idle");
  const [result, setResult] = useState<PredictionResponse | null>(null);
  const [error, setError] = useState<{ title: string; message: string } | null>(
    null,
  );
  const [smilesError, setSmilesError] = useState<string | null>(null);
  const [showEndpoint, setShowEndpoint] = useState(false);
  const [endpoint, setEndpoint] = useState(() =>
    typeof window === "undefined" ? DEFAULT_BACKEND_URL : getBackendUrl(),
  );

  const tempNumber = temperature.trim() === "" ? null : Number(temperature);
  const tempOutOfRange =
    tempNumber !== null &&
    (Number.isNaN(tempNumber) ||
      tempNumber < MODEL_METADATA.temperatureMinK ||
      tempNumber > MODEL_METADATA.temperatureMaxK);

  async function run(event: React.FormEvent) {
    event.preventDefault();
    setError(null);

    if (!smiles.trim()) {
      setSmilesError("Enter a SMILES string to describe the molecule.");
      return;
    }
    setSmilesError(null);

    if (
      !SUPPORTED_SOLVENTS.includes(
        solvent as (typeof SUPPORTED_SOLVENTS)[number],
      )
    ) {
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
      setError({ title: pe.title, message: pe.message });
      setStatus("error");
    }
  }

  return (
    <div className="grid gap-10 lg:grid-cols-[minmax(0,1fr)_minmax(0,0.95fr)] lg:gap-14">
      {/* ---------------- experiment bench ---------------- */}
      <form onSubmit={run} className="space-y-9" noValidate>
        <fieldset className="space-y-3">
          <legend className="label-caps">Molecular structure</legend>
          <label
            htmlFor="smiles"
            className="block text-sm text-muted-foreground"
          >
            SMILES notation for the solute
          </label>
          <input
            id="smiles"
            name="smiles"
            value={smiles}
            onChange={(e) => {
              setSmiles(e.target.value);
              setSmilesError(null);
            }}
            spellCheck={false}
            autoComplete="off"
            placeholder="CCO"
            aria-invalid={smilesError ? true : undefined}
            aria-describedby={smilesError ? "smiles-error" : "smiles-help"}
            className="w-full border-b-2 border-input bg-transparent pb-2 font-mono text-2xl text-foreground transition-colors placeholder:text-muted-foreground/60 focus:border-primary focus:outline-none aria-[invalid=true]:border-destructive sm:text-3xl"
          />
          {smilesError ? (
            <p id="smiles-error" className="text-sm text-destructive">
              ✕ {smilesError}
            </p>
          ) : (
            <p
              id="smiles-help"
              className="flex flex-wrap items-center gap-2 text-sm"
            >
              <span className="text-muted-foreground">Example molecules:</span>
              {EXAMPLE_MOLECULES.map((m) => (
                <button
                  key={m.smiles}
                  type="button"
                  onClick={() => setSmiles(m.smiles)}
                  className="border border-border px-2 py-0.5 font-mono text-xs text-foreground transition-colors hover:border-primary hover:text-primary"
                >
                  {m.smiles}
                  <span className="ml-1.5 font-sans text-[0.65rem] text-muted-foreground">
                    {m.label}
                  </span>
                </button>
              ))}
            </p>
          )}
        </fieldset>

        <fieldset className="space-y-4">
          <legend className="label-caps">
            Solvent · {MODEL_METADATA.totalSolvents} supported
          </legend>
          <div className="max-h-72 space-y-5 overflow-y-auto border border-border bg-card p-4">
            {SOLVENT_GROUPS.map((group) => (
              <div key={group.label}>
                <p className="font-mono text-[0.65rem] uppercase tracking-[0.14em] text-muted-foreground">
                  {group.label}
                </p>
                <div className="mt-2 flex flex-wrap gap-1.5">
                  {group.solvents.map((s) => {
                    const active = s === solvent;
                    return (
                      <button
                        key={s}
                        type="button"
                        onClick={() => setSolvent(s)}
                        aria-pressed={active}
                        className={
                          active
                            ? "border border-primary bg-primary px-2.5 py-1 text-xs text-primary-foreground"
                            : "border border-border px-2.5 py-1 text-xs text-foreground transition-colors hover:border-primary/60 hover:text-primary"
                        }
                      >
                        {active && <span aria-hidden="true">✓ </span>}
                        {s}
                      </button>
                    );
                  })}
                </div>
              </div>
            ))}
          </div>
        </fieldset>

        <fieldset className="space-y-2">
          <legend className="label-caps">Temperature — optional</legend>
          <div className="flex flex-wrap items-baseline gap-3">
            <input
              id="temperature"
              type="number"
              inputMode="decimal"
              step="0.01"
              value={temperature}
              onChange={(e) => setTemperature(e.target.value)}
              placeholder={String(MODEL_METADATA.defaultTemperatureK)}
              aria-describedby="temp-help"
              className="w-36 border-b border-input bg-transparent pb-1 font-mono text-lg focus:border-primary focus:outline-none"
            />
            <span className="font-mono text-sm text-muted-foreground">K</span>
          </div>
          <p
            id="temp-help"
            className={
              tempOutOfRange
                ? "text-sm text-accent-foreground"
                : "text-sm text-muted-foreground"
            }
          >
            Training data spans {MODEL_METADATA.temperatureMinK}–
            {MODEL_METADATA.temperatureMaxK} K. Left blank, the backend uses{" "}
            {MODEL_METADATA.defaultTemperatureK} K; values outside the range are
            clamped by the backend.
          </p>
        </fieldset>

        <div className="flex flex-wrap items-center gap-4 border-t border-border pt-6">
          <button
            type="submit"
            disabled={status === "running"}
            className="border border-primary bg-primary px-6 py-3 font-mono text-xs uppercase tracking-[0.18em] text-primary-foreground transition-opacity hover:opacity-90 disabled:opacity-60"
          >
            {status === "running" ? "Running…" : "Run prediction"}
          </button>
          <button
            type="button"
            onClick={() => setShowEndpoint((v) => !v)}
            aria-expanded={showEndpoint}
            className="font-mono text-xs uppercase tracking-[0.14em] text-muted-foreground underline decoration-border underline-offset-4 hover:text-primary"
          >
            Prediction service
          </button>
        </div>

        {showEndpoint && (
          <div className="border border-border bg-secondary/50 p-4">
            <label htmlFor="endpoint" className="label-caps">
              Flask backend base URL
            </label>
            <input
              id="endpoint"
              value={endpoint}
              onChange={(e) => setEndpoint(e.target.value)}
              onBlur={() => setBackendUrl(endpoint)}
              spellCheck={false}
              className="mt-2 w-full border-b border-input bg-transparent pb-1 font-mono text-sm focus:border-primary focus:outline-none"
            />
            <p className="mt-3 text-sm leading-relaxed text-muted-foreground">
              Predictions are computed by this project&apos;s own Flask service
              (
              <code className="font-mono">
                solvent_project/model_training_BigSolDB/app.py
              </code>
              ), which loads the trained Random Forest model and solvent
              encoder. Run it locally and this interface posts to{" "}
              <code className="font-mono">{endpoint}/predict</code>. No values
              are generated in the browser.
            </p>
          </div>
        )}
      </form>

      {/* ---------------- structure + output ---------------- */}
      <div className="space-y-6">
        <div>
          <p className="label-caps mb-3">
            Structure preview — rendered from your SMILES
          </p>
          <div className="relative">
            <MoleculeStructure smiles={smiles} />
            {status === "running" && (
              <div className="pointer-events-none absolute inset-0 overflow-hidden">
                <div className="animate-scan h-px w-full bg-primary/70" />
              </div>
            )}
          </div>
        </div>

        <AnimatePresence mode="wait">
          {status === "running" && (
            <motion.div
              key="running"
              initial={{ opacity: 0, y: 8 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0 }}
              className="paper-panel p-6"
            >
              <p className="label-caps">Request in flight</p>
              <ul className="mt-4 space-y-2 text-sm text-muted-foreground">
                <li className="flex items-center gap-3">
                  <Dot /> RDKit parses the structure and computes descriptors
                </li>
                <li className="flex items-center gap-3">
                  <Dot /> Solvent one-hot encoding is appended to the feature
                  row
                </li>
                <li className="flex items-center gap-3">
                  <Dot /> Random Forest regressor estimates solubility
                </li>
              </ul>
              <p className="mt-4 font-mono text-xs text-muted-foreground">
                Awaiting the backend response — stages are the pipeline in
                app.py, not a progress simulation.
              </p>
            </motion.div>
          )}

          {status === "error" && error && (
            <motion.div
              key="error"
              initial={{ opacity: 0, y: 8 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0 }}
              role="alert"
              className="border border-destructive/40 bg-destructive/5 p-6"
            >
              <p className="font-display text-xl text-foreground">
                {error.title}
              </p>
              <p className="mt-2 text-sm leading-relaxed text-muted-foreground">
                {error.message}
              </p>
            </motion.div>
          )}

          {status === "done" && result && (
            <motion.div
              key="result"
              initial={{ opacity: 0, y: 10 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.4, ease: [0.22, 1, 0.36, 1] }}
              className="space-y-6"
            >
              <ResultPanel result={result} />
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </div>
  );
}

function Dot() {
  return (
    <span
      aria-hidden="true"
      className="h-1.5 w-1.5 shrink-0 rounded-full bg-primary"
    />
  );
}

function ResultPanel({ result }: { result: PredictionResponse }) {
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

  return (
    <>
      <section className="paper-panel p-6 sm:p-8" aria-live="polite">
        <p className="label-caps">Predicted solubility</p>
        <p className="mt-2 font-mono text-5xl tracking-tight text-foreground sm:text-6xl">
          {value.toFixed(3)}
        </p>
        <p className="mt-1 font-mono text-sm text-muted-foreground">
          log mol/L
        </p>

        <div className="mt-8">
          <div className="relative h-8 border border-border bg-secondary/60">
            <div
              className="absolute inset-y-0 left-0 bg-primary/20"
              style={{ width: `${position}%` }}
            />
            <div
              className="absolute inset-y-0 w-px bg-primary"
              style={{ left: `${position}%` }}
              aria-hidden="true"
            />
            <div
              className="absolute inset-y-0 w-px border-l border-dashed border-accent-foreground/60"
              style={{ left: `${thresholdPos}%` }}
              aria-hidden="true"
            />
          </div>
          <div className="mt-2 flex justify-between font-mono text-[0.65rem] text-muted-foreground">
            <span>{lower.toFixed(1)}</span>
            <span>
              threshold {MODEL_METADATA.solubilityThreshold.toFixed(1)}
            </span>
            <span>{upper.toFixed(1)}</span>
          </div>
          <p className="mt-2 text-xs text-muted-foreground">
            Output range is clipped by the backend to [{lower}, {upper}] log
            mol/L; the dashed mark is the repository&apos;s poor-solubility
            threshold.
          </p>
        </div>

        <dl className="mt-8 grid gap-5 border-t border-border pt-6 sm:grid-cols-2">
          <div>
            <dt className="label-caps">Compound</dt>
            <dd className="mt-1 break-words text-sm text-foreground">
              {result.solute_name}
            </dd>
          </div>
          <div>
            <dt className="label-caps">Solvent</dt>
            <dd className="mt-1 text-sm text-foreground">
              {result.solvent_name}
            </dd>
          </div>
          <div>
            <dt className="label-caps">Temperature used</dt>
            <dd className="mt-1 font-mono text-sm text-foreground">
              {result.temperature_used_K} K
            </dd>
          </div>
          <div>
            <dt className="label-caps">Model</dt>
            <dd className="mt-1 text-sm text-foreground">
              Random Forest Regressor
            </dd>
          </div>
        </dl>

        {result.note && (
          <p className="mt-6 border-l-2 border-accent-foreground/50 bg-accent/40 px-4 py-3 text-sm text-accent-foreground">
            {result.note}
          </p>
        )}
        {result.message && (
          <p className="mt-4 text-sm leading-relaxed text-muted-foreground">
            {result.message}
          </p>
        )}
      </section>

      {suggestions.length > 0 && (
        <section className="paper-panel p-6 sm:p-8">
          <h3 className="font-display text-2xl text-foreground">
            Looking for a better solvent?
          </h3>
          <p className="mt-2 text-sm leading-relaxed text-muted-foreground">
            The backend re-ran the same molecule through every remaining solvent
            in the encoder and returned those that clear its solubility
            threshold.
          </p>
          <ul className="mt-6 space-y-3">
            {suggestions.map((s) => {
              const width = Math.min(
                100,
                Math.max(
                  2,
                  ((s.predicted_solubility - lower) / (upper - lower)) * 100,
                ),
              );
              return (
                <li key={s.solvent}>
                  <div className="flex items-baseline justify-between gap-4">
                    <span className="text-sm text-foreground">{s.solvent}</span>
                    <span className="font-mono text-sm text-foreground">
                      {s.predicted_solubility.toFixed(3)}
                      <span className="ml-1 text-xs text-muted-foreground">
                        log mol/L
                      </span>
                    </span>
                  </div>
                  <div className="mt-1.5 h-1.5 bg-secondary">
                    <div
                      className="h-full bg-primary/70"
                      style={{ width: `${width}%` }}
                    />
                  </div>
                </li>
              );
            })}
          </ul>
        </section>
      )}
    </>
  );
}
