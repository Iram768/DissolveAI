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

const DEFAULT_BASE_URL = "http://127.0.0.1:5000";
const STORAGE_KEY = "dissolveai.backend-url";

export function getBackendUrl(): string {
  if (typeof window !== "undefined") {
    const stored = window.localStorage.getItem(STORAGE_KEY);
    if (stored && stored.trim()) return stored.trim().replace(/\/+$/, "");
  }
  const fromEnv = import.meta.env["VITE_DISSOLVEAI_API_URL"] as
    string | undefined;
  return (fromEnv?.trim() || DEFAULT_BASE_URL).replace(/\/+$/, "");
}

export function setBackendUrl(url: string) {
  if (typeof window === "undefined") return;
  if (url.trim()) window.localStorage.setItem(STORAGE_KEY, url.trim());
  else window.localStorage.removeItem(STORAGE_KEY);
}

export const DEFAULT_BACKEND_URL = DEFAULT_BASE_URL;

export type SuggestedSolvent = {
  solvent: string;
  predicted_solubility: number;
};

export type PredictionResponse = {
  solute_name?: string;
  solvent_name?: string;
  temperature_used_K?: number;
  predicted_solubility?: number;
  message?: string;
  note?: string;
  suggested_solvents?: SuggestedSolvent[];
};

export type PredictionRequest = {
  smiles: string;
  solvent: string;
  temperatureK?: number | null;
};

export class PredictionError extends Error {
  title: string;
  detail?: string | undefined;

  constructor(title: string, message: string, detail?: string) {
    super(message);
    this.name = "PredictionError";
    this.title = title;
    this.detail = detail;
  }
}

export async function predictSolubility(
  input: PredictionRequest,
): Promise<PredictionResponse> {
  const url = `${getBackendUrl()}/predict`;

  let response: Response;
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

  let payload: PredictionResponse & { error?: string; detail?: string };
  try {
    payload = (await response.json()) as typeof payload;
  } catch {
    throw new PredictionError(
      "Unreadable response",
      "The backend replied with a payload that is not valid JSON.",
    );
  }

  if (!response.ok || payload.error) {
    const error =
      payload.error ?? `Request failed with status ${response.status}`;
    if (/invalid smiles/i.test(error)) {
      throw new PredictionError(
        "Unable to analyze molecule",
        "The provided SMILES could not be parsed by RDKit. Please check the molecular structure and try again.",
      );
    }
    if (/no solute smiles/i.test(error)) {
      throw new PredictionError(
        "No molecular structure supplied",
        "Enter a SMILES string before running a prediction.",
      );
    }
    throw new PredictionError("Prediction failed", error, payload.detail);
  }

  // Unknown-solvent case: backend returns 200 with only a message.
  if (payload.predicted_solubility === undefined) {
    throw new PredictionError(
      "Solvent not covered by the model",
      payload.message ??
        "The selected solvent is not present in the trained solvent encoder.",
    );
  }

  return payload;
}

/** Presentation helper only — mirrors the backend's own metadata threshold. */
export function solubilityBandLabel(value: number, threshold: number): string {
  return value > threshold
    ? "Above dataset threshold"
    : "Below dataset threshold";
}
