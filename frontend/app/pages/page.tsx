"use client";

import { useState } from "react";

export default function Project() {
  const [smiles, setSmiles] = useState("");
  const [solvent, setSolvent] = useState("");
  const [temperature, setTemperature] = useState("");
  const [result, setResult] = useState<string | null>(null);

  const handleSubmit = async () => {
    const response = await fetch("http://127.0.0.1:5000/predict", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ smiles, solvent, temperature }),
    });
    const data = await response.json();
    setResult(data.prediction);
  };

  return (
    <main className="flex flex-col min-h-screen items-center p-8 bg-gray-50 text-gray-900">
      <h1 className="text-3xl font-bold mb-6">Solubility Prediction</h1>

      <div className="flex flex-col space-y-4 w-full max-w-md">
        <input
          type="text"
          placeholder="Enter Solute SMILES"
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          className="p-2 border rounded-lg"
        />
        <input
          type="text"
          placeholder="Enter Solvent"
          value={solvent}
          onChange={(e) => setSolvent(e.target.value)}
          className="p-2 border rounded-lg"
        />
        <input
          type="number"
          placeholder="Enter Temperature (K, optional)"
          value={temperature}
          onChange={(e) => setTemperature(e.target.value)}
          className="p-2 border rounded-lg"
        />
        <button
          onClick={handleSubmit}
          className="bg-blue-500 text-white py-2 rounded-lg hover:bg-blue-600"
        >
          Predict
        </button>
      </div>

      {result && (
        <p className="mt-6 text-xl font-semibold">Prediction: {result}</p>
      )}
    </main>
  );
}
