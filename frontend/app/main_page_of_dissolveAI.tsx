"use client";

import { useState } from "react";

export default function Home() {
  const [smiles, setSmiles] = useState("");
  const [solvent, setSolvent] = useState("");
  const [temperature, setTemperature] = useState("");
  const [result, setResult] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setResult(null);

    try {
      const res = await fetch("http://127.0.0.1:5000/predict", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          solute_smiles: smiles,
          solvent_name: solvent,
          temperature: temperature ? parseFloat(temperature) : null,
        }),
      });

      const data = await res.json();
      setResult(data.prediction || data.message);
    } catch (err) {
      setResult("Error connecting to backend.");
    }
    setLoading(false);
  };

  return (
    <div className="min-h-screen flex flex-col items-center justify-between bg-gray-50 text-gray-800">
      {/* Header */}
      <header className="w-full p-6 bg-indigo-600 text-white text-center text-2xl font-bold shadow-md">
        DissolveAI
      </header>

      {/* Hero Image */}
      <div className="w-full flex justify-center mt-8">
        <img
          src="/chemistry.jpg"
          alt="Chemistry AI"
          className="w-2/3 rounded-xl shadow-lg"
        />
      </div>

      {/* About Section */}
      <section className="max-w-3xl text-center my-10 px-6">
        <h2 className="text-2xl font-bold mb-4">Transforming Chemistry with AI</h2>
        <p className="leading-relaxed">
          Chemistry is full of experiments, but not every answer needs hours in
          the lab. DissolveAI makes it easier for students and researchers to
          understand which solvents work best for a compound before they even
          start mixing.
        </p>
        <p className="leading-relaxed mt-4">
          This tool uses machine learning to study patterns in how compounds
          dissolve in different solvents. We trained our model with a dataset
          that includes thousands of examples of solubility behavior, giving it
          the ability to recognize hidden rules that are hard for humans to see.
        </p>
        <p className="leading-relaxed mt-4">
          With DissolveAI, you don’t just get a yes or no answer. You get a
          prediction that guides your choices, saving time, reducing waste, and
          building confidence in your results.
        </p>
      </section>

      {/* Tools Section */}
      <section className="w-full max-w-xl p-6 bg-white rounded-xl shadow-lg mb-10">
        <h2 className="text-xl font-semibold mb-4 text-center">Solubility Prediction Tool</h2>
        <form onSubmit={handleSubmit} className="space-y-4">
          <input
            type="text"
            placeholder="Enter solute SMILES"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            className="w-full border p-2 rounded"
            required
          />
          <input
            type="text"
            placeholder="Enter solvent name"
            value={solvent}
            onChange={(e) => setSolvent(e.target.value)}
            className="w-full border p-2 rounded"
            required
          />
          <input
            type="number"
            placeholder="Enter temperature (K, optional)"
            value={temperature}
            onChange={(e) => setTemperature(e.target.value)}
            className="w-full border p-2 rounded"
          />
          <button
            type="submit"
            className="w-full bg-indigo-600 text-white py-2 rounded hover:bg-indigo-700 transition"
          >
            {loading ? "Predicting..." : "Predict Solubility"}
          </button>
        </form>

        {result && (
          <div className="mt-4 p-3 bg-gray-100 rounded text-center">
            <strong>Result:</strong> {result}
          </div>
        )}
      </section>

      {/* Footer */}
      <footer className="w-full bg-gray-800 text-white text-center py-4">
        Contact: your.email@gmail.com | LinkedIn: /in/yourprofile
      </footer>
    </div>
  );
}
