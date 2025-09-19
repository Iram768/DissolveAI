// frontend/frontend/pages/project.js
import { useState } from "react";

export default function Project() {
  const [smiles, setSmiles] = useState("");
  const [solvent, setSolvent] = useState("");
  const [temperature, setTemperature] = useState("");
  const [result, setResult] = useState(null);

  const handleSubmit = async (e) => {
    e.preventDefault();
    const response = await fetch("http://127.0.0.1:5000/predict", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ solute: smiles, solvent, temperature }),
    });
    const data = await response.json();
    setResult(data);
  };

  return (
    <div className="min-h-screen flex flex-col bg-gray-50">
      {/* Header */}
      <header className="bg-blue-600 text-white p-4 flex justify-between items-center">
        <h1 className="text-xl font-bold">DissolveAI - Prediction</h1>
      </header>

      {/* Form */}
      <main className="flex-1 p-8">
        <form onSubmit={handleSubmit} className="space-y-4 max-w-md mx-auto bg-white shadow-md p-6 rounded-lg">
          <div>
            <label className="block mb-1">Solute SMILES:</label>
            <input
              type="text"
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              className="w-full border p-2 rounded"
              required
            />
          </div>
          <div>
            <label className="block mb-1">Solvent:</label>
            <input
              type="text"
              value={solvent}
              onChange={(e) => setSolvent(e.target.value)}
              className="w-full border p-2 rounded"
              required
            />
          </div>
          <div>
            <label className="block mb-1">Temperature (K, optional):</label>
            <input
              type="number"
              value={temperature}
              onChange={(e) => setTemperature(e.target.value)}
              className="w-full border p-2 rounded"
            />
          </div>
          <button
            type="submit"
            className="bg-blue-600 text-white px-4 py-2 rounded hover:bg-blue-700"
          >
            Predict
          </button>
        </form>
        {result && (
            <div className="mt-6 max-w-md mx-auto bg-green-100 p-4 rounded">
            <p>
            {result.user_provided_temp
                ? `Predicted solubility: ${result.predicted_solubility}`
                : `The solute ${result.solute_name} is soluble in ${result.solvent} at room temperature (298K).`}
            </p>
            </div>
        )}

{/* this section is replaced with above just to get the name of solutes insted of smiles in result  */}

        {/* {result && (
          <div className="mt-6 max-w-md mx-auto bg-green-100 p-4 rounded">
            <p>
              {temperature
                ? `Predicted solubility: ${result.predicted_solubility}`
                : `The solute ${result.solute_name} is soluble in ${solvent} at room temperature (298K).`}
            </p>
          </div> 
        )} */}
      </main>

      {/* Footer */}
      <footer className="bg-gray-800 text-white p-4 text-center">
        <p>Contact: yourgmail@gmail.com | LinkedIn: linkedin.com/in/yourprofile</p>
      </footer>
    </div>
  );
}
