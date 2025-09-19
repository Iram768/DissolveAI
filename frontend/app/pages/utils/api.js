// utils/api.js
export async function predictSolubility(solute, solvent, temperature) {
    const res = await fetch("http://127.0.0.1:5000/predict", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        solute_smiles: solute,
        solvent_name: solvent,
        temperature: temperature ? Number(temperature) : null, // allow empty
      }),
    });
  
    if (!res.ok) {
      throw new Error("Prediction request failed");
    }
  
    return res.json();
  }
  