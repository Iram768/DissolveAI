// frontend/frontend/pages/index.js
import Link from "next/link";

export default function Home() {
  return (
    <div className="min-h-screen flex flex-col bg-gray-50">
      {/* Header */}
      <header className="bg-blue-600 text-white p-4 flex justify-between items-center">
        <h1 className="text-xl font-bold">DissolveAI</h1>
        <nav className="space-x-4">
          <Link href="/" className="hover:underline">About</Link>
          <Link href="/project" className="hover:underline">Check Solubility Prediction</Link>
        </nav>
      </header>

      {/* Main content */}
      <main className="flex flex-1 items-center justify-between p-8 select-none">
        <div className="max-w-lg">
          <h2 className="text-2xl font-semibold mb-4">About the Project</h2>
          <p>
            DissolveAI is a machine learning tool that predicts solubility of solutes
            in different solvents using SMILES representations and experimental data.
            This project combines chemistry and artificial intelligence.
          </p>
        </div>
        <div>
          <img src="/logo.png" alt="Project Logo" className="w-48 h-48" />
        </div>
      </main>

      {/* Footer */}
      <footer className="bg-gray-800 text-white p-4 text-center">
        <p>Contact: yourgmail@gmail.com | LinkedIn: linkedin.com/in/yourprofile</p>
      </footer>
    </div>
  );
}



"use client";
import { useEffect, useState } from "react";
import CountUp from "react-countup";

export default function Home() {
  const [visits, setVisits] = useState(0);

  useEffect(() => {
    // Call your Flask backend
    fetch("http://127.0.0.1:5000/visit-count")
      .then((res) => res.json())
      .then((data) => setVisits(data.visits))
      .catch((err) => console.error("Error fetching visits:", err));
  }, []);

  return (
    <div className="flex flex-col items-center justify-center min-h-screen bg-gray-100">
      <h1 className="text-3xl font-bold mb-4">Welcome to DissolveAI</h1>

      <p className="text-lg mb-6">Solubility prediction made easy.</p>

      <div className="text-2xl font-semibold text-blue-600">
        Visitors:{" "}
        <CountUp
          end={visits}
          duration={1.5}
          separator=","
        />
      </div>
    </div>
  );
}
