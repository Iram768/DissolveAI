"use client";

import { useEffect, useState } from "react";
import CountUp from "react-countup";
import Image from "next/image";

export default function About() {
  const [visits, setVisits] = useState<number | null>(null);

  useEffect(() => {
    fetch("http://127.0.0.1:5000/visit-count")
      .then((res) => res.json())
      .then((data) => setVisits(data.visits))
      .catch(() => setVisits(null));
  }, []);

  return (
    <main className="flex flex-col min-h-screen bg-gray-50 text-gray-900 select-none">
      {/* Header */}
      <header className="w-full flex justify-between items-center p-4 shadow-md bg-white">
        <h1 className="text-2xl font-bold">DissolveAI</h1>
        <nav className="space-x-4">
          <a href="/about" className="hover:underline">About</a>
          <a href="/project" className="hover:underline">Check Solubility Prediction</a>
        </nav>
      </header>

      {/* Main Content */}
      <section className="flex flex-1 justify-between items-center p-12">
        {/* Left Side */}
        <div className="max-w-lg">
          <h2 className="text-3xl font-semibold mb-4">About the Project</h2>
          <p className="text-lg leading-relaxed">
            DissolveAI predicts solubility of solutes in different solvents using
            machine learning and cheminformatics. Enter solute SMILES and solvent
            to get predictions.
          </p>
          {visits !== null && (
            <p className="mt-6 text-xl">
              Visitors: <CountUp end={visits} duration={1.5} />
            </p>
          )}
        </div>

        {/* Right Side (Logo) */}
        <div>
          <Image
            src="/DissolveAIlogo.png"
            alt="DissolveAI Logo"
            width={200}
            height={200}
          />
        </div>
      </section>

      {/* Footer */}
      <footer className="w-full text-center p-4 border-t border-gray-300">
        Contact:{" "}
        <a href="mailto:iramjaved751@gmail.com" className="underline">iramjaved751@gmail.com</a>{" "}
        |{" "}
        <a href="https://linkedin.com/in/iram-javed" target="_blank" className="underline">
          LinkedIn
        </a>
      </footer>
    </main>
  );
}

