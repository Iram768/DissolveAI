import { Link } from "@tanstack/react-router";
import { PROFILE_URL, REPO_URL } from "@/lib/dissolveai";

const NAV = [
  { to: "#home", label: "Home" },
  { to: "#predictor", label: "Predictor" },
  { to: "#how-it-works", label: "How It Works" },
  { to: "#model", label: "Model" },
  { to: "#about", label: "About" },
] as const;

const CONTACT = [
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
  { label: "GitHub", href: PROFILE_URL, text: "github.com/Iram768" },
  {
    label: "Repository",
    href: REPO_URL,
    text: "github.com/Iram768/dissolveai",
  },
];

export function SiteFooter() {
  return (
    <footer className="mt-24 border-t border-border bg-secondary/60">
      <div className="mx-auto grid max-w-6xl gap-10 px-5 py-14 sm:px-8 md:grid-cols-3">
        <div>
          <p className="font-display text-xl text-foreground">DissolveAI</p>
          <p className="label-caps mt-1">Molecular Solubility Prediction</p>
          <p className="mt-4 max-w-xs text-sm leading-relaxed text-muted-foreground">
            A supervised machine-learning system for estimating compound
            solubility across solvents using RDKit-derived features and a Random
            Forest regressor.
          </p>
        </div>

        <nav aria-label="Footer">
          <h2 className="label-caps">Navigation</h2>
          <ul className="mt-4 space-y-2">
            {NAV.map((item) => (
              <li key={item.to}>
                <Link
                  to="/"
                  hash={item.to.replace("#", "")}
                  className="text-sm text-muted-foreground transition-colors hover:text-primary"
                >
                  {item.label}
                </Link>
              </li>
            ))}
          </ul>
        </nav>

        <div>
          <h2 className="label-caps">Contact</h2>
          <ul className="mt-4 space-y-2">
            {CONTACT.map((item) => (
              <li key={item.label} className="text-sm">
                <span className="mr-2 font-mono text-xs uppercase tracking-[0.12em] text-muted-foreground">
                  {item.label}
                </span>
                <a
                  href={item.href}
                  target={
                    item.href.startsWith("mailto:") ? undefined : "_blank"
                  }
                  rel="noreferrer noopener"
                  className="break-all text-foreground underline decoration-border underline-offset-4 transition-colors hover:text-primary hover:decoration-primary"
                >
                  {item.text}
                </a>
              </li>
            ))}
          </ul>
        </div>
      </div>
      <div className="border-t border-border">
        <p className="mx-auto max-w-6xl px-5 py-5 font-mono text-xs text-muted-foreground sm:px-8">
          DissolveAI — research prototype. Predictions are computational
          estimates and are not a substitute for experimental measurement.
        </p>
      </div>
    </footer>
  );
}
