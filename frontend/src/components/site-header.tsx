import { Link } from "@tanstack/react-router";
import { useEffect, useState } from "react";
import { REPO_URL } from "@/lib/dissolveai";

const NAV = [
  { to: "#home", label: "Home" },
  { to: "#predictor", label: "Predictor" },
  { to: "#how-it-works", label: "How It Works" },
  { to: "#model", label: "Model" },
  { to: "#about", label: "About" },
] as const;

export function SiteHeader() {
  const [scrolled, setScrolled] = useState(false);
  const [open, setOpen] = useState(false);

  useEffect(() => {
    const onScroll = () => setScrolled(window.scrollY > 12);
    onScroll();
    window.addEventListener("scroll", onScroll, { passive: true });
    return () => window.removeEventListener("scroll", onScroll);
  }, []);

  return (
    <header
      className={`sticky top-0 z-50 transition-colors duration-300 ${
        scrolled
          ? "border-b border-border bg-background/85 backdrop-blur-md"
          : "border-b border-transparent"
      }`}
    >
      <div className="mx-auto flex max-w-6xl items-center gap-6 px-5 py-4 sm:px-8">
        <Link
          to="/"
          hash="home"
          className="group flex items-center gap-3"
          aria-label="DissolveAI home"
        >
          <FlaskMark />
          <span className="leading-tight">
            <span className="block font-display text-lg tracking-tight text-foreground">
              DissolveAI
            </span>
            <span className="label-caps hidden sm:block">
              Molecular Solubility Prediction
            </span>
          </span>
        </Link>

        <nav
          className="ml-auto hidden items-center gap-7 md:flex"
          aria-label="Main"
        >
          {NAV.map((item) => (
            <Link
              key={item.to}
              to="/"
              hash={item.to.replace("#", "")}
              className="text-sm text-muted-foreground transition-colors hover:text-foreground"
            >
              {item.label}
            </Link>
          ))}
          <a
            href={REPO_URL}
            target="_blank"
            rel="noreferrer noopener"
            className="border border-input px-3 py-1.5 font-mono text-xs uppercase tracking-[0.14em] text-foreground transition-colors hover:border-primary hover:text-primary"
          >
            GitHub
          </a>
        </nav>

        <button
          type="button"
          onClick={() => setOpen((v) => !v)}
          aria-expanded={open}
          aria-controls="mobile-nav"
          className="ml-auto flex h-9 w-9 flex-col items-center justify-center gap-1.5 border border-input md:hidden"
        >
          <span className="sr-only">{open ? "Close menu" : "Open menu"}</span>
          <span className="block h-px w-4 bg-foreground" />
          <span className="block h-px w-4 bg-foreground" />
        </button>
      </div>

      {open && (
        <nav
          id="mobile-nav"
          aria-label="Main, mobile"
          className="border-t border-border bg-background md:hidden"
        >
          <ul className="mx-auto max-w-6xl px-5 py-2 sm:px-8">
            {NAV.map((item) => (
              <li
                key={item.to}
                className="border-b border-border last:border-b-0"
              >
                <Link
                  to="/"
                  hash={item.to.replace("#", "")}
                  onClick={() => setOpen(false)}
                  className="block py-3 text-sm text-foreground"
                >
                  {item.label}
                </Link>
              </li>
            ))}
            <li className="py-3">
              <a
                href={REPO_URL}
                target="_blank"
                rel="noreferrer noopener"
                className="font-mono text-xs uppercase tracking-[0.14em] text-primary"
              >
                GitHub repository
              </a>
            </li>
          </ul>
        </nav>
      )}
    </header>
  );
}

function FlaskMark() {
  return (
    <svg
      width="26"
      height="26"
      viewBox="0 0 24 24"
      aria-hidden="true"
      className="text-primary"
      fill="none"
      stroke="currentColor"
      strokeWidth="1.1"
    >
      <path d="M9 3h6M10 3v5.2L4.9 18.4A2 2 0 0 0 6.7 21.4h10.6a2 2 0 0 0 1.8-3L14 8.2V3" />
      <path d="M7.4 15.4h9.2" />
      <circle cx="10.6" cy="17.9" r="0.9" fill="currentColor" stroke="none" />
      <circle cx="13.9" cy="18.6" r="0.6" fill="currentColor" stroke="none" />
    </svg>
  );
}
