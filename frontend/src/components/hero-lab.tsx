/**
 * Hero visual: a restrained scientific illustration — a slowly rotating
 * molecular skeleton inside measurement rules, with a few drifting solvent
 * particles. Pure SVG/CSS, no WebGL, ~20 animated nodes.
 */
export function HeroLab() {
  return (
    <div className="relative">
      <div className="rule-grid border border-border bg-card p-4 sm:p-6">
        <div className="flex items-center justify-between border-b border-border pb-3">
          <p className="label-caps">Fig. 1 — solute in solvent field</p>
          <p className="font-mono text-[0.65rem] text-muted-foreground">
            298 K
          </p>
        </div>

        <svg
          viewBox="0 0 320 300"
          role="img"
          aria-label="Illustration of a molecular skeleton surrounded by solvent particles and measurement guides"
          className="mt-4 h-auto w-full"
        >
          {/* measurement guides */}
          <g stroke="var(--rule-strong)" strokeWidth="1" fill="none">
            <path d="M24 34 H296" strokeDasharray="2 6" />
            <path d="M24 266 H296" strokeDasharray="2 6" />
            <path d="M24 34 V44 M296 34 V44" />
            <path d="M24 256 V266 M296 256 V266" />
          </g>
          <text
            x="160"
            y="26"
            textAnchor="middle"
            className="fill-muted-foreground"
            style={{ font: "9px var(--font-mono)", letterSpacing: "0.14em" }}
          >
            SOLVENT SHELL
          </text>
          <text
            x="160"
            y="286"
            textAnchor="middle"
            className="fill-muted-foreground"
            style={{ font: "9px var(--font-mono)", letterSpacing: "0.14em" }}
          >
            log mol/L ESTIMATE
          </text>

          {/* solvent particles */}
          <g
            className="animate-drift"
            style={{ transformOrigin: "160px 150px" }}
          >
            {[
              [52, 78],
              [268, 96],
              [70, 214],
              [258, 210],
              [160, 58],
              [104, 128],
              [220, 176],
              [186, 240],
            ].map(([cx, cy], i) => (
              <g key={`${cx}-${cy}`}>
                <circle
                  cx={cx}
                  cy={cy}
                  r={i % 3 === 0 ? 4 : 2.6}
                  fill="var(--teal-soft)"
                />
                <circle
                  cx={cx}
                  cy={cy}
                  r={i % 3 === 0 ? 4 : 2.6}
                  fill="none"
                  stroke="var(--teal)"
                  strokeWidth="0.6"
                />
              </g>
            ))}
          </g>

          {/* rotating molecular skeleton */}
          <g
            style={{
              transformOrigin: "160px 150px",
              animation: "dissolve-spin 64s linear infinite",
            }}
          >
            <g stroke="var(--ink-soft)" strokeWidth="1.4" fill="none">
              <path d="M110 150 L136 106 L188 106 L214 150 L188 194 L136 194 Z" />
              <path d="M140 112 L184 112 M140 188 L184 188" strokeWidth="1" />
              <path d="M214 150 L252 150" />
              <path d="M252 150 L272 118" />
              <path d="M110 150 L78 150" />
              <path d="M136 106 L120 70" />
            </g>
            {[
              [110, 150],
              [136, 106],
              [188, 106],
              [214, 150],
              [188, 194],
              [136, 194],
              [252, 150],
            ].map(([cx, cy]) => (
              <circle
                key={`${cx}-${cy}`}
                cx={cx}
                cy={cy}
                r="3.4"
                fill="var(--ink)"
              />
            ))}
            <circle
              cx="272"
              cy="118"
              r="5"
              fill="var(--paper)"
              stroke="var(--teal)"
              strokeWidth="1.4"
            />
            <text
              x="280"
              y="112"
              className="fill-primary"
              style={{ font: "9px var(--font-mono)" }}
            >
              OH
            </text>
            <circle
              cx="78"
              cy="150"
              r="5"
              fill="var(--paper)"
              stroke="var(--ochre)"
              strokeWidth="1.4"
            />
            <circle
              cx="120"
              cy="70"
              r="4"
              fill="var(--paper)"
              stroke="var(--ink-soft)"
              strokeWidth="1.2"
            />
          </g>
        </svg>

        <div className="mt-4 grid grid-cols-3 gap-3 border-t border-border pt-3">
          {[
            ["RDKit", "descriptors"],
            ["Morgan", "radius 2 · 1024 bit"],
            ["Random Forest", "500 trees"],
          ].map(([k, v]) => (
            <div key={k}>
              <p className="font-mono text-[0.68rem] text-foreground">{k}</p>
              <p className="font-mono text-[0.6rem] text-muted-foreground">
                {v}
              </p>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}
