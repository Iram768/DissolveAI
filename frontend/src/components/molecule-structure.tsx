import { useEffect, useRef, useState } from "react";

type Props = {
  smiles: string;
  height?: number;
  className?: string;
};

/**
 * 2D structure rendered client-side from the SMILES the user typed, using
 * smiles-drawer. Nothing decorative: if the string cannot be parsed we say so.
 */
export function MoleculeStructure({ smiles, height = 260, className }: Props) {
  const svgRef = useRef<SVGSVGElement | null>(null);
  const [parseError, setParseError] = useState(false);

  useEffect(() => {
    let cancelled = false;
    const svg = svgRef.current;
    if (!svg) return;

    if (!smiles.trim()) {
      svg.innerHTML = "";
      setParseError(false);
      return;
    }

    void (async () => {
      const mod = await import("smiles-drawer");
      if (cancelled || !svgRef.current) return;
      const SmilesDrawer =
        (mod as unknown as { default?: unknown }).default ?? mod;
      const lib = SmilesDrawer as {
        SvgDrawer: new (options: Record<string, unknown>) => {
          draw: (
            tree: unknown,
            target: SVGSVGElement | string,
            theme: string,
          ) => void;
        };
        parse: (
          smiles: string,
          success: (tree: unknown) => void,
          error: (err: unknown) => void,
        ) => void;
      };

      const drawer = new lib.SvgDrawer({
        width: 480,
        height,
        bondThickness: 1.1,
        bondLength: 18,
        atomVisualization: "default",
        fontSizeLarge: 6,
        fontSizeSmall: 4,
        padding: 18,
        themes: {
          light: {
            C: "#3a3f47",
            O: "#0f6b6b",
            N: "#1f4f7a",
            F: "#2f7a4f",
            CL: "#2f7a4f",
            BR: "#8a5a1f",
            I: "#6b3f7a",
            P: "#a06a1f",
            S: "#a08a1f",
            B: "#7a5a3f",
            SI: "#7a7a7a",
            H: "#5b6068",
            BACKGROUND: "transparent",
          },
        },
      });

      lib.parse(
        smiles.trim(),
        (tree) => {
          if (cancelled || !svgRef.current) return;
          try {
            svgRef.current.innerHTML = "";
            drawer.draw(tree, svgRef.current, "light");
            setParseError(false);
          } catch {
            setParseError(true);
          }
        },
        () => {
          if (cancelled) return;
          if (svgRef.current) svgRef.current.innerHTML = "";
          setParseError(true);
        },
      );
    })();

    return () => {
      cancelled = true;
    };
  }, [smiles, height]);

  return (
    <figure className={className}>
      <div className="relative overflow-hidden border border-border bg-card rule-grid">
        <svg
          ref={svgRef}
          role="img"
          aria-label={
            smiles.trim()
              ? `Two-dimensional structure drawn from SMILES ${smiles}`
              : "No molecular structure entered yet"
          }
          viewBox="0 0 480 260"
          className="h-auto w-full"
          style={{ maxHeight: height }}
        />
        {!smiles.trim() && (
          <figcaption className="absolute inset-0 flex items-center justify-center px-6 text-center text-sm text-muted-foreground">
            The structure appears here as you type a SMILES string.
          </figcaption>
        )}
        {parseError && smiles.trim() && (
          <figcaption className="absolute inset-0 flex items-center justify-center px-6 text-center text-sm text-destructive">
            This SMILES string could not be parsed into a structure.
          </figcaption>
        )}
      </div>
    </figure>
  );
}
