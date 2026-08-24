#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"

AB="$ROOT/figures/publication/figure_3ab_pca.pdf"
C="$ROOT/figures/publication/figure_3c_extent_direction_differential_expression.pdf"
D="$ROOT/results/publication/figure3/raincloud_candidates/figure3D_raincloud_integrated_v5.pdf"

OUTDIR="$ROOT/figures/publication/final"
mkdir -p "$OUTDIR"

for f in "$AB" "$C" "$D"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing input: $f" >&2
        exit 1
    fi
done

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

cat > "$TMPDIR/figure3.tex" <<EOF
\documentclass{article}

\usepackage[paperwidth=180mm,paperheight=220mm,
            margin=0mm]{geometry}
\usepackage{graphicx}
\usepackage{grffile}
\pagestyle{empty}

\begin{document}
\noindent
\begin{minipage}[t]{180mm}

\includegraphics[width=\linewidth,height=78mm,keepaspectratio]{$AB}

\vspace{2mm}

\includegraphics[width=\linewidth,height=60mm,keepaspectratio]{$C}

\vspace{2mm}

\includegraphics[width=\linewidth,height=72mm,keepaspectratio]{$D}

\end{minipage}
\end{document}
EOF

cd "$TMPDIR"

pdflatex \
    -interaction=nonstopmode \
    -halt-on-error \
    figure3.tex >/dev/null

cp figure3.pdf \
   "$OUTDIR/Figure3_global_transcriptional_response.pdf"

echo
echo "============================================================"
echo "FIGURE 3 ASSEMBLY COMPLETE"
echo "============================================================"
echo
echo "Input AB:"
echo "  $AB"
echo
echo "Input C:"
echo "  $C"
echo
echo "Input D:"
echo "  $D"
echo
echo "Final figure:"
echo "  $OUTDIR/Figure3_global_transcriptional_response.pdf"
