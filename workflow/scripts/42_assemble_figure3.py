#!/usr/bin/env python3

from pathlib import Path
import fitz

ROOT = Path.cwd()

AB = ROOT / "figures/publication/figure_3ab_pca.pdf"
C = ROOT / "figures/publication/figure_3c_extent_direction_differential_expression.pdf"
D = ROOT / "results/publication/figure3/raincloud_candidates/figure3D_raincloud_integrated_v5.pdf"

OUTDIR = ROOT / "figures/publication/final"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTPDF = OUTDIR / "Figure3_global_transcriptional_response.pdf"
OUTPNG = OUTDIR / "Figure3_global_transcriptional_response_preview.png"

for path in (AB, C, D):
    if not path.exists():
        raise FileNotFoundError(f"Missing input: {path}")

MM_TO_PT = 72 / 25.4

page_width = 180 * MM_TO_PT
left_margin = 5 * MM_TO_PT
right_margin = 5 * MM_TO_PT
top_margin = 5 * MM_TO_PT
bottom_margin = 5 * MM_TO_PT
gap = 3 * MM_TO_PT

content_width = page_width - left_margin - right_margin


def natural_height(pdf_path, target_width):
    doc = fitz.open(pdf_path)
    rect = doc[0].rect
    doc.close()
    return target_width * rect.height / rect.width


h_ab = natural_height(AB, content_width)
h_c = natural_height(C, content_width)
h_d = natural_height(D, content_width)

page_height = (
    top_margin
    + h_ab
    + gap
    + h_c
    + gap
    + h_d
    + bottom_margin
)

out = fitz.open()
page = out.new_page(
    width=page_width,
    height=page_height
)

y = top_margin

for pdf_path, height in (
    (AB, h_ab),
    (C, h_c),
    (D, h_d),
):
    src = fitz.open(pdf_path)

    target = fitz.Rect(
        left_margin,
        y,
        left_margin + content_width,
        y + height,
    )

    page.show_pdf_page(
        target,
        src,
        0,
        keep_proportion=True,
        overlay=True,
    )

    src.close()
    y += height + gap

out.save(
    OUTPDF,
    garbage=4,
    deflate=True,
)
out.close()

preview = fitz.open(OUTPDF)
pix = preview[0].get_pixmap(
    matrix=fitz.Matrix(2.5, 2.5),
    alpha=False,
)
pix.save(OUTPNG)
preview.close()

print("=" * 64)
print("FIGURE 3 MANUSCRIPT ASSEMBLY COMPLETE")
print("=" * 64)
print()
print(f"Final PDF: {OUTPDF.relative_to(ROOT)}")
print(f"Preview:   {OUTPNG.relative_to(ROOT)}")
print(f"Width:     180 mm")
print(f"Height:    {page_height / MM_TO_PT:.1f} mm")
print()
print("PASS: source PDFs embedded without rasterisation.")
