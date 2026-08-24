#!/usr/bin/env python3

from pathlib import Path
import fitz  # PyMuPDF

ROOT = Path.cwd()

HEATMAP = ROOT / (
    "figures/figure6_ophiostoma_spatial_response/"
    "Figure6_Ophiostoma_spatial_response_heatmap.pdf"
)

ENRICHMENT = ROOT / (
    "figures/publication/final/"
    "Figure6b_spatial_extracellular_enrichment.pdf"
)

OUTDIR = ROOT / "figures/publication/final"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTPDF = OUTDIR / "Figure6_spatial_transcriptional_response.pdf"
OUTPNG = OUTDIR / "Figure6_spatial_transcriptional_response_preview.png"

for f in (HEATMAP, ENRICHMENT):
    if not f.exists():
        raise FileNotFoundError(f"Missing source PDF: {f}")

# ----------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------

MM_TO_PT = 72.0 / 25.4

page_width_mm = 180.0

left_margin_mm = 3.0
right_margin_mm = 3.0
top_margin_mm = 3.5
bottom_margin_mm = 3.5

content_width_mm = (
    page_width_mm
    - left_margin_mm
    - right_margin_mm
)

panel_gap_mm = 4.0

# ----------------------------------------------------------------------
# Read source PDFs
# ----------------------------------------------------------------------

heat_doc = fitz.open(HEATMAP)
enrich_doc = fitz.open(ENRICHMENT)

if len(heat_doc) != 1:
    raise RuntimeError("Expected heatmap PDF to contain exactly one page.")

if len(enrich_doc) != 1:
    raise RuntimeError("Expected enrichment PDF to contain exactly one page.")

heat_page = heat_doc[0]
enrich_page = enrich_doc[0]

heat_ratio = heat_page.rect.height / heat_page.rect.width
enrich_ratio = enrich_page.rect.height / enrich_page.rect.width

# Scale both to identical final width.
heat_height_mm = content_width_mm * heat_ratio
enrich_height_mm = content_width_mm * enrich_ratio

page_height_mm = (
    top_margin_mm
    + heat_height_mm
    + panel_gap_mm
    + enrich_height_mm
    + bottom_margin_mm
)

# ----------------------------------------------------------------------
# New figure canvas
# ----------------------------------------------------------------------

out = fitz.open()

page = out.new_page(
    width=page_width_mm * MM_TO_PT,
    height=page_height_mm * MM_TO_PT
)

x0 = left_margin_mm * MM_TO_PT
x1 = (page_width_mm - right_margin_mm) * MM_TO_PT

heat_y0 = top_margin_mm * MM_TO_PT
heat_y1 = (
    top_margin_mm + heat_height_mm
) * MM_TO_PT

enrich_y0 = (
    top_margin_mm
    + heat_height_mm
    + panel_gap_mm
) * MM_TO_PT

enrich_y1 = (
    top_margin_mm
    + heat_height_mm
    + panel_gap_mm
    + enrich_height_mm
) * MM_TO_PT

heat_rect = fitz.Rect(
    x0,
    heat_y0,
    x1,
    heat_y1
)

enrich_rect = fitz.Rect(
    x0,
    enrich_y0,
    x1,
    enrich_y1
)

# Embed source PDF pages directly.
page.show_pdf_page(
    heat_rect,
    heat_doc,
    0,
    keep_proportion=True
)

page.show_pdf_page(
    enrich_rect,
    enrich_doc,
    0,
    keep_proportion=True
)

# ----------------------------------------------------------------------
# Panel a label
#
# Panel b already contains its own lower-case b.
# Use standard Helvetica-equivalent PDF font here; visually it matches
# the Nimbus Sans house style closely at this small size.
# ----------------------------------------------------------------------

page.insert_text(
    fitz.Point(
        5.0 * MM_TO_PT,
        8.0 * MM_TO_PT
    ),
    "a",
    fontsize=10,
    fontname="helv",
    color=(0, 0, 0)
)

# ----------------------------------------------------------------------
# Save vector-preserving PDF
# ----------------------------------------------------------------------

out.save(
    OUTPDF,
    garbage=4,
    deflate=True
)

out.close()
heat_doc.close()
enrich_doc.close()

# ----------------------------------------------------------------------
# Preview
# ----------------------------------------------------------------------

check = fitz.open(OUTPDF)
preview_page = check[0]

pix = preview_page.get_pixmap(
    matrix=fitz.Matrix(2.0, 2.0),
    alpha=False
)

pix.save(OUTPNG)

final_width_mm = preview_page.rect.width / MM_TO_PT
final_height_mm = preview_page.rect.height / MM_TO_PT

check.close()

print("=" * 64)
print("FIGURE 6 MANUSCRIPT ASSEMBLY COMPLETE")
print("=" * 64)
print()
print(f"Final PDF: {OUTPDF.relative_to(ROOT)}")
print(f"Preview:   {OUTPNG.relative_to(ROOT)}")
print(f"Width:     {final_width_mm:.1f} mm")
print(f"Height:    {final_height_mm:.1f} mm")
print()
print("PASS: source PDFs embedded without rasterising the figure.")
