# Vendored third-party JavaScript

Bundled inline into `bulk2spot_report.html` by `generate_report.py` (same
mechanism as Plotly's own JS bundle, read via `plotly.offline.get_plotlyjs()`)
so the report stays a single, self-contained, offline file -- these libraries
are never fetched over the network, at build time or when the report is
later opened. Used only for the per-figure "Download PDF" button (vector SVG
-> PDF, client-side, no server render).

| File | Library | Version | License | Source |
|---|---|---|---|---|
| `jspdf.umd.min.js` | [jsPDF](https://github.com/parallax/jsPDF) | 2.5.1 | MIT | https://cdn.jsdelivr.net/npm/jspdf@2.5.1/dist/jspdf.umd.min.js |
| `svg2pdf.umd.min.js` | [svg2pdf.js](https://github.com/yWorks/svg2pdf.js) | 2.2.1 | MIT | https://cdn.jsdelivr.net/npm/svg2pdf.js@2.2.1/dist/svg2pdf.umd.min.js |

Each file's own header carries its full copyright/license text. Load order
matters: `svg2pdf.umd.min.js` extends `jsPDF.API` and expects the global
`window.jspdf` (created by `jspdf.umd.min.js`) to already exist, so jsPDF
must be inlined first.

To update either pin: download the new version's UMD build from the URL
above with the new version number, replace the file, and update this table.
