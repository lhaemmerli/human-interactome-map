# Human Interactome Map

Interactive map of known physical protein–protein interactions in humans.
The most connected proteins sit at the centre; proteins can be filtered by
function and by the body system (tissues) where they are enriched.

**Live site:** https://lhaemmerli.github.io/human-interactome-map/

## Run locally

    python3 -m http.server 8000 -d docs
    # open http://localhost:8000

The site is the static `docs/` folder, published by GitHub Pages from the `main` branch.

## Rebuild the data

    ./pipeline/download_data.sh            # raw inputs (~150 MB) into data/raw/ (not in git)
    python3 pipeline/build_data.py         # default STRING confidence >= 700
    python3 pipeline/build_data.py 900     # stricter: fewer, more certain interactions

Raw inputs live in `data/raw/` (STRING v12 physical links + info, UniProt reviewed
human proteome, Human Protein Atlas `proteinatlas.tsv`).

## How it is built

- **Interactions**: STRING v12 *physical* subnetwork, combined score >= 700.
- **Position**: radius = sqrt(degree rank / N), so hubs are central and rings have equal area;
  angle = Louvain interaction module (sector), smoothed so partners sit near each other.
- **Function**: UniProt keywords + HPA protein classes grouped into 22 classes (`FUNCTIONS` in the pipeline).
- **Body system**: HPA RNA tissue specificity (enriched/enhanced tissues) mapped to organ systems;
  low-specificity genes are "Ubiquitous". Immune-cell-enriched genes also count as "Immune & blood".
- **Disease**: HPA "Disease involvement" terms (UniProt disease keywords + cancer-related genes)
  grouped into well-known conditions (`DISEASES` in the pipeline); the exact terms are shown per protein.
- **Search**: matches gene symbol, UniProt protein name and alternative names, and HPA gene synonyms.
- **3D structure** (detail panel, fetched live): AlphaFold DB model coloured by pLDDT confidence, or the
  best experimental PDB chain chosen by the PDBe "best structures" API and downloaded as a single chain from the
  RCSB model server. Rendered with 3Dmol.js (loaded on first use). Needs an internet connection.
- **Domains**: Pfam families per protein (UniProt cross-references + `pfam_clans.tsv.gz` for names).
  The detail panel ranks domains over-represented among a protein's partners
  (hypergeometric test against all annotated network proteins, Benjamini–Hochberg FDR < 5%).
- **Binding interfaces**: 3did (`data/raw/3did/3did_flat.gz`, Mar 2025) lists Pfam domain pairs seen in contact
  in PDB structures. Only inter-chain contacts are kept. An interaction is "explained" when a domain of each
  protein forms such a pair (~26% of interactions). This is an inferred interface, not proof for that exact pair.
- **Complex viewer** ("3D" next to each partner): intersects both proteins' PDBe best-structure lists to find a
  solved complex containing both; otherwise shows the 3did example structure for their domain pair. Chains are
  fetched individually from the RCSB model server; interface residues (heavy atoms within 4.5 Å) are drawn as sticks.
- **Guided tour**: opens automatically on the first visit (remembered in the browser's localStorage) and can be
  replayed with the **i** button. Steps are defined in the `TOUR` array in `docs/index.html`.

## Features

- **Map**: radial layout, hubs in the centre, slices = interaction modules. Dark theme by default.
- **Filters**: function, body system, disease; minimum number of partners.
- **Search**: gene symbol, protein name or synonym.
- **Protein panel**: AlphaFold / PDB 3D structure (AlphaFold model colourable by confidence or by
  AlphaMissense mutation sensitivity), binding interfaces, partners' common domains, partner list with
  "3D" complex viewer.
- **Gene list analysis**: paste genes, see their connections (vs. a degree-matched expectation), enriched
  functions / systems / diseases / domains, and candidate linker proteins.
- **Path finder**: shortest chain of interactions between two proteins (most confident path first).
- **Share & export**: the URL keeps the current view (`#p=TP53&d=Cancer`, `#genes=…`, `#path=INS,AKT1`,
  `#tour`); CSV downloads and PNG snapshot of the map.
- **Pathways**: Reactome pathway filter (search box), pathway list per protein, pathway enrichment for gene
  lists, and a path finder that can stay inside the selected pathway(s).
- **Activation / inhibition**: OmniPath directed, signed interactions (4,716 of 86,519) drawn as arrows
  (→ activation, ⊣ inhibition) on the map, as badges and filters in the partner list, with the net effect of
  paths and effect counts for gene lists.
- **Drugs**: approved drugs and clinical candidates per protein, loaded live from Open Targets.
- **Disease mutations**: UniProt curated pathogenic missense variants (`docs/variants.json`, loaded on demand)
  shown on the AlphaFold structure.
- **Methods & citation** (`docs/methods.html`): data versions, methods, statistics, limitations; `CITATION.cff`.
- **Beginner's guide** (`docs/learn.html`): interactive introduction to proteins, interactions and hubs.
- **English / French**: chosen from the browser language, switchable with the FR/EN button (remembered),
  or forced with `?lang=fr` / `?lang=en`. Protein descriptions stay in English (UniProt).
- **Easy picking**: while a filter, gene list, path or protein is highlighted, grey proteins ignore the mouse
  and the pointer snaps to the nearest coloured protein.

Contact: leo.hammerli@epfl.ch

## Data sources and licences

- [STRING](https://string-db.org) v12.0 — CC BY 4.0
- [UniProt](https://www.uniprot.org) — CC BY 4.0
- [Human Protein Atlas](https://www.proteinatlas.org) — CC BY-SA 3.0
- [Pfam / InterPro](https://www.ebi.ac.uk/interpro/) — CC0
- [3did](https://3did.irbbarcelona.org) — domain–domain interactions from PDB structures
- [Reactome](https://reactome.org) v98 — CC BY 4.0
- [OmniPath](https://omnipathdb.org) (Türei et al., Mol Syst Biol 2021) — academic-licence set, mixed licences
- [Open Targets Platform](https://platform.opentargets.org) 26.09 (live) — CC0
- [AlphaMissense](https://alphamissense.hegelab.org) (Cheng et al., Science 2023), via AlphaFold DB — CC BY 4.0
- Structures are loaded live from [AlphaFold DB](https://alphafold.ebi.ac.uk) (CC BY 4.0),
  [PDBe](https://www.ebi.ac.uk/pdbe/) and [RCSB PDB](https://www.rcsb.org)

`docs/data.json` is derived from these sources; please cite them if you reuse it.
