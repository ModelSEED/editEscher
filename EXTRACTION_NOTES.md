# Extraction notes

This repo holds the Escher-editing code extracted from
`ABX_mouse_gut/Notebooks/escher_API_mapping.ipynb`. The extraction was not a
verbatim copy — the following changes were made while moving the code out of
the notebook.

## Structural

- Wrapped top-level script blocks into functions:
  `filter_escher_map`, `build_node_lookups`, `build_direction_tracking`,
  `categorize_segments_and_nodes`, `build_name_abbrev_table`.
- Added `if __name__ == "__main__":` blocks that replicate the notebook's
  concrete call.
- Split the code across four files with cross-module imports
  (e.g. `from escher_clean_json import CPD_ID_ABBRV`).

## Cosmetic

- Renamed `cpdID_abbrv` → `CPD_ID_ABBRV`, `aa_names` → `AA_NAMES` to match
  Python constant conventions.
- Added module- and function-level docstrings that did not exist in the
  notebook.
- Reformatted long lines, normalized to double quotes, cleaned indentation.
- Dropped most commented-out debug code (the long commented block inside the
  filter loop, `# display(...)` / `# print(...)` lines).
- Removed a few stray variables that weren't used downstream
  (`og_ele`, `nodes_to_remove` inside the filter loop).
- Moved the `shapely` import into the helper functions that use it and left
  it commented at module top.

## Small logic touch-ups

- In `filter_escher_map`: the original did `idName[content["bigg_id"]]`,
  which would raise `KeyError` on nodes that don't carry a `bigg_id`.
  Changed to `idName.get(content.get("bigg_id", ""))`.
- `nodeNumCoef.update({...})` is now guarded with `if bid in metCoefs`,
  since not every node referenced by a reaction's segments is in that
  reaction's metabolite list.

If a verbatim extraction is preferred — same variable names, same
commented-out code, same top-level script style — these files can be
regenerated from the notebook directly.

## File map

| File | Purpose (notebook section) |
| --- | --- |
| `src/escher_edit/build_map.py` | *(new, not from the notebook)* builds an Escher map from community exchange fluxes — see below |
| `src/escher_edit/filter_map.py` | "filtering the Escher map for fewer reactions based on the metabolite fluxes" |
| `src/escher_edit/clean_json.py` | "editing the SVG Escher Map" > "updating the JSON file" |
| `src/escher_edit/model_mapping.py` | builds `modelSVG_mapping.json` + consumption-edge list, categorizes segments/nodes |
| `src/escher_edit/svg_editor.py` | `EscherSVG_processing` + helpers (the SVG post-processor) |

## Map generation

`build_map.py` is not an extraction — it closes the gap the notebook left
open. Everything else in this package edits a map that was drawn by hand in
the Escher web editor; `build_map.py` generates one from the per-member
exchange fluxes computed in
`MicrobiomeNotebooks/NewWesternDiet/ASVCommunityModeling.ipynb` (cell 49) and
sliced by `processing_CSVs.ipynb` (cells 51-52).

Each community member becomes one net organism reaction — consumed compounds
(negative flux) as reactants, excreted compounds (positive flux) as products.
Each reaction keeps the marker topology of `metabolite_focused_map_IDs.json`
(a midmarker, two multimarkers 20 px away, `multimarker -> metabolite`
segments), and the reactions are arranged as three columns:

- member reactions stacked in one column with their midmarkers on a single
  vertical axis;
- compounds only ever consumed in a left input column;
- compounds only ever produced in a right output column;
- compounds consumed by one member and produced by another in two narrow
  lanes flanking the member column, level with the members that exchange
  them.

Every compound is collapsed to one shared node per condition, so a compound
several members draw on is drawn once. Nodes are placed at the mean height of
the members they connect to; the outer columns are then spread over the full
height of the member column, while the exchange lanes stay level with their
members.

Run against `ASVMetaboliteInteractions.csv` it reproduces the reference map's
115 reactions with identical IDs, metabolite sets, and coefficients, on a
generated layout instead of hand-placed coordinates. Adding
`--min-flux 0.05 --skip-amino-acids` reproduces the reference figure's
compound selection.

Install with `pip install -e .` (or `uv pip install -e .`) from the repo
root. Console scripts `escher-edit-filter`, `escher-edit-clean`,
`escher-edit-map`, `escher-edit-svg` run each module's notebook-equivalent
entrypoint; `escher-edit-build` runs the map generator.
