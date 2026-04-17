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
| `escher_filter_map.py` | "filtering the Escher map for fewer reactions based on the metabolite fluxes" |
| `escher_clean_json.py` | "editing the SVG Escher Map" > "updating the JSON file" |
| `escher_model_mapping.py` | builds `modelSVG_mapping.json` + consumption-edge list, categorizes segments/nodes |
| `escher_svg_editor.py` | `EscherSVG_processing` + helpers (the SVG post-processor) |
