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
- compounds consumed by one member and produced by another in two lanes
  flanking the member column, level with the members that exchange them.

Every compound is collapsed to one shared node per condition, so a compound
several members draw on is drawn once.

### Ordered by connectivity

The member column is ordered by how many compounds each member exchanges:
the busiest member in the middle, the quietest at the two ends. The input and
output columns take the same shape — the compounds with the most edges in the
middle, the one-off compounds at the ends — and a compound goes in whichever
half of its column, top or bottom, its own members mostly occupy. The busiest
traffic therefore crosses the middle of the map on short edges, and the long
thin edges are left to the ends where they cross little.

Exchange-lane nodes are the exception: they still sit at the mean height of
the members they link, because reading cross-feeding level with the members
doing it is the point of the lanes.

### Curved edges

`marker -> metabolite` segments are drawn as S-curves (`MapStyle.edge_curve`,
`"s"` by default): each Bezier handle keeps its own endpoint's height, so an
edge leaves the marker and reaches the compound horizontally and climbs in
between. Edges sharing a marker or a compound bundle instead of fanning out
as straight diagonals, and arriving horizontally puts them on the node's
inner side, clear of the outward-running labels. `edge_curve="chord"` puts
both handles on the straight line between the endpoints, which is how the
edges were drawn before.

A steep edge has no room to make that turn — the exchange lanes are only
`mixed_lane_dx` across but can span the whole member column — so an edge is
blended back towards its chord as `|dx| / |dy|` falls below
`curve_steepness`; without that, several lane edges would run vertically at
the same x and hide each other.

### Scaling with the community

Two things scale with how large the community is, so a map of forty members
reads like a map of five rather than like a ribbon:

- the exchange lanes move out to `lane_dx_fraction` of the tallest member
  column (`mixed_lane_dx` is the floor), which keeps lane edges slanted
  enough to be told apart;
- `MapStyle.max_aspect` (5 by default) bounds the exported canvas at 1:5 and
  5:1. Blocks are tiled `n` across into a grid instead of stacked in one
  column — the fewest per row that fits — and a map still too tall has its
  lanes and columns pushed further out, which cannot collide with anything
  and leaves every height where it was. A map that is somehow too wide has
  its vertical spacings opened up instead. `max_aspect=None`
  (`--no-aspect-limit`) restores the single stacked column.

### What only the SVG can carry

Escher's JSON schema has no per-segment style and no box node — marker
circles and a floating reaction label are all it can say about a member — so
two things happen on the rendered SVG instead:

- `cross_feeding_segments` returns the ids of every segment touching a
  cross-fed compound, and `svg_editor.dash_segments` (or
  `EscherSVG_processing(dashedEdges=...)`) dashes them;
- `svg_editor.draw_member_boxes` (on by default in `EscherSVG_processing`)
  draws a labelled box on each reaction anchor — the point where its two
  short marker segments meet — and moves the reaction's own label inside it,
  so a member reads as a node instead of as the bare vertex where its edges
  happen to meet. This restores the notebook's ASV rectangles, which found
  the same segments but placed the box with offsets tuned to one figure's
  viewBox.

```python
dashed = cross_feeding_segments(json.load(open("map.json")))
EscherSVG_processing("map.svg", dashedEdges=dashed)
```

Run against `ASVMetaboliteInteractions.csv` it reproduces the reference map's
115 reactions with identical IDs, metabolite sets, and coefficients, on a
generated layout instead of hand-placed coordinates. Adding
`--min-flux 0.05 --skip-amino-acids` reproduces the reference figure's
compound selection.

Install with `pip install -e .` (or `uv pip install -e .`) from the repo
root. Console scripts `escher-edit-filter`, `escher-edit-clean`,
`escher-edit-map`, `escher-edit-svg` run each module's notebook-equivalent
entrypoint; `escher-edit-build` runs the map generator.
