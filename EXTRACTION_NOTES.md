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
| `src/escher_edit/render.py` | `render_map_svg`: map JSON -> finished SVG figure, no browser |
| `src/escher_edit/palette.py` | member colours, by member or by group (phylogeny) |
| `src/escher_edit/interactive.py` | the interactive HTML figure written beside each rendered SVG |

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
- compounds consumed by one member and produced by another in two lanes —
  secondary columns between the member column and the input and output
  columns — level with the members that exchange them.

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

`marker -> metabolite` segments are drawn as circular arcs (`MapStyle.
edge_curve`, `"arc"` by default): the handles are placed geometrically rather
than at a fraction of the run, so an edge turns at the same rate for its whole
length instead of running straight and then bending. The arc wanted is the
circle tangent to the horizontal at the marker that passes through the
compound — leaving level is what makes the edges sharing a marker bundle
rather than fan out as straight diagonals.

An edge that climbs much more than it runs cannot have that arc: all of its
turning has to happen in the narrow strip between the member and the lane, so
the circle bows out well past the compound's own column and comes back. So the
bow is capped. `arc_bulge` is how far an arc may leave the straight chord
between its ends, as a fraction of the horizontal run; past it the edge is
still a circular arc through the same two points, but a shallower one, which
tilts its departure off level. Steep edges straighten towards their chord
rather than kinking, and at the default cap an arc overshoots the column it is
heading for by at most a few percent of its own run. Lowering `arc_bulge`
flattens everything; 0 draws straight chords.

Keeping the second handle off the compound in the shallow case is not only
cosmetic: Escher's `displacedCoords` pulls a segment's endpoint back from the
node along the `b2 -> end` direction and rotates the arrowhead by it, so a
handle sitting on the node's own centre line gives a vertical arrow on a
level edge, and one closer than the displacement makes the curve overshoot
into a cusp.

`edge_curve="s"` is the earlier two-bend shape — level at both ends, with
`curve_steepness` blending steep edges back towards their chord so lane edges
do not collapse into coincident vertical runs. `edge_curve="chord"` puts both
handles on the straight line between the endpoints, which is how the edges
were drawn before either.

### Which face an edge uses

What a member consumes joins its node on the left half; what it excretes
leaves on the right half. For most edges the layout is what enforces it: the
compound is simply drawn on the side its direction asks for, and the plain arc
above joins it.

A cross-fed compound is one shared node in a lane, so it can only be on one
side, and the edges approaching from the other side cannot be satisfied by
placement. The rule says a compound sits right of every member that produces
it and left of every member that consumes it, so every cross-feeding edge
orders its two members along x, and a consistent set of orderings exists only
if the cross-feeding graph is acyclic — two members that feed each other are a
cycle. The lane rule puts the compound where the fewest edges are stranded,
and `_loop_handles` turns each stranded edge back through 180 degrees: out by
the required face, round, and across to the compound. Its partner edge needs
no turn, so the two together read as one S through the shared node.

The turn is grown smallest-first and accepted as soon as it leaves by the
*face* — not merely on the correct half, which an edge climbing over a corner
also satisfies while putting its arrowhead under the box — and turns nowhere
tighter than `_LOOP_MIN_RADIUS` (40 px against a 10 px stroke). The box it has
to clear comes from `MapStyle.member_node_half`, so `member_node_pad` and
`member_node_min_half` must describe whatever finally draws the member: fit a
turn to a smaller box than the one drawn and the arrowhead lands under its
corner, which is exactly how the previous attempt failed.

The alternative is drawing the compound once in each lane, which satisfies the
rule with no turns at all but puts the same compound on the map twice and
breaks the visible producer -> node -> consumer link.

### Scaling with the community

Two things scale with how large the community is, so a map of forty members
reads like a map of five rather than like a ribbon:

- the map widens with the tallest member column: the input and output
  columns are fitted `lane_dx_fraction` of its height out (`mixed_lane_dx`
  is the floor), plus the room the lane labels take, which keeps edges
  slanted enough to be told apart;
- the lanes then sit `lane_position` (half, by default) of the way out to
  the input and output columns, so the cross-fed compounds read as a column
  of their own, as in the 3H11/R12 SynCom figure. Before, they stayed a fixed
  label gap inside the outer columns, so as the map widened they drifted out
  to sit beside those columns — 324 px from them and 3,500 px from the
  members in a 104-member community. Where the columns are only just far
  enough out for the lane labels (a small community, or large labels) the
  lanes stay at the closest the labels allow, which is why the SynCom map
  itself is unchanged;
- `MapStyle.max_aspect` (5 by default) bounds the exported canvas at 1:5 and
  5:1. Blocks are tiled `n` across into a grid instead of stacked in one
  column — the fewest per row that fits — and a map still too tall has its
  columns pushed further out, and the lanes with them to stay between, which
  cannot collide with anything and leaves every height where it was. A map that is somehow too wide has
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
  draws a labelled box on each member and moves the reaction's own label
  inside it, so a member reads as a node instead of as the bare vertex where
  its edges happen to meet. This restores the notebook's ASV rectangles.
  The builder records on every member reaction the size of the box its
  edges were fitted around (`reaction["member_box"]`, centred on the
  midmarker; `member_box_rects` turns it into rectangles), and given the map
  (`json_path`) the boxes are drawn exactly there, labelled with the member
  name. Escher ignores the extra key on load but drops it when a map is
  re-saved from its editor; a map without it gets a box sized to its label
  on the point where its two short marker segments meet.

- `svg_editor.color_member_edges` (on by default: `member_colors=True`)
  gives every member its own colour — its edges, their arrowheads and its
  box, with the box label in white or ink, whichever reads on that fill —
  while metabolite nodes stay ink. This is the 3H11/R12 SynCom figure's
  scheme made the default. Colours come from `escher_edit.palette`: the
  first eight are a validated categorical palette (the SynCom's blue and
  orange are slots 1 and 2); past eight, each further member gets the
  in-band colour furthest from those already used, with a warning that the
  labels have to carry identity. Members are named from the map JSON when
  `json_path` is given, so a member keeps its colour in every condition
  block; pass a `palette.member_colors(...)` mapping instead of `True` to
  colour a series of per-condition maps alike. When `largeEdgeLabels` fades
  the map, a faded edge is tinted from its own member's colour; segments and
  reactions in `colorElements` keep their highlight colour.

- Colouring by group, e.g. by phylogeny: `member_groups` (`{member: group}`
  or a function of the member) makes every member of a group share one
  colour. Groups take the palette slots largest first (ties by name), so the
  bulk of the map gets the validated colours; members without a group are
  grey (`palette.UNGROUPED`, with a warning naming them) so a taxonomy gap
  cannot pass for a taxon. `palette.taxon_groups(taxonomy, rank)` reads the
  groups off a taxonomy — `{member: {rank: name}}` in any key case, with the
  names stripped (the ` Firmicutes` in `asvname_taxonomy.json`), or GTDB
  `d__...;p__...` strings. A legend below the map names the groups
  (`draw_member_legend`; the canvas grows to hold it): `member_legend` titles
  it when a string, drops it when `False`. `group_colors` (`{group: colour}`)
  fixes colours outright; build it once with `palette.group_palette(
  all_groups)` to colour a series of maps with different members alike — slot
  order otherwise depends on each map's own members.

```python
from escher_edit import render_map_svg, taxon_groups
phylum = taxon_groups(json.load(open("asvname_taxonomy.json")), "phylum")
render_map_svg("map.json", member_groups=phylum, member_legend="Phylum")
```

From the command line (`escher-edit-render`, and `escher-edit-build` for the
figures it renders): `--member-groups taxonomy.json --group-rank phylum`, or
`--member-groups groups.json` for a flat `{member: group}` object;
`--group-colors`, `--legend-title`, `--no-legend`.

```python
escher_map = json.load(open("map.json"))
EscherSVG_processing("map.svg", json_path="map.json",
                     dashedEdges=cross_feeding_segments(escher_map))
```

Escher itself can never draw the boxes, so the figures are generated here
instead: `render.render_map_svg("map.json")` (or `escher-edit-render
map.json`) writes the SVG Escher would export straight from the JSON and
runs it through `EscherSVG_processing` — member boxes, member colours and
dashed cross-feeding, all on by default. `escher-edit-build` renders a
figure beside every map it writes (`--no-figure` to skip), with the labels
drawn at the sizes the layout reserved.

Every rendered figure also gets an interactive HTML version beside it
(`<figure>.html`; `html=False` or `--no-html` to skip), made by
`interactive.write_interactive_html` from the finished SVG and the map:

- hovering an edge thickens it and lifts it over the rest, and enlarges the
  compound and member box at its two ends with their labels; hovering a
  compound or its label does that for all its edges and the members at their
  other ends; hovering a member box, for all its edges and compounds.
  Everything else fades to 12%. An edge lights only its own two ends (hover
  the compound to follow cross-feeding through it), and a highlight stays in
  its condition block, since nodes are per block.
- a tooltip names what is under the pointer: member, compound and flux for
  an edge; who produces and consumes a compound; how many compounds a member
  consumes and produces.
- scroll zooms, drag pans, `Fit` resets; click or tap pins a highlight, and
  the address then links to it (`figure.html#n36`).
- at the whole-map zoom, lit labels are brought up to 14 px on screen while
  few are lit (one member, or up to 12 compounds); larger sets are enlarged
  modestly instead of piling up, and the tooltip names them.

How it is put together, for anyone changing it: the SVG is inlined with no
`width`/`height`; each segment group gets `data-member`, `data-node` and
`data-flux` from the JSON's `from_node_id`/`to_node_id` and metabolite
coefficients; each drawn edge gets a transparent copy with a 12-screen-px
`vector-effect: non-scaling-stroke` hit stroke, in a layer under everything;
and the page lays a fixed transparent hover area over each compound (circle
plus label) and each member (box plus label). Only those hit shapes take the
pointer, so a label that grows when lit cannot grow out from under it and
flicker. Every change a highlight makes is undone exactly, style attributes
restored verbatim. The page needs no network: no scripts or fonts are
fetched. The only mark on the SVG itself is `data-reaction` on each member
box and its label.

Sizes, and how they hang together:

- labels: `MapStyle.label_font_px` 19.6 and `reaction_label_font_px` 29.4
  (the member labels, which also size the member boxes), with
  `EscherStyle` defaults to match for SVGs processed on their own; block
  captions `block_label_font_px` 84. The layout reserves room for these, so
  a figure that draws larger text has to lay the map out with larger values
  too;
- metabolite nodes: the compounds only consumed or only produced are drawn
  at `metabolite_radius` 78, and the cross-fed ones at `crossfed_radius` 60 —
  30% smaller. The builder marks them with Escher's own `node_is_primary`
  (True and False), so Escher's viewer draws the same distinction at its own
  radii; `render.render_map_svg` draws the style's. A column never packs its
  nodes closer than their diameter plus two label gaps;
- the gap between a node and its label, `label_pad`, is half the label size
  unless set, so it stays visible however large the labels are drawn;
- edges: 14 px (`EscherStyle.rxn_edge_px`). On a figure with member boxes,
  every member edge is cut back to where it crosses the box border
  (`svg_editor._trim_to_box`) instead of running on under the box to the
  markers at its centre, and the marker links wholly inside the box are
  hidden.

Run against `ASVMetaboliteInteractions.csv` it reproduces the reference map's
115 reactions with identical IDs, metabolite sets, and coefficients, on a
generated layout instead of hand-placed coordinates. Adding
`--min-flux 0.05 --skip-amino-acids` reproduces the reference figure's
compound selection.

Install with `pip install -e .` (or `uv pip install -e .`) from the repo
root. Console scripts `escher-edit-filter`, `escher-edit-clean`,
`escher-edit-map`, `escher-edit-svg` run each module's notebook-equivalent
entrypoint; `escher-edit-build` runs the map generator and
`escher-edit-render` renders maps to figures.
