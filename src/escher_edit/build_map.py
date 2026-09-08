"""Build an Escher map directly from community exchange-flux data.

Every community member becomes one Escher *reaction* — its **net organism
reaction** — where the compounds the member consumes are the reactants
(negative flux) and the compounds it excretes are the products (positive
flux). Each reaction keeps the marker topology of the clusters in
``metabolite_focused_map_IDs.json``:

* one ``midmarker`` node at the reaction centre;
* two ``multimarker`` nodes at ``±marker_offset`` px — the left one carries
  every reactant, the right one every product;
* segments ``multimarker -> midmarker`` (``b1``/``b2`` = ``None``) plus one
  ``multimarker -> metabolite`` segment per compound, with bezier handles.

The members are stacked in a single column with their midmarkers aligned on
one vertical axis, and every compound is collapsed to **one shared node per
condition**:

* compounds only ever consumed sit in a left input column;
* compounds only ever produced sit in a right output column;
* compounds consumed by one member and produced by another — the cross-fed
  ones — sit in two lanes flanking the member column, level with the members
  that exchange them. The lanes stand off as far as the member column is tall
  (``lane_dx_fraction`` of it), so their edges keep enough slant to be told
  apart however many members feed into them.

Everything is ordered by connectivity. The members that exchange the most
compounds sit in the middle of the member column and the quietest at its two
ends, and the input and output columns take the same shape — the compounds
with the most edges in the middle, the one-off compounds at the ends — with
each compound in the half of its column, top or bottom, that its own members
mostly occupy. So the busiest traffic crosses the middle of the map on short
edges, and the long thin edges are left to the ends where they cross little.
An exchange-lane node still sits at the mean height of the members it links,
because reading cross-feeding level with the members doing it is the whole
point of the lanes.

Edges are drawn as single arcs: one bend, no flat run at either end. The
first Bezier handle stays level with the marker, so an edge leaves the member
horizontally and the edges sharing a marker bundle instead of fanning out as
straight diagonals; the second handle sits short of the compound's own
height, which bends the arc once on its way there and leaves it running into
the compound at an angle. ``MapStyle.arc_arrival`` sets how far short, and so
how much the arc bows.

Which face of a member an edge uses is fixed: what the member consumes
arrives on the left of its node, what it excretes leaves on the right. A
cross-fed compound is one shared node, though, so the members on the far
side of it are reaching across, and their edges have to leave by their own
face and come back. Those — and only those — are drawn with a longer first
handle, swinging out of the face and round the node in one sweep. Every
other edge on the map is the plain curve above.

Blocks are tiled into a grid rather than stacked in one column, and a block
that is still too tall gets its input and output columns pushed outwards, so
the exported canvas stays inside ``MapStyle.max_aspect`` (5:1 by default) no
matter how many members or conditions the community has.

Reaction ``bigg_id`` is ``f"{member}{model_id}"`` (e.g.
``Acetivibrio.1RC-ABX_12.5``) and ``name`` is the bare member name, matching
the reference map so that :mod:`escher_edit.model_mapping`,
:mod:`escher_edit.clean_json`, and :mod:`escher_edit.filter_map` all keep
working on the output.

Input
-----
The per-member exchange fluxes produced by
``MicrobiomeNotebooks/NewWesternDiet/ASVCommunityModeling.ipynb`` (cell 49,
"Printing ASV results for analysis"), which sums each member's ``_e0``
transport reactions with consumption negative and excretion positive:

* ``ASVMetaboliteInteractions.csv`` — wide, one row per member, columns
  ``<cpd>_<model_id>`` (plus a ``<cpd>_ave`` block). This is the preferred
  input because its columns already carry ModelSEED compound IDs.
* ``<diet>_<day>_communityFluxes.csv`` — the per-condition slices written by
  ``processing_CSVs.ipynb`` (cells 51-52). Columns are compound *names*, so a
  name->ID mapping is needed to emit ``cpd`` bigg_ids.

An empty cell means the member is absent from that condition; an explicit
``0`` means it is present but does not exchange that compound. Both are left
out of the map, and the two are counted separately in the log.
"""
import argparse
import copy
import csv
import logging
import math
import uuid
from json import dump
from pathlib import Path

log = logging.getLogger(__name__)

ESCHER_SCHEMA = "https://escher.github.io/escher/jsonschema/1-0-0#"

#: Sentinel condition whose reaction ``bigg_id``s do not carry the
#: ``<diet>-ABX_<day>`` suffix, so ``model_mapping.build_direction_tracking``
#: cannot parse a map built from it.
AVERAGE_CONDITION = "ave"

#: Radii Escher renders nodes at, used to keep a node's own circle inside the
#: canvas when the block extents are measured.
NODE_RADIUS = 20.0
MARKER_RADIUS = 5.0

#: How many times :func:`build_escher_map` may re-lay a map out while fitting
#: it to ``MapStyle.max_aspect``. Widening a block hits the target exactly, so
#: the extra passes only cover the grid changing shape underneath it.
ASPECT_PASSES = 5

#: Slack on the aspect test that ends the fitting loop, so a canvas landing a
#: rounding error short of the limit is not sent round again. What is left is
#: made up with canvas padding, which is exact.
ASPECT_TOLERANCE = 1e-6

#: How much of a dimension that padding may quietly add. More than this and
#: the layout genuinely would not stretch, which is worth a warning.
ASPECT_PAD_WARNING = 1e-3


# --------------------------------------------------------------- parsing


def _to_float(text):
    """Return ``float(text)`` or None for blank/``nan``/non-numeric cells."""
    text = (text or "").strip()
    if not text or text.lower() in ("nan", "na", "none"):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_interaction_matrix(csv_path):
    """Parse a wide ``ASVMetaboliteInteractions.csv``.

    Columns after the first are named ``<compound_id>_<model_id>``; the
    model id itself contains underscores (``RC-ABX_-1.5``), so only the
    first underscore separates the two fields.

    Returns
    -------
    dict[str, dict[str, dict[str, float]]]
        ``{model_id: {member: {compound_id: flux}}}``. Members that are
        absent from a condition (all-blank row) are omitted from it, and
        blank / zero cells are omitted from a member's flux dict.
    """
    csv_path = Path(csv_path)
    with csv_path.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        fields = [(i, *col.split("_", 1))
                  for i, col in enumerate(header[1:], start=1)
                  if "_" in col]

        data = {}
        absent = present_zero = kept = 0
        for row in reader:
            if not row or not row[0].strip():
                continue
            member = row[0].strip()
            seen_any = {}
            for i, cpd, model_id in fields:
                value = _to_float(row[i]) if i < len(row) else None
                bucket = data.setdefault(model_id, {})
                if value is None:
                    absent += 1
                    continue
                seen_any[model_id] = True
                if value == 0:
                    present_zero += 1
                    continue
                bucket.setdefault(member, {})[cpd] = value
                kept += 1
            for model_id in seen_any:
                data.setdefault(model_id, {}).setdefault(member, {})

    log.info("%s: %d conditions; %d exchange fluxes kept, "
             "%d present-but-zero, %d absent (blank)",
             csv_path.name, len(data), kept, present_zero, absent)
    return data


def load_compound_names(id_csv, names_csv):
    """Map compound ID -> display name by zipping two header rows.

    ``ASVMetaboliteInteractions.csv`` and
    ``ASVMetaboliteInteractions_names.csv`` are the same matrix with the
    headers written as IDs and as names respectively, so the mapping is
    positional. Raises ``ValueError`` if the headers disagree in length or
    if a compound maps to two different names.
    """
    def _header(path):
        with Path(path).open(newline="") as fh:
            return next(csv.reader(fh))[1:]

    ids, names = _header(id_csv), _header(names_csv)
    if len(ids) != len(names):
        raise ValueError(
            f"header length mismatch: {len(ids)} ids vs {len(names)} names")

    mapping = {}
    for id_col, name_col in zip(ids, names):
        cpd = id_col.split("_", 1)[0]
        # pandas de-duplicates repeated headers as "Name.1", "Name.2", ...
        name = name_col.rsplit(".", 1)[0] if name_col.rsplit(".", 1)[-1].isdigit() else name_col
        if mapping.get(cpd, name) != name:
            raise ValueError(
                f"{cpd} maps to both {mapping[cpd]!r} and {name!r}")
        mapping[cpd] = name
    log.info("resolved %d compound names", len(mapping))
    return mapping


def parse_condition_matrix(csv_path, name_to_id=None):
    """Parse one ``<diet>_<day>_communityFluxes.csv`` slice.

    Columns are compound *names*. When ``name_to_id`` is given, names are
    translated to compound IDs (unmapped names are kept verbatim so nothing
    is silently dropped).

    Returns
    -------
    dict[str, dict[str, float]]
        ``{member: {compound: flux}}``
    """
    name_to_id = name_to_id or {}
    csv_path = Path(csv_path)
    with csv_path.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        cols = [(i, name_to_id.get(col, col))
                for i, col in enumerate(header[1:], start=1)]

        members = {}
        for row in reader:
            if not row or not row[0].strip():
                continue
            fluxes = {}
            for i, cpd in cols:
                value = _to_float(row[i]) if i < len(row) else None
                if value:
                    fluxes[cpd] = value
            members[row[0].strip()] = fluxes
    log.info("%s: %d members", csv_path.name, len(members))
    return members


# ------------------------------------------------------- net reactions


def build_member_reactions(fluxes_by_member, model_id="",
                           compound_names=None,
                           min_abs_flux=0.0, skip_names=None,
                           skip_ids=None, drop_empty=True):
    """Turn ``{member: {compound: flux}}`` into net organism reactions.

    Parameters
    ----------
    fluxes_by_member : dict[str, dict[str, float]]
        Net exchange flux per member per compound; negative = consumed,
        positive = excreted.
    model_id : str
        Condition suffix appended to each member name to form the reaction
        ``bigg_id`` (e.g. ``"RC-ABX_12.5"`` -> ``Acetivibrio.1RC-ABX_12.5``).
        Pass ``""`` to use the bare member name.
    compound_names : dict[str, str] or None
        Compound ID -> display name, used for node ``name`` fields and for
        matching ``skip_names``.
    min_abs_flux : float
        Drop compounds whose ``abs(flux)`` is below this. The reference
        figure used 0.05.
    skip_names : iterable[str] or None
        Display names to exclude, e.g.
        ``escher_edit.filter_map.DEFAULT_SKIP_NAMES``.
    skip_ids : iterable[str] or None
        Compound IDs to exclude.
    drop_empty : bool
        Omit members left with no compounds after filtering.

    Returns
    -------
    list[dict]
        ``[{"name", "bigg_id", "fluxes": {compound: flux}}, ...]`` sorted by
        member name.
    """
    compound_names = compound_names or {}
    skip_names = set(skip_names or ())
    skip_ids = set(skip_ids or ())

    members = []
    for member in sorted(fluxes_by_member):
        fluxes = {}
        for cpd, flux in fluxes_by_member[member].items():
            if not flux or abs(flux) < min_abs_flux:
                continue
            if cpd in skip_ids or compound_names.get(cpd, cpd) in skip_names:
                continue
            fluxes[cpd] = flux
        if not fluxes and drop_empty:
            continue
        members.append({
            "name": member,
            "bigg_id": f"{member}{model_id}",
            "fluxes": fluxes,
        })
    log.info("condition %r: %d member reactions "
             "(%d reactants, %d products total)",
             model_id or "<none>", len(members),
             sum(1 for m in members for f in m["fluxes"].values() if f < 0),
             sum(1 for m in members for f in m["fluxes"].values() if f > 0))
    return members


# ------------------------------------------------------------- geometry


def _single_bend_limit(bezier_fracs):
    """The largest ``arc_arrival`` that still bends the edge only once.

    An arc is one bend for as long as its control polygon keeps turning the
    same way. With the handles at ``bezier_fracs`` along the run and the
    first held level with the marker, that holds while the second has risen
    less than this much of the way to the compound; past it the edge turns
    back on itself near the compound and the arc becomes a shallow S.
    """
    f1, f2 = bezier_fracs
    return (f2 - f1) / (1.0 - f1)


def _label_width(text, font_px):
    """Rough on-map width of a label, in px.

    Escher anchors label text at ``label_x`` and runs it to the right, so this
    is how far past that anchor a label reaches. 0.6 em per character is the
    usual approximation for the sans-serif face Escher renders with.
    """
    return 0.6 * font_px * len(str(text))


class MapStyle:
    """Column layout for a community exchange map.

    Member reactions are stacked in one vertical column with their midmarkers
    aligned on a single axis. Every compound is collapsed to one node per
    condition block, placed in one of three regions:

    * **inputs** — consumed by at least one member and produced by none — in
      a left column at ``input_column_dx``;
    * **outputs** — produced but never consumed — in a right column at
      ``output_column_dx``;
    * **exchanged** — consumed by one member and produced by another — in two
      lanes flanking the member column, level with the members that consume
      and produce them.

    ``mixed_lane_dx`` is a floor, not the offset: a lane is only as slanted
    as its width lets it be, so a lane that keeps 240 px while the member
    column grows to twenty thousand ends up with every one of its edges
    running vertically on top of the others. The offset therefore scales with
    the tallest member column in the map — ``lane_dx_fraction`` of it — which
    keeps the lane edges as separable in a forty-member community as in a
    five-member one (see :meth:`fitted_lane_dx`). Set ``lane_dx_fraction`` to
    0 to stop the lanes growing with the column; the aspect fit below may
    still move them.

    Members are ordered by connectivity — the ones exchanging the most
    compounds in the middle of the column, the quietest at its two ends — and
    the input and output columns take the same shape, ``min_node_spacing``
    apart at the closest (see :func:`_connectivity_column`). Which half of
    its column a compound goes in, top or bottom, is decided by where its own
    members lie, so it still sits on the side its edges come from.

    An exchange-lane node instead sits at the mean height of the members it
    links to, pushed apart as far as ``lane_node_spacing`` requires: being
    level with those members is the point of drawing it there. The lanes get
    the wider spacing default because they hold the busiest compounds, whose
    barycentres all fall near the middle of the member column.
    ``lane_span_fraction`` keeps a lane from growing into a fourth
    full-height column: when the packed lane would exceed that fraction of
    the member column's height, its spacing is compressed to fit, still
    centred on the lane's barycentre.

    ``span_columns`` makes the input and output columns run the full height of
    the member column — they read as columns rather than as a knot in the
    middle. Without it each half packs tight around the middle instead.

    ``input_column_dx`` and ``output_column_dx`` default to ``None``, which
    fits the columns as close to the member column as the exchange-lane
    labels allow (see :meth:`fitted_column_dx`). Pass numbers to place them
    by hand.

    Edges
    -----
    ``edge_curve`` picks the shape of the ``marker -> metabolite`` segments.
    Both handles sit at ``bezier_fracs`` along the run; what differs is how
    far each has risen towards the compound:

    * ``"arc"`` (the default) holds the first handle level with the marker
      and stops the second ``arc_arrival`` of the way up, which bends the
      edge exactly once — a single bow with no flat run at either end.
      Leaving level means the edges sharing a marker still bundle rather than
      fanning out as straight diagonals; stopping the second handle short is
      what keeps the bend single, and lowering ``arc_arrival`` deepens the
      bow. Above :func:`_single_bend_limit` the arc picks up a second bend
      and stops being an arc, so that is the ceiling.
    * ``"s"`` holds each handle at its own endpoint's height: the edge leaves
      the marker horizontally, climbs in the middle, and arrives at the
      compound horizontally — two bends. A steep edge has no room for that
      turn, so ``curve_steepness`` is the ``|dx|/|dy|`` at which one still
      gets the full curve; below it the handles blend back towards the chord,
      which leaves the steep exchange-lane edges nearly straight.
    * ``"chord"`` puts both handles on the straight line between the
      endpoints, reproducing the flat edges of earlier maps.

    Which face of the member node an edge uses is not a matter of style: a
    consumed compound joins on the left, an excreted one leaves on the right,
    whatever ``edge_curve`` says. Only an edge whose compound lies the other
    way is drawn differently, and only in its first handle, which is swung
    out of the required face far enough to carry the curve clear of the node
    (see :func:`_hook`). ``member_node_pad`` is how much room the member's
    label is given inside its box, and should match the ``member_box_pad``
    given to ``svg_editor.draw_member_boxes``, because that box is the node
    those edges have to clear.

    Aspect ratio
    ------------
    ``max_aspect`` bounds the shape of the exported canvas: the figure is
    never narrower than ``1:max_aspect`` nor wider than ``max_aspect:1``,
    however many members or conditions it holds. The bound is symmetric, so
    it can be given either way up — 5 and 0.2 ask for the same thing. Two
    mechanisms get it there, in order:

    * blocks are tiled ``n`` across into a grid rather than stacked in one
      column (see :func:`_choose_n_cols`), with ``block_column_gap`` between
      grid columns and ``block_gap`` between rows;
    * a map that is still too tall — one tall condition on its own, say — has
      its lanes and its input and output columns pushed outwards together,
      which keeps the room the labels between them need, can never collide
      with anything, and leaves every height untouched. A block whose
      compounds are all cross-fed has no column nodes to move, which is why
      the lanes move too. Positions given by hand through
      ``input_column_dx``, ``output_column_dx`` or ``mixed_lane_dx`` are
      starting points for this, not exemptions from it.

    A map that is somehow too *wide* has its vertical spacings opened up
    instead, so the height it gains carries nodes rather than blank canvas.
    Set ``max_aspect`` to ``None`` to place blocks in a single column and
    leave the proportions alone.
    """

    def __init__(self, member_pitch=420.0, input_column_dx=None,
                 output_column_dx=None, mixed_lane_dx=240.0,
                 lane_dx_fraction=0.08, column_clearance=160.0,
                 min_node_spacing=90.0, lane_node_spacing=200.0,
                 lane_span_fraction=0.6,
                 span_columns=True, marker_offset=20.0,
                 member_label_gap=60.0, label_pad=22.0, label_font_px=20.0,
                 reaction_label_font_px=30.0, bezier_fracs=(0.25, 0.75),
                 edge_curve="arc", arc_arrival=0.6, curve_steepness=0.5,
                 member_node_pad=12.0,
                 max_aspect=5.0, block_gap=900.0, block_column_gap=700.0,
                 block_label_gap=280.0, block_label_font_px=60.0,
                 canvas_margin=400.0):
        if edge_curve not in ("arc", "s", "chord"):
            raise ValueError(
                f"edge_curve must be 'arc', 's' or 'chord', not "
                f"{edge_curve!r}")
        limit = _single_bend_limit(bezier_fracs)
        if not 0 < arc_arrival < limit:
            raise ValueError(
                f"arc_arrival must be between 0 and {limit:.3g} for "
                f"bezier_fracs {bezier_fracs}: at or above that the arc picks "
                f"up a second bend, which is the shape it exists to avoid")
        if max_aspect is not None and max_aspect <= 0:
            raise ValueError(f"max_aspect must be positive: {max_aspect!r}")
        if max_aspect is not None and max_aspect < 1:
            # the bound is symmetric, so "1:5" written as 0.2 asks for the
            # same band as 5 and is taken to mean it
            max_aspect = 1.0 / max_aspect
        self.member_pitch = member_pitch
        self.span_columns = span_columns
        self.input_column_dx = input_column_dx
        self.output_column_dx = output_column_dx
        self.mixed_lane_dx = mixed_lane_dx
        self.lane_dx_fraction = lane_dx_fraction
        self.column_clearance = column_clearance
        self.min_node_spacing = min_node_spacing
        self.lane_node_spacing = lane_node_spacing
        self.lane_span_fraction = lane_span_fraction
        self.marker_offset = marker_offset
        self.member_label_gap = member_label_gap
        self.label_pad = label_pad
        self.label_font_px = label_font_px
        self.reaction_label_font_px = reaction_label_font_px
        self.bezier_fracs = bezier_fracs
        self.edge_curve = edge_curve
        self.arc_arrival = arc_arrival
        self.curve_steepness = curve_steepness
        self.member_node_pad = member_node_pad
        self.max_aspect = max_aspect
        self.block_gap = block_gap
        self.block_column_gap = block_column_gap
        self.block_label_gap = block_label_gap
        self.block_label_font_px = block_label_font_px
        self.canvas_margin = canvas_margin

    def member_node_half(self, name):
        """Half the width and height a member is drawn at, in px.

        A member reaction is rendered as a labelled box by
        ``svg_editor.draw_member_boxes``, so the layout has to know how much
        room that box takes: its own label plus ``member_node_pad``, which
        should match the ``member_box_pad`` the SVG step is given. Edges use
        it to leave by the correct face of the box rather than crossing it.
        """
        return (_label_width(name, self.reaction_label_font_px) / 2
                + self.member_node_pad,
                self.reaction_label_font_px / 2 + self.member_node_pad)

    def fitted_lane_dx(self, column_height):
        """Offset for the exchange lanes given the tallest member column.

        A lane edge crosses ``lane_dx`` horizontally however far it climbs, so
        the taller the member column the steeper — and the more overplotted —
        every edge in the lane becomes. Scaling the offset with the column
        keeps their slant, and so their separation, roughly constant as the
        community grows; ``mixed_lane_dx`` is the floor for the small
        communities that do not need the room.
        """
        if not self.lane_dx_fraction:
            return self.mixed_lane_dx
        return max(self.mixed_lane_dx, self.lane_dx_fraction * column_height)

    def fitted_column_dx(self, compounds, lane_dx=None):
        """Tightest offset for the input/output columns, in px.

        An exchange-lane node's label runs outwards from the lane, so the
        outer columns can come in no closer than the far edge of those labels
        plus ``column_clearance``. Fitting to that keeps the figure compact
        without the lane labels colliding with the column nodes.
        """
        widest = max((_label_width(compound, self.label_font_px)
                      for compound in compounds),
                     default=_label_width("x" * 8, self.label_font_px))
        lane_label_edge = ((self.mixed_lane_dx if lane_dx is None else lane_dx)
                           + self.label_pad + widest)
        return lane_label_edge + self.column_clearance

    def column_positions(self, compounds, lane_dx=None):
        """``(input_x, output_x)`` — the fitted offsets, or the overrides."""
        fitted = self.fitted_column_dx(compounds, lane_dx)
        return (self.input_column_dx if self.input_column_dx is not None else -fitted,
                self.output_column_dx if self.output_column_dx is not None else fitted)

    def vertically_scaled(self, factor):
        """A copy with every height-setting spacing multiplied by ``factor``.

        Only the spacings that decide how tall a block comes out are touched —
        the member pitch, the two node spacings, and the gap between rows of
        blocks. The label gaps are left alone because they are set by the size
        of the text, not by how much room the figure has.
        """
        scaled = copy.copy(self)
        scaled.member_pitch *= factor
        scaled.min_node_spacing *= factor
        scaled.lane_node_spacing *= factor
        scaled.block_gap *= factor
        return scaled


def classify_compounds(members):
    """Split a block's compounds by how the community uses them.

    Returns
    -------
    tuple[list[str], list[str], list[str]]
        ``(inputs, outputs, exchanged)`` — consumed-only, produced-only, and
        both-consumed-and-produced compound IDs. ``exchanged`` is the set
        that is cross-fed between members, so those nodes are drawn beside
        the member column rather than out at the edges.
    """
    consumed, produced = set(), set()
    for member in members:
        for compound, flux in member["fluxes"].items():
            if flux < 0:
                consumed.add(compound)
            elif flux > 0:
                produced.add(compound)
    return (sorted(consumed - produced), sorted(produced - consumed),
            sorted(consumed & produced))


def cross_feeding_segments(escher_map, prefix="s"):
    """Segment ids for the edges that carry cross-feeding.

    A metabolite node is cross-fed when the reactions touching it disagree on
    sign — at least one consumes it and at least one produces it. Those are
    exactly the nodes drawn in the exchange lanes, so every segment ending on
    one is an edge along which one member feeds another.

    Ids come back ``s``-prefixed to match what Escher writes into the rendered
    SVG, ready to hand to ``svg_editor.EscherSVG_processing(dashedEdges=...)``
    or ``svg_editor.dash_segments``. Escher's JSON schema has no per-segment
    style, so dashing can only happen on the rendered SVG.

    Works on any Escher map whose nodes are shared between reactions,
    including hand-drawn ones. Direction comes from the coefficient signs —
    Escher's segments run marker -> metabolite for reactants and products
    alike, so the topology alone does not say which is which — and that means
    a node's ``bigg_id`` has to match the ids in its reactions'
    ``metabolites`` lists. Maps cleaned by the notebook's original
    ``cleanEscherJSON`` renamed the nodes but not the metabolite lists; such a
    map is logged as a warning rather than quietly returning nothing.
    """
    nodes = escher_map[1]["nodes"]
    reactions = escher_map[1]["reactions"]

    node_ids = {node.get("bigg_id") for node in nodes.values()
                if node.get("node_type") == "metabolite"}
    metabolite_ids = {met["bigg_id"] for reaction in reactions.values()
                      for met in reaction.get("metabolites", [])}
    if node_ids and not (node_ids & metabolite_ids):
        log.warning(
            "no metabolite node bigg_id matches any reaction metabolite id "
            "(e.g. node %r vs reaction metabolite %r) — this map's node names "
            "were rewritten without updating its metabolites lists, so "
            "cross-feeding cannot be read from it; use the map from before "
            "that rewrite",
            sorted(node_ids)[0], sorted(metabolite_ids)[0] if metabolite_ids else None)

    signs = {}
    for reaction in reactions.values():
        coefficients = {met["bigg_id"]: met["coefficient"]
                        for met in reaction.get("metabolites", [])}
        for segment in reaction.get("segments", {}).values():
            for node_id in (segment["from_node_id"], segment["to_node_id"]):
                node = nodes.get(node_id) or {}
                if node.get("node_type") != "metabolite":
                    continue
                if node.get("bigg_id") in coefficients:
                    signs.setdefault(node_id, set()).add(
                        coefficients[node["bigg_id"]] > 0)

    cross_fed = {node_id for node_id, seen in signs.items() if len(seen) > 1}
    found = set()
    for reaction in reactions.values():
        for segment_id, segment in reaction.get("segments", {}).items():
            if (segment["from_node_id"] in cross_fed
                    or segment["to_node_id"] in cross_fed):
                found.add(segment_id)

    log.info("%d cross-fed compound node(s) on %d segment(s)",
             len(cross_fed), len(found))
    return [f"{prefix}{seg}" for seg in
            sorted(found, key=lambda s: (len(s), s))]


def _chord_blend(dx, dy, steepness):
    """How much of the straight chord a segment keeps, in ``[0, 1]``.

    0 is a full S-curve — the handles stay level with their own endpoints —
    and 1 is the flat chord. An edge with room to turn, ``|dx| >= steepness *
    |dy|``, gets the full curve; from there the blend climbs to 1 as the run
    goes vertical, because edges that are nearly vertical would otherwise all
    curve onto the same x and hide each other.
    """
    reach = steepness * abs(dy)
    if reach <= 0:
        return 0.0
    return max(0.0, 1.0 - abs(dx) / reach)


def _sign(value):
    """-1, 0 or 1."""
    return (value > 0) - (value < 0)


def _cubic(p0, p1, p2, p3, t):
    """The point at ``t`` on the cubic Bezier through those four."""
    u = 1.0 - t
    return (u * u * u * p0[0] + 3 * u * u * t * p1[0]
            + 3 * u * t * t * p2[0] + t * t * t * p3[0],
            u * u * u * p0[1] + 3 * u * u * t * p1[1]
            + 3 * u * t * t * p2[1] + t * t * t * p3[1])


def _leaves_by(p0, p1, p2, p3, node, side, steps=96):
    """Does this curve cross the member node's edge on ``side`` and stay out?

    The member is drawn as a box, so an edge obeys the in-on-the-left,
    out-on-the-right rule by where it crosses that box's outline — not by
    where its marker sits, which is buried inside the box either way.
    """
    cx, cy, half_w, half_h = node

    def inside(point):
        return abs(point[0] - cx) <= half_w and abs(point[1] - cy) <= half_h

    out = False
    for step in range(1, steps + 1):
        point = _cubic(p0, p1, p2, p3, step / steps)
        if not out:
            if not inside(point):
                if (point[0] - cx) * side < 0:
                    return False        # left by the wrong face
                out = True
        elif inside(point):
            return False                # came back in through the node
    return out


def _hook(p0, p3, handle2, side, node, f1, growth=(1, 2, 4, 8, 16, 32)):
    """First handle for an edge whose compound is on the far side of its member.

    A cross-fed compound is one shared node, so whichever lane it sits in,
    the members on the other side of it are reaching across — and they still
    have to leave by their own face: consumed on the left, excreted on the
    right. Pointing the handle out of that face is not enough on its own,
    because a cubic only travels about a third of the way towards its handle:
    a compound level with the member would see the edge turn back through the
    box it just left. So the hook is grown — out and towards the compound —
    until the curve really does clear the node and stay clear.

    Growing one handle keeps the edge a single curve. Bending it around the
    node through waypoint markers would put the turn where it is wanted, but
    the corner it makes there is tight enough to read as a kink; one long
    sweep is both smoother and the same shape the rest of the map is drawn in.
    """
    (ax, ay), (bx, by) = p0, p3
    _, _, half_w, half_h = node
    reach = max(abs(bx - ax) * f1, 2 * half_w)
    tilt = (_sign(by - ay) or 1) * 2 * half_h
    for scale in growth:
        handle = (ax + side * reach * scale, ay + tilt * scale)
        if _leaves_by(p0, handle, handle2, p3, node, side):
            return handle
    log.warning("edge from (%.0f, %.0f) to (%.0f, %.0f) will not clear its "
                "member node on the %s — drawn with the widest hook",
                ax, ay, bx, by, "left" if side < 0 else "right")
    return handle


def _bezier(ax, ay, bx, by, style, exit_side=0, node=None):
    """Bezier handles for the segment from (ax, ay) to (bx, by).

    ``style.bezier_fracs`` places both handles along the run in x; how far
    each one has risen towards the compound is what ``style.edge_curve``
    decides — see :class:`MapStyle`.

    ``exit_side`` is the face of the member node this edge has to leave by:
    -1 for a compound the member consumes, +1 for one it excretes. Together
    with ``node`` — the member's ``(x, y, half width, half height)`` as it
    will be drawn — it keeps consumption arriving on the left of the box and
    excretion leaving on the right. Only a compound on the far side is drawn
    any differently, and only in its first handle (see :func:`_hook`); every
    other edge on the map is the plain curve it always was.
    """
    f1, f2 = style.bezier_fracs
    dx, dy = bx - ax, by - ay
    if style.edge_curve == "chord":
        rise1, rise2 = f1, f2
    elif style.edge_curve == "s":
        blend = _chord_blend(dx, dy, style.curve_steepness)
        rise1, rise2 = f1 * blend, 1.0 - (1.0 - f2) * blend
    else:
        rise1, rise2 = 0.0, style.arc_arrival
    handle2 = (ax + dx * f2, ay + dy * rise2)
    if exit_side and node and dx * exit_side < 0:
        handle1 = _hook((ax, ay), (bx, by), handle2, exit_side, node, f1)
    else:
        handle1 = (ax + dx * f1, ay + dy * rise1)
    return ({"x": handle1[0], "y": handle1[1]},
            {"x": handle2[0], "y": handle2[1]})


class _Counter:
    """Monotonic integer-string key allocator (Escher keys are '0', '1', …)."""

    def __init__(self):
        self.value = 0

    def next(self):
        key = str(self.value)
        self.value += 1
        return key


def _pack(desired, min_spacing):
    """Place items at their preferred heights, pushed apart to ``min_spacing``.

    ``desired`` maps item -> preferred y. Items keep the order of their
    preferences; the packed run is then re-centred on the mean preference so
    the column does not drift downwards.
    """
    if not desired:
        return {}
    order = sorted(desired, key=lambda key: (desired[key], str(key)))
    ys = []
    for key in order:
        y = desired[key]
        if ys and y < ys[-1] + min_spacing:
            y = ys[-1] + min_spacing
        ys.append(y)
    shift = sum(desired.values()) / len(desired) - sum(ys) / len(ys)
    return {key: y + shift for key, y in zip(order, ys)}


def _organ_pipe(items, weight):
    """Order items with the heaviest in the middle and the lightest at the ends.

    The ranked items are dealt alternately to the two sides of the centre, so
    the run reads 5th, 3rd, 1st, 2nd, 4th heaviest from one end to the other.
    Items of equal weight keep the order they arrived in.
    """
    ranked = sorted(items, key=weight, reverse=True)
    above, below = [], []
    for rank, item in enumerate(ranked):
        (below if rank % 2 else above).append(item)
    return list(reversed(above)) + below


def _column_halves(compounds, members_at, middle):
    """Split a column's compounds into ``(top, bottom)`` by where their members
    are.

    A compound goes with the half of the member column holding most of the
    members it touches. An even split is settled on their mean height, and a
    compound whose members balance exactly around the middle — one that only
    the centre member touches, most often — has no preference either way, so
    it is dealt to the emptier half and the column stays balanced.
    """
    top, bottom, undecided = [], [], []
    for compound in compounds:
        heights = members_at[compound]
        above = sum(1 for y in heights if y < middle)
        below = sum(1 for y in heights if y > middle)
        mean = sum(heights) / len(heights)
        if above != below:
            (top if above > below else bottom).append(compound)
        elif mean != middle:
            (top if mean < middle else bottom).append(compound)
        else:
            undecided.append(compound)
    for compound in undecided:
        (top if len(top) <= len(bottom) else bottom).append(compound)
    return top, bottom


def _connectivity_column(compounds, members_at, top, bottom,
                         min_spacing, span):
    """Place one compound column by how many members each compound touches.

    The compounds with the most edges sit in the middle of the column and the
    one-off compounds at its two ends. That puts them level with the members
    that have the most edges — who are in the middle of the member column for
    the same reason — so the busiest traffic crosses the middle of the map on
    short edges and the long thin ones are left to the ends, where there is
    room for them to cross nothing.

    Which end a compound works out from is decided by its members: one whose
    members lie mostly in the top half of the member column goes in the top
    half of the column, and likewise for the bottom (see
    :func:`_column_halves`), so a compound still sits on the side its edges
    come from. Within a half the compounds run from the busiest, next to the
    middle, out to the quietest at the end; compounds with the same number of
    edges are ordered by which sits closer to the middle already.

    ``members_at`` maps each compound to the heights of the members it
    touches — how many of them there are is its connectivity, and where they
    are is its half.

    With ``span`` each half fills its half of the member column; without it
    the halves pack together at ``min_spacing`` around the middle.
    """
    middle = (top + bottom) / 2
    placed = {}
    for half, direction in zip(_column_halves(compounds, members_at, middle),
                               (-1, 1)):
        if not half:
            continue
        order = sorted(half, key=lambda c: (
            -len(members_at[c]),
            direction * sum(members_at[c]) / len(members_at[c]),
            c))
        room = (bottom - top) / 2 if span else 0.0
        step = max(room / len(order), min_spacing)
        for index, compound in enumerate(order):
            placed[compound] = middle + direction * (index + 0.5) * step
    return placed


def _layout_block(members, style, compound_names, node_ids, segment_ids,
                  columns, lane_dx):
    """Lay out one condition as a three-column block, in local coordinates.

    The member column runs down ``x = 0`` starting at ``y = 0``, ordered by
    connectivity rather than by name; compound columns may extend above and
    below it. ``columns`` is the
    ``(input_x, output_x)`` pair and ``lane_dx`` the exchange-lane offset,
    both shared by every block so the map's columns line up. Returns
    ``(nodes, reactions)`` with ``reactions`` as a list in member order.
    """
    input_x, output_x = columns
    # The busiest members go in the middle of the column, which is where the
    # busiest compounds end up too, so most of the traffic is short and
    # horizontal and the quiet members at the two ends have room to spread.
    member_y = {index: slot * style.member_pitch for slot, index in enumerate(
        _organ_pipe(range(len(members)),
                    lambda index: len(members[index]["fluxes"])))}
    inputs, outputs, exchanged = classify_compounds(members)

    # Where a compound's members lie, and how many of them there are: the
    # first decides which half of its column it goes in, the second how close
    # to the middle of that column it sits.
    users = {}
    for index, member in enumerate(members):
        for compound in member["fluxes"]:
            users.setdefault(compound, []).append(member_y[index])
    desired = {c: sum(v) / len(v) for c, v in users.items()}

    # An exchanged compound goes in the lane holding most of its edges:
    # consumers attach to a member's left marker and producers to its right,
    # so the majority side is the one whose edges avoid crossing the column.
    # A tie carries no such preference, so it goes to the emptier lane; that
    # keeps the two lanes balanced instead of piling everything on the left.
    left_lane, right_lane = [], []
    for compound in exchanged:
        consumers = sum(1 for m in members if m["fluxes"].get(compound, 0) < 0)
        producers = sum(1 for m in members if m["fluxes"].get(compound, 0) > 0)
        if consumers == producers:
            go_left = len(left_lane) <= len(right_lane)
        else:
            go_left = consumers > producers
        (left_lane if go_left else right_lane).append(compound)

    top, bottom = 0.0, (len(members) - 1) * style.member_pitch
    nodes, compound_node = {}, {}
    for compounds, x, span in ((inputs, input_x, True),
                               (outputs, output_x, True),
                               (left_lane, -lane_dx, False),
                               (right_lane, lane_dx, False)):
        if span:
            column = _connectivity_column(
                compounds, users, top, bottom,
                style.min_node_spacing, style.span_columns)
        else:
            # the exchange lanes stay level with their members: compress the
            # spacing whenever the packed lane would exceed lane_span_fraction
            # of the member column's height (_pack re-centres on the
            # barycentre, so the compressed lane stays in the middle)
            wanted = {c: desired[c] for c in compounds}
            lane_spacing = style.lane_node_spacing
            if (style.lane_span_fraction is not None and len(wanted) > 1
                    and bottom > top):
                max_span = style.lane_span_fraction * (bottom - top)
                lane_spacing = min(lane_spacing, max_span / (len(wanted) - 1))
            column = _pack(wanted, lane_spacing)
        for compound, y in column.items():
            # Escher renders bigg_id as the on-map label (``name`` is only
            # tooltip metadata), so labels left of the axis are shifted by the
            # width of the compound ID.
            if x < 0:
                label_x = x - style.label_pad - _label_width(
                    compound, style.label_font_px)
            else:
                label_x = x + style.label_pad
            node_id = node_ids.next()
            compound_node[compound] = node_id
            nodes[node_id] = {
                "node_type": "metabolite",
                "x": x,
                "y": y,
                "bigg_id": compound,
                "name": compound_names.get(compound, compound),
                "label_x": label_x,
                "label_y": y + style.label_font_px / 4,
                "node_is_primary": True,
            }

    reactions = []
    for index, member in enumerate(members):
        cy = member_y[index]
        mid_id = node_ids.next()
        nodes[mid_id] = {"node_type": "midmarker", "x": 0.0, "y": cy}
        reactant_marker = node_ids.next()
        nodes[reactant_marker] = {"node_type": "multimarker",
                                  "x": -style.marker_offset, "y": cy}
        product_marker = node_ids.next()
        nodes[product_marker] = {"node_type": "multimarker",
                                 "x": style.marker_offset, "y": cy}

        half_width, half_height = style.member_node_half(member["name"])
        node = (0.0, cy, half_width, half_height)
        segments = {
            segment_ids.next(): {"from_node_id": reactant_marker,
                                 "to_node_id": mid_id, "b1": None, "b2": None},
            segment_ids.next(): {"from_node_id": product_marker,
                                 "to_node_id": mid_id, "b1": None, "b2": None},
        }
        metabolites = []
        for compound, flux in sorted(member["fluxes"].items()):
            marker = reactant_marker if flux < 0 else product_marker
            target = nodes[compound_node[compound]]
            b1, b2 = _bezier(nodes[marker]["x"], cy,
                             target["x"], target["y"], style,
                             exit_side=_sign(flux), node=node)
            segments[segment_ids.next()] = {
                "from_node_id": marker,
                "to_node_id": compound_node[compound],
                "b1": b1, "b2": b2,
            }
            metabolites.append({"bigg_id": compound, "coefficient": flux})

        reactions.append({
            "name": member["name"],
            "bigg_id": member["bigg_id"],
            "reversibility": False,
            # centred over the midmarker, clear of both exchange lanes
            "label_x": -_label_width(
                member["name"], style.reaction_label_font_px) / 2,
            "label_y": cy - style.member_label_gap,
            "gene_reaction_rule": "",
            "genes": [],
            "metabolites": metabolites,
            "segments": segments,
        })

    log.debug("block: %d members, %d inputs, %d outputs, %d exchanged",
              len(members), len(inputs), len(outputs), len(exchanged))
    return nodes, reactions


# ------------------------------------------------------------- assembly


def _block_extent(nodes, reactions, style, caption="", caption_x=0.0):
    """``(left, right, top, bottom)`` of one laid-out block.

    Labels count towards the extent, not just the nodes they belong to: a
    right-hand label reaches past its node by its own width, and the caption a
    block carries in a multi-block map sits ``block_label_gap`` above
    everything else. Node radii are included so a circle on the edge of the
    map is not sliced in half by the canvas.
    """
    lefts, rights, tops, bottoms = [], [], [], []

    def add(left, right, top, bottom):
        lefts.append(left)
        rights.append(right)
        tops.append(top)
        bottoms.append(bottom)

    for node in nodes.values():
        radius = (NODE_RADIUS if node.get("node_type") == "metabolite"
                  else MARKER_RADIUS)
        add(node["x"] - radius, node["x"] + radius,
            node["y"] - radius, node["y"] + radius)
        if "label_x" in node:
            width = _label_width(node.get("bigg_id", ""), style.label_font_px)
            add(node["label_x"], node["label_x"] + width,
                node["label_y"] - style.label_font_px, node["label_y"])
    for reaction in reactions:
        width = _label_width(reaction.get("name", ""),
                             style.reaction_label_font_px)
        add(reaction["label_x"], reaction["label_x"] + width,
            reaction["label_y"] - style.reaction_label_font_px,
            reaction["label_y"])
    if caption:
        caption_y = _caption_y(nodes, style)
        add(caption_x,
            caption_x + _label_width(caption, style.block_label_font_px),
            caption_y - style.block_label_font_px, caption_y)
    return min(lefts), max(rights), min(tops), max(bottoms)


def _caption_y(nodes, style):
    """Height of a block's caption — ``block_label_gap`` above its top node."""
    return min(node["y"] for node in nodes.values()) - style.block_label_gap


def _tile(extents, n_cols, style):
    """Place blocks in a grid ``n_cols`` across, row by row.

    A grid column is as wide as its blocks' bounding boxes together, and each
    of them is placed against the column's left edge rather than its own, so
    their member axes stay on one vertical line; a row is likewise as tall as
    its blocks together, aligned on the top of the row.

    Returns ``(offsets, width, height)`` — the ``(dx, dy)`` to move each block
    by, in block order, and the size of the grid they end up filling with its
    top-left corner at the origin.
    """
    rows = [extents[i:i + n_cols] for i in range(0, len(extents), n_cols)]

    x_offsets, cursor = [], 0.0
    for column in range(n_cols):
        cells = [row[column] for row in rows if column < len(row)]
        left = min(left for left, _, _, _ in cells)
        right = max(right for _, right, _, _ in cells)
        x_offsets.append(cursor - left)
        cursor += (right - left) + style.block_column_gap
    width = max(cursor - style.block_column_gap, 0.0)

    offsets, cursor = [], 0.0
    for row in rows:
        top = min(top for _, _, top, _ in row)
        bottom = max(bottom for _, _, _, bottom in row)
        for column in range(len(row)):
            offsets.append((x_offsets[column], cursor - top))
        cursor += (bottom - top) + style.block_gap
    height = max(cursor - style.block_gap, 0.0)
    return offsets, width, height


def _canvas_size(width, height, style):
    """The exported canvas — the tiled grid plus its margin on every side."""
    return (width + 2 * style.canvas_margin, height + 2 * style.canvas_margin)


def _within_aspect(width, height, max_aspect):
    """Is ``width:height`` inside ``1:max_aspect`` .. ``max_aspect:1``?

    The corrections aim at the limit exactly, so the comparison is made with a
    hair of tolerance — landing a rounding error short of the target is not
    worth another pass over the whole map.
    """
    if max_aspect is None or not width or not height:
        return True
    ratio = width / height
    return ((1.0 / max_aspect) * (1 - ASPECT_TOLERANCE) <= ratio
            <= max_aspect * (1 + ASPECT_TOLERANCE))


def _choose_n_cols(extents, style):
    """How many blocks to put in a row.

    Blocks per row is the one control that trades the map's height for its
    width, so it is chosen first: width grows and height falls as it rises,
    which makes the fewest blocks per row that land inside ``max_aspect`` also
    the tallest, most column-like arrangement that qualifies. When no count
    fits — a single very tall block, most often — the closest one is used and
    :func:`build_escher_map` widens the blocks themselves from there.
    """
    if style.max_aspect is None or len(extents) < 2:
        return 1

    def _miss(n_cols):
        _, grid_width, grid_height = _tile(extents, n_cols, style)
        width, height = _canvas_size(grid_width, grid_height, style)
        ratio = math.log(width / height) if width and height else 0.0
        limit = math.log(style.max_aspect)
        return max(0.0, -limit - ratio, ratio - limit)

    best, best_miss = 1, None
    for n_cols in range(1, len(extents) + 1):
        miss = _miss(n_cols)
        if best_miss is None or miss < best_miss:
            best, best_miss = n_cols, miss
        if not miss:
            break
    return best


def _shift_block(nodes, reactions, dx, dy):
    """Move a laid-out block by ``(dx, dy)``.

    Bezier handles hold absolute coordinates, so they have to travel with the
    segment they shape — a handle left behind drags its edge back towards
    wherever the block was laid out.
    """
    for node in nodes.values():
        node["x"] += dx
        node["y"] += dy
        if "label_x" in node:
            node["label_x"] += dx
            node["label_y"] += dy
    for reaction in reactions:
        reaction["label_x"] += dx
        reaction["label_y"] += dy
        for segment in reaction["segments"].values():
            for handle in (segment["b1"], segment["b2"]):
                if handle:
                    handle["x"] += dx
                    handle["y"] += dy


def build_escher_map(blocks, compound_names=None,
                     map_name="community_exchange_map",
                     map_description=None, style=None):
    """Assemble a complete Escher map from one or more blocks of members.

    Blocks are tiled into a grid and, if that is not enough, the lanes and
    the input/output columns are pushed outwards, until the canvas is inside
    ``style.max_aspect`` (see :class:`MapStyle`) — so the figure keeps a
    usable shape however large the community is. Every block is laid out
    against the same lane and column positions, and blocks sharing a grid
    column are placed against the same left edge, so their member axes line
    up down the map.

    Parameters
    ----------
    blocks : list[tuple[str, list[dict]]] or list[dict]
        Either ``[(block_label, members), ...]`` — one block per condition,
        each captioned with a ``text_labels`` entry in the style of the
        reference map's "Days 5 to 7" captions — or a bare member list for a
        single unlabelled block. Compounds are collapsed within a block,
        never across blocks, so each condition keeps its own node set.
    compound_names : dict[str, str] or None
        Compound ID -> display name for metabolite node ``name`` fields.
    map_name, map_description : str
        Header metadata.
    style : MapStyle or None
        Layout geometry; defaults to ``MapStyle()``.

    Returns
    -------
    list
        The two-element ``[header, body]`` Escher map structure.
    """
    if blocks and isinstance(blocks[0], dict):
        blocks = [("", blocks)]
    style = style or MapStyle()
    compound_names = compound_names or {}

    drawable = []
    for block_label, members in blocks:
        if members:
            drawable.append((block_label, members))
        else:
            log.warning("block %r has no members; skipped", block_label)
    if not drawable:
        raise ValueError("no members to draw")

    # One geometry for the whole map — lanes fitted to its tallest member
    # column, columns to the widest compound label anywhere in it — so the
    # blocks stay aligned with each other and read at the same scale.
    lane_dx = style.fitted_lane_dx(
        max((len(members) - 1) * style.member_pitch for _, members in drawable))
    columns = style.column_positions(
        {c for _, members in drawable for m in members for c in m["fluxes"]},
        lane_dx)

    def lay_out(style, columns, lane_dx):
        """Lay every block out in its own local coordinates."""
        node_ids, segment_ids = _Counter(), _Counter()
        laid = [(label,) + _layout_block(members, style, compound_names,
                                         node_ids, segment_ids, columns,
                                         lane_dx)
                for label, members in drawable]
        extents = [_block_extent(nodes, reactions, style, label, columns[0])
                   for label, nodes, reactions in laid]
        return laid, extents, segment_ids.value

    laid, extents, n_segments = lay_out(style, columns, lane_dx)
    n_cols = _choose_n_cols(extents, style)
    offsets, width, height = _tile(extents, n_cols, style)

    # Tiling alone cannot always reach the target — one condition with forty
    # members is a ribbon whatever else shares its row — so the blocks
    # themselves take up the rest. Either correction changes the extents
    # underneath the grid, hence the re-layout and the further passes.
    for _ in range(ASPECT_PASSES):
        canvas_width, canvas_height = _canvas_size(width, height, style)
        if _within_aspect(canvas_width, canvas_height, style.max_aspect):
            break
        ratio = canvas_width / canvas_height
        if ratio * style.max_aspect < 1.0:
            # Too tall. Moving the lanes and the columns out by the same delta
            # widens every block by exactly 2 * delta whichever of them is
            # outermost — a block whose compounds are all cross-fed has no
            # column nodes to move — and leaves the gap between them, and so
            # the room its labels need, exactly as it was. No height changes,
            # so this converges in one pass.
            delta = (canvas_height / style.max_aspect - canvas_width) / (2 * n_cols)
            columns = (columns[0] - delta, columns[1] + delta)
            lane_dx += delta
            log.info("aspect %.3f is narrower than 1:%g — moving the columns "
                     "and lanes %.0f px further out",
                     ratio, style.max_aspect, delta)
        else:
            # too wide: gain the height by opening the vertical spacings, so
            # it carries nodes rather than blank canvas
            factor = max(1.0, (canvas_width / style.max_aspect
                               - 2 * style.canvas_margin) / max(height, 1.0))
            log.info("aspect %.3f is wider than %g:1 — opening the vertical "
                     "spacing %.2fx", ratio, style.max_aspect, factor)
            style = style.vertically_scaled(factor)
        laid, extents, n_segments = lay_out(style, columns, lane_dx)
        n_cols = _choose_n_cols(extents, style)
        before = (width, height)
        offsets, width, height = _tile(extents, n_cols, style)
        if (width, height) == before:
            # the correction had nothing to bite on — a block with a single
            # member and a single compound per side has no spacing to open up
            break

    nodes, reactions, text_labels = {}, {}, {}
    reaction_ids, label_ids = _Counter(), _Counter()
    for (block_label, block_nodes, block_reactions), (dx, dy) in zip(laid, offsets):
        _shift_block(block_nodes, block_reactions, dx, dy)
        if block_label:
            text_labels[label_ids.next()] = {
                "x": columns[0] + dx,
                "y": _caption_y(block_nodes, style),
                "text": block_label,
            }
        nodes.update(block_nodes)
        for reaction in block_reactions:
            reactions[reaction_ids.next()] = reaction

    # _tile puts the grid's top-left corner at the origin, extents and all.
    canvas_width, canvas_height = _canvas_size(width, height, style)
    canvas_x = canvas_y = -style.canvas_margin
    if style.max_aspect is not None:
        # Whatever the layout could not deliver is made up with blank canvas.
        # Usually that is the last fraction of a pixel — the corrections
        # converge on the limit from outside, and finishing the job here is
        # cheaper than another pass over the map. A real shortfall means the
        # layout would not stretch at all, which is worth saying out loud: the
        # figure then keeps the promise the ratio makes about the canvas, but
        # the map inside it is the shape it always was.
        short_width = canvas_height / style.max_aspect - canvas_width
        short_height = canvas_width / style.max_aspect - canvas_height
        padded = 0.0
        if short_width > 0:
            canvas_x -= short_width / 2
            canvas_width += short_width
            padded = short_width / canvas_width
        elif short_height > 0:
            canvas_y -= short_height / 2
            canvas_height += short_height
            padded = short_height / canvas_height
        if padded > ASPECT_PAD_WARNING:
            log.warning("the layout would not stretch to %g:1 — padding the "
                        "canvas out to %.0f x %.0f with blank space instead",
                        style.max_aspect, canvas_width, canvas_height)
    canvas = {
        "x": canvas_x,
        "y": canvas_y,
        "width": canvas_width,
        "height": canvas_height,
    }

    log.info("built map: %d reactions, %d nodes, %d segments, "
             "%d block(s) %d across, canvas %.0f x %.0f (%.2f:1)",
             len(reactions), len(nodes), n_segments, len(laid), n_cols,
             canvas_width, canvas_height, canvas_width / canvas_height)
    return [
        {
            "map_name": map_name,
            "map_id": uuid.uuid4().hex[:12],
            "map_description": map_description or
            f"Net organism exchange reactions for {len(reactions)} members",
            "homepage": "https://escher.github.io",
            "schema": ESCHER_SCHEMA,
        },
        {
            "reactions": reactions,
            "nodes": nodes,
            "text_labels": text_labels,
            "canvas": canvas,
        },
    ]


# ------------------------------------------------------------ top level


def build_map_from_interactions(csv_path, output_path=None, names_csv=None,
                                conditions=None, min_abs_flux=0.0,
                                skip_names=None, skip_ids=None,
                                map_name=None, style=None,
                                separate_maps=False):
    """Read ``ASVMetaboliteInteractions.csv`` and write Escher map(s).

    Parameters
    ----------
    csv_path : str or Path
        The wide interaction matrix (compound-ID columns).
    output_path : str or Path or None
        Output JSON. Defaults to ``<stem>_membermap.json`` beside the input;
        with ``separate_maps`` the condition is inserted before the suffix.
    names_csv : str or Path or None
        The name-headed twin of ``csv_path``, used for display names and for
        matching ``skip_names``.
    conditions : iterable[str] or None
        Condition (model) ids to include; default is every condition in the
        file except :data:`AVERAGE_CONDITION`. Note that a map built from
        ``"ave"`` has reaction ``bigg_id``s without the ``<diet>-ABX_<day>``
        suffix, so ``model_mapping.build_direction_tracking`` cannot parse it.
    separate_maps : bool
        Write one map per condition instead of stacking the conditions as
        captioned blocks in a single map.

    Returns
    -------
    pathlib.Path or list[pathlib.Path]
        The path(s) written.
    """
    csv_path = Path(csv_path)
    data = parse_interaction_matrix(csv_path)
    compound_names = load_compound_names(csv_path, names_csv) if names_csv else {}

    if conditions is None:
        conditions = [c for c in data if c != AVERAGE_CONDITION]
    else:
        unknown = [c for c in conditions if c not in data]
        if unknown:
            raise ValueError(
                f"unknown condition(s) {unknown}; {csv_path.name} has "
                f"{sorted(data)}")
        conditions = list(conditions)
    if not conditions:
        raise ValueError(f"no conditions to build from {csv_path.name}")

    def _members(condition):
        return build_member_reactions(
            data[condition], model_id=condition,
            compound_names=compound_names, min_abs_flux=min_abs_flux,
            skip_names=skip_names, skip_ids=skip_ids,
        )

    def _write(escher_map, path):
        path = Path(path)
        with path.open("w") as handle:
            dump(escher_map, handle, indent=3)
        log.info("wrote Escher map to %s", path)
        return path

    default_out = output_path or csv_path.with_name(csv_path.stem + "_membermap.json")
    default_out = Path(default_out)

    if separate_maps:
        written = []
        for condition in conditions:
            escher_map = build_escher_map(
                _members(condition), compound_names=compound_names,
                map_name=map_name or f"{csv_path.stem}_{condition}",
                style=style,
            )
            safe = condition.replace("/", "-")
            written.append(_write(escher_map, default_out.with_name(
                f"{default_out.stem}_{safe}{default_out.suffix}")))
        return written

    blocks = [(condition, _members(condition)) for condition in conditions]
    escher_map = build_escher_map(
        blocks, compound_names=compound_names,
        map_name=map_name or csv_path.stem, style=style,
    )
    return _write(escher_map, default_out)


def main():
    # The CLI takes its defaults from MapStyle rather than repeating them, so
    # the two cannot drift apart; None means "whatever MapStyle says".
    defaults = MapStyle()
    parser = argparse.ArgumentParser(
        description="Build an Escher map of per-member net exchange "
                    "reactions from community flux data.")
    parser.add_argument(
        "csv_path",
        help="ASVMetaboliteInteractions.csv (compound-ID column headers)")
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output JSON path (default: <stem>_membermap.json)")
    parser.add_argument(
        "-n", "--names", default=None,
        help="ASVMetaboliteInteractions_names.csv, for compound display names")
    parser.add_argument(
        "-c", "--condition", action="append", dest="conditions", default=None,
        help="Condition/model id to include; repeatable "
             "(default: every condition except 'ave')")
    parser.add_argument(
        "--separate-maps", action="store_true",
        help="Write one map per condition instead of one map of captioned "
             "blocks tiled into a grid")
    parser.add_argument(
        "--min-flux", type=float, default=0.0,
        help="Drop exchanges with abs(flux) below this (reference used 0.05)")
    parser.add_argument(
        "--skip-amino-acids", action="store_true",
        help="Exclude the amino acids and extras in config/filter.json, "
             "as the reference figure does (requires --names)")
    parser.add_argument(
        "--member-pitch", type=float, default=None,
        help=f"Vertical spacing between member reactions "
             f"(default: {defaults.member_pitch:g})")
    parser.add_argument(
        "--column-dx", type=float, default=None,
        help="Horizontal distance from the member column to the input and "
             "output columns (default: fitted just clear of the lane labels)")
    parser.add_argument(
        "--lane-dx", type=float, default=None,
        help=f"Smallest horizontal distance from the member column to the "
             f"exchanged-compound lanes (default: {defaults.mixed_lane_dx:g})")
    parser.add_argument(
        "--lane-dx-fraction", type=float, default=None,
        help=f"Grow that distance to this fraction of the tallest member "
             f"column, so lane edges stay slanted in a big community "
             f"(default: {defaults.lane_dx_fraction:g}; 0 stops the growth, "
             f"though the aspect fit may still move the lanes)")
    parser.add_argument(
        "--node-spacing", type=float, default=None,
        help=f"Minimum vertical spacing between nodes in the input/output "
             f"columns (default: {defaults.min_node_spacing:g})")
    parser.add_argument(
        "--lane-spacing", type=float, default=None,
        help=f"Minimum vertical spacing between nodes in the exchange lanes "
             f"(default: {defaults.lane_node_spacing:g})")
    parser.add_argument(
        "--edge-curve", choices=("arc", "s", "chord"), default=None,
        help=f"'arc' bends each edge once, level out of its marker and into "
             f"its compound at an angle; 's' bends it twice, level at both "
             f"ends; 'chord' draws it flat. None of them changes which face "
             f"of a member an edge uses, nor the detour a compound on the far "
             f"side needs (default: {defaults.edge_curve})")
    parser.add_argument(
        "--arc-arrival", type=float, default=None,
        help=f"How far up to the compound the second handle of an arc sits: "
             f"lower bows the arc more (default: {defaults.arc_arrival:g}, "
             f"ceiling {_single_bend_limit(defaults.bezier_fracs):.3g})")
    parser.add_argument(
        "--max-aspect", type=float, default=None,
        help=f"Keep the canvas within this width:height ratio either way, by "
             f"tiling the blocks wider and moving the lanes and columns out — "
             f"which overrides --column-dx and --lane-dx if it has to "
             f"(default: {defaults.max_aspect:g})")
    parser.add_argument(
        "--no-aspect-limit", action="store_true",
        help="Stack the blocks in one column and let the figure end up "
             "whatever shape it wants")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    skip_names = None
    if args.skip_amino_acids:
        from .filter_map import DEFAULT_SKIP_NAMES
        skip_names = DEFAULT_SKIP_NAMES
        if not args.names:
            log.warning("--skip-amino-acids matches display names; "
                        "pass --names or nothing will be skipped")

    geometry = {
        "member_pitch": args.member_pitch,
        "mixed_lane_dx": args.lane_dx,
        "lane_dx_fraction": args.lane_dx_fraction,
        "min_node_spacing": args.node_spacing,
        "lane_node_spacing": args.lane_spacing,
        "edge_curve": args.edge_curve,
        "arc_arrival": args.arc_arrival,
        "max_aspect": args.max_aspect,
    }
    geometry = {name: value for name, value in geometry.items()
                if value is not None}
    if args.column_dx is not None:
        geometry["input_column_dx"] = -abs(args.column_dx)
        geometry["output_column_dx"] = abs(args.column_dx)
    if args.no_aspect_limit:
        geometry["max_aspect"] = None

    build_map_from_interactions(
        args.csv_path,
        output_path=args.output,
        names_csv=args.names,
        conditions=args.conditions,
        min_abs_flux=args.min_flux,
        skip_names=skip_names,
        style=MapStyle(**geometry),
        separate_maps=args.separate_maps,
    )


if __name__ == "__main__":
    main()
