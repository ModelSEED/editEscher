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
  ones — sit in two narrow lanes flanking the member column, level with the
  members that exchange them.

Within each region a node is placed at the mean height of the members it
connects to, so edges stay short and cross-feeding reads off the middle of
the map.

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
import csv
import logging
import uuid
from json import dump
from pathlib import Path

log = logging.getLogger(__name__)

ESCHER_SCHEMA = "https://escher.github.io/escher/jsonschema/1-0-0#"

#: Sentinel condition whose reaction ``bigg_id``s do not carry the
#: ``<diet>-ABX_<day>`` suffix, so ``model_mapping.build_direction_tracking``
#: cannot parse a map built from it.
AVERAGE_CONDITION = "ave"


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
      narrow lanes flanking the member column at ``±mixed_lane_dx``, level
      with the members that consume and produce them.

    Within each region a node sits at the mean height of the members it links
    to, pushed apart as far as ``min_node_spacing`` (outer columns) or
    ``lane_node_spacing`` (exchange lanes) requires. The lanes get the wider
    default because they hold the busiest compounds, whose barycentres all
    fall near the middle of the member column. ``lane_span_fraction`` keeps a
    lane from growing into a fourth full-height column: when the packed lane
    would exceed that fraction of the member column's height, its spacing is
    compressed to fit, still centred on the lane's barycentre.

    ``span_columns`` makes the input and output columns run the full height of
    the member column — they read as columns rather than as a knot in the
    middle — while the exchange lanes always stay level with their members,
    which is the point of drawing them there.

    ``input_column_dx`` and ``output_column_dx`` default to ``None``, which
    fits the columns as close to the member column as the exchange-lane
    labels allow (see :meth:`fitted_column_dx`). Pass numbers to place them
    by hand.
    """

    def __init__(self, member_pitch=420.0, input_column_dx=None,
                 output_column_dx=None, mixed_lane_dx=240.0,
                 column_clearance=160.0,
                 min_node_spacing=90.0, lane_node_spacing=200.0,
                 lane_span_fraction=0.6,
                 span_columns=True, marker_offset=20.0,
                 member_label_gap=60.0, label_pad=22.0, label_font_px=20.0,
                 reaction_label_font_px=30.0, bezier_fracs=(0.25, 0.75),
                 block_gap=900.0, block_label_gap=280.0, canvas_margin=400.0):
        self.member_pitch = member_pitch
        self.span_columns = span_columns
        self.input_column_dx = input_column_dx
        self.output_column_dx = output_column_dx
        self.mixed_lane_dx = mixed_lane_dx
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
        self.block_gap = block_gap
        self.block_label_gap = block_label_gap
        self.canvas_margin = canvas_margin

    def fitted_column_dx(self, compounds):
        """Tightest offset for the input/output columns, in px.

        An exchange-lane node's label runs outwards from the lane, so the
        outer columns can come in no closer than the far edge of those labels
        plus ``column_clearance``. Fitting to that keeps the figure compact
        without the lane labels colliding with the column nodes.
        """
        widest = max((len(compound) for compound in compounds), default=8)
        lane_label_edge = (self.mixed_lane_dx + self.label_pad
                           + 0.6 * self.label_font_px * widest)
        return lane_label_edge + self.column_clearance

    def column_positions(self, compounds):
        """``(input_x, output_x)`` — the fitted offsets, or the overrides."""
        fitted = self.fitted_column_dx(compounds)
        return (self.input_column_dx if self.input_column_dx is not None else -fitted,
                self.output_column_dx if self.output_column_dx is not None else fitted)


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


def _bezier(ax, ay, bx, by, fracs):
    """Bezier handles placed along the chord from (ax, ay) to (bx, by)."""
    f1, f2 = fracs
    return ({"x": ax + (bx - ax) * f1, "y": ay + (by - ay) * f1},
            {"x": ax + (bx - ax) * f2, "y": ay + (by - ay) * f2})


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


def _spread(desired, top, bottom, min_spacing):
    """Distribute items evenly between ``top`` and ``bottom``.

    Items keep the order of their preferred heights — that ordering is what
    keeps edges from crossing — but are spaced evenly so the column runs the
    full height of the member column instead of bunching at the barycentre.
    """
    if not desired:
        return {}
    order = sorted(desired, key=lambda key: (desired[key], str(key)))
    if len(order) == 1:
        return {order[0]: (top + bottom) / 2}
    step = max((bottom - top) / (len(order) - 1), min_spacing)
    start = (top + bottom) / 2 - step * (len(order) - 1) / 2
    return {key: start + index * step for index, key in enumerate(order)}


def _layout_block(members, style, compound_names, node_ids, segment_ids,
                  columns):
    """Lay out one condition as a three-column block, in local coordinates.

    The member column runs down ``x = 0`` starting at ``y = 0``; compound
    columns may extend above and below it. ``columns`` is the
    ``(input_x, output_x)`` pair, shared by every block so the columns line
    up. Returns ``(nodes, reactions)`` with ``reactions`` as a list in member
    order.
    """
    input_x, output_x = columns
    member_y = {index: index * style.member_pitch
                for index in range(len(members))}
    inputs, outputs, exchanged = classify_compounds(members)

    # Barycentre: a compound wants to sit level with the members it links to.
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

    top, bottom = 0.0, member_y[len(members) - 1]
    nodes, compound_node = {}, {}
    for compounds, x, span in ((inputs, input_x, True),
                               (outputs, output_x, True),
                               (left_lane, -style.mixed_lane_dx, False),
                               (right_lane, style.mixed_lane_dx, False)):
        wanted = {c: desired[c] for c in compounds}
        if span and style.span_columns:
            column = _spread(wanted, top, bottom, style.min_node_spacing)
        elif span:
            column = _pack(wanted, style.min_node_spacing)
        else:
            # the exchange lanes stay level with their members: compress the
            # spacing whenever the packed lane would exceed lane_span_fraction
            # of the member column's height (_pack re-centres on the
            # barycentre, so the compressed lane stays in the middle)
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
                label_x = x - style.label_pad - 0.6 * style.label_font_px * len(compound)
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
            b1, b2 = _bezier(nodes[marker]["x"], cy, target["x"], target["y"],
                             style.bezier_fracs)
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
            "label_x": -0.3 * style.reaction_label_font_px * len(member["name"]),
            "label_y": cy - style.member_label_gap,
            "gene_reaction_rule": "",
            "genes": [],
            "metabolites": metabolites,
            "segments": segments,
        })

    log.debug("block: %d members, %d inputs, %d outputs, %d exchanged",
              len(members), len(inputs), len(outputs), len(exchanged))
    return nodes, reactions


def build_escher_map(blocks, compound_names=None,
                     map_name="community_exchange_map",
                     map_description=None, style=None):
    """Assemble a complete Escher map from one or more blocks of members.

    Parameters
    ----------
    blocks : list[tuple[str, list[dict]]] or list[dict]
        Either ``[(block_label, members), ...]`` — one block per condition,
        stacked vertically and captioned with a ``text_labels`` entry, in
        the style of the reference map's "Days 5 to 7" captions — or a bare
        member list for a single unlabelled block. Compounds are collapsed
        within a block, never across blocks, so each condition keeps its own
        node set.
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

    if not any(members for _, members in blocks):
        raise ValueError("no members to draw")

    nodes, reactions, text_labels = {}, {}, {}
    node_ids, segment_ids = _Counter(), _Counter()
    reaction_ids, label_ids = _Counter(), _Counter()

    # One column position for the whole map, fitted to the widest compound
    # label anywhere in it, so blocks stay aligned with each other.
    columns = style.column_positions(
        {c for _, members in blocks for m in members for c in m["fluxes"]})

    cursor = 0.0
    for block_label, members in blocks:
        if not members:
            log.warning("block %r has no members; skipped", block_label)
            continue
        block_nodes, block_reactions = _layout_block(
            members, style, compound_names, node_ids, segment_ids, columns)

        top = min(node["y"] for node in block_nodes.values())
        bottom = max(node["y"] for node in block_nodes.values())
        shift = cursor - top
        for node in block_nodes.values():
            node["y"] += shift
            if "label_y" in node:
                node["label_y"] += shift
        for reaction in block_reactions:
            reaction["label_y"] += shift

        if block_label:
            text_labels[label_ids.next()] = {
                "x": columns[0],
                "y": cursor - style.block_label_gap,
                "text": block_label,
            }

        nodes.update(block_nodes)
        for reaction in block_reactions:
            reactions[reaction_ids.next()] = reaction
        cursor += (bottom - top) + style.block_gap

    xs = [n["x"] for n in nodes.values()]
    ys = [n["y"] for n in nodes.values()]
    xs += [n["label_x"] for n in nodes.values() if "label_x" in n]
    ys += [n["label_y"] for n in nodes.values() if "label_y" in n]
    pad = style.canvas_margin
    canvas = {
        "x": min(xs) - pad,
        "y": min(ys) - pad,
        "width": (max(xs) - min(xs)) + 2 * pad,
        "height": (max(ys) - min(ys)) + 2 * pad,
    }

    log.info("built map: %d reactions, %d nodes, %d segments",
             len(reactions), len(nodes), segment_ids.value)
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
        help="Write one map per condition instead of stacked blocks")
    parser.add_argument(
        "--min-flux", type=float, default=0.0,
        help="Drop exchanges with abs(flux) below this (reference used 0.05)")
    parser.add_argument(
        "--skip-amino-acids", action="store_true",
        help="Exclude the amino acids and extras in config/filter.json, "
             "as the reference figure does (requires --names)")
    parser.add_argument(
        "--member-pitch", type=float, default=420.0,
        help="Vertical spacing between member reactions")
    parser.add_argument(
        "--column-dx", type=float, default=None,
        help="Horizontal distance from the member column to the input and "
             "output columns (default: fitted just clear of the lane labels)")
    parser.add_argument(
        "--lane-dx", type=float, default=300.0,
        help="Horizontal distance from the member column to the exchanged-"
             "compound lanes")
    parser.add_argument(
        "--node-spacing", type=float, default=90.0,
        help="Minimum vertical spacing between nodes in the input/output columns")
    parser.add_argument(
        "--lane-spacing", type=float, default=200.0,
        help="Minimum vertical spacing between nodes in the exchange lanes")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    skip_names = None
    if args.skip_amino_acids:
        from .filter_map import DEFAULT_SKIP_NAMES
        skip_names = DEFAULT_SKIP_NAMES
        if not args.names:
            log.warning("--skip-amino-acids matches display names; "
                        "pass --names or nothing will be skipped")

    build_map_from_interactions(
        args.csv_path,
        output_path=args.output,
        names_csv=args.names,
        conditions=args.conditions,
        min_abs_flux=args.min_flux,
        skip_names=skip_names,
        style=MapStyle(
            member_pitch=args.member_pitch,
            input_column_dx=-args.column_dx if args.column_dx else None,
            output_column_dx=args.column_dx,
            mixed_lane_dx=args.lane_dx,
            min_node_spacing=args.node_spacing,
            lane_node_spacing=args.lane_spacing,
        ),
        separate_maps=args.separate_maps,
    )


if __name__ == "__main__":
    main()
