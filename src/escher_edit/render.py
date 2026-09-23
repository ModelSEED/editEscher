"""Render a community map to a finished SVG figure, without a browser.

Escher's JSON carries all of a map's geometry, so the SVG the Escher web
editor would export can be written directly: :func:`write_escher_svg` does
that, mirroring Escher's own drawing — the same stylesheet rules and
``#reactions`` / ``#nodes`` groups, ``r<key>`` / ``s<key>`` ids, segments
pulled back from their metabolites the way Escher's ``displacedCoords``
does, and arrowheads only where a product is made.
:func:`render_map_svg` then runs that through
:func:`escher_edit.svg_editor.EscherSVG_processing`, so a rendered figure and
an exported-then-processed one come out alike: every member a coloured box
node, with its edges in its colour and the cross-feeding ones dashed. Beside
the SVG goes an interactive HTML version (:mod:`escher_edit.interactive`),
where pointing at an edge, a compound or a member picks out everything it
touches.

The member boxes are the reason this exists. Escher's schema has no box
node, so a map opened in Escher itself — or re-saved from its editor, which
drops the box size the builder records — can only show a member as the
vertex its edges meet at. A figure generated here always has them.

    escher-edit-render map.json [-o figure.svg]
    escher-edit-render map.json --member-groups taxonomy.json --group-rank phylum
"""
import argparse
import json
import logging
import math
import shutil
import tempfile
from pathlib import Path
from xml.sax.saxutils import escape

from .build_map import MapStyle, cross_feeding_segments
from .interactive import write_interactive_html
from .palette import taxon_groups
from .svg_editor import EscherStyle, EscherSVG_processing, _load_map

log = logging.getLogger(__name__)

# The Escher stylesheet rules svg_editor rewrites, verbatim, so its string
# replacements find them; the canvas is a plain white ground here, where the
# editor draws it with a grey border.
CSS = (
    "svg.escher-svg #canvas{stroke:none;fill:#fff}"
    "svg.escher-svg .label{font-family:sans-serif;font-style:italic;"
    "font-weight:700;font-size:8px;fill:#000;stroke:none}"
    "svg.escher-svg .reaction-label{font-size:30px;fill:#202078;"
    "text-rendering:optimizelegibility}"
    "svg.escher-svg .node-label{font-size:20px}"
    "svg.escher-svg .text-label .label,svg.escher-svg .text-label-input"
    "{font-size:50px}"
    "svg.escher-svg .node-circle{stroke-width:2px}"
    "svg.escher-svg .midmarker-circle,svg.escher-svg .multimarker-circle"
    "{fill:#fff;fill-opacity:.2;stroke:#323232}"
    "svg.escher-svg .metabolite-circle{stroke:#a24510;fill:#e0865b}"
    "svg.escher-svg .segment{stroke:#334E75;stroke-width:10px;fill:none}"
    "svg.escher-svg .arrowhead{fill:#334E75}"
)

# Escher's defaults: arrowhead width and length, marker radius, and the
# daylight left between an arrowhead and its node. Metabolite radii come from
# the layout (MapStyle.node_radius).
ARROW_W, ARROW_H = 20.0, 13.0
MARKER_R = 5.0
ARROW_GAP = 10.0


def _displacement(coefficient, reversibility, radius):
    """Escher's ``get_disp``: how far a segment stops short of a metabolite —
    its radius, plus room for an arrowhead where one is drawn."""
    drawn = reversibility or coefficient > 0
    return radius + (ARROW_H if drawn else 0.0) + ARROW_GAP


def _displaced(start, end, length, which):
    """Escher's ``displacedCoords``: move ``length`` along ``start -> end``,
    from the start or back from the end."""
    hyp = math.dist(start, end)
    if not length or not hyp:
        log.warning("no room to pull a segment back from its node at %s", start)
        return start if which == "start" else end
    ux, uy = (end[0] - start[0]) / hyp, (end[1] - start[1]) / hyp
    if which == "start":
        return start[0] + length * ux, start[1] + length * uy
    return end[0] - length * ux, end[1] - length * uy


def write_escher_svg(escher_map, out_path, layout=None):
    """Write ``escher_map`` as the SVG Escher's editor would export for it.

    Unstyled beyond Escher's defaults — :func:`render_map_svg` is the finished
    figure. ``escher_map`` is a ``[header, body]`` map or a path to one.
    Metabolite nodes are drawn at the radii of ``layout``, the
    :class:`~escher_edit.build_map.MapStyle` the map was built with (a default
    one if not given), where Escher would use its own settings: the labels
    were placed clear of those circles. Returns the path written.
    """
    layout = layout or MapStyle()
    body = _load_map(escher_map)[1]
    nodes, canvas = body["nodes"], body["canvas"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" class="escher-svg" '
        f'viewBox="{canvas["x"]} {canvas["y"]} {canvas["width"]} {canvas["height"]}" '
        f'width="{canvas["width"]}" height="{canvas["height"]}">',
        f'<defs><style type="text/css">{CSS}</style></defs>',
        f'<rect id="canvas" x="{canvas["x"]}" y="{canvas["y"]}" '
        f'width="{canvas["width"]}" height="{canvas["height"]}"/>',
        '<g id="reactions">',
    ]
    for key, reaction in body["reactions"].items():
        parts.append(f'<g id="r{key}" class="reaction">')
        parts.append(
            f'<g class="reaction-label-group" transform="translate('
            f'{reaction.get("label_x", 0)},{reaction.get("label_y", 0)})">'
            f'<text class="reaction-label label">'
            f'{escape(str(reaction.get("bigg_id", "")))}</text></g>')
        coefficients = {met["bigg_id"]: met["coefficient"]
                        for met in reaction.get("metabolites", [])}
        reversible = reaction.get("reversibility", False)
        for segment_key, segment in reaction["segments"].items():
            a, b = nodes[segment["from_node_id"]], nodes[segment["to_node_id"]]
            start, end = (a["x"], a["y"]), (b["x"], b["y"])
            b1 = (segment["b1"]["x"], segment["b1"]["y"]) if segment.get("b1") else None
            b2 = (segment["b2"]["x"], segment["b2"]["y"]) if segment.get("b2") else None
            heads = []
            for node, which in ((a, "start"), (b, "end")):
                if node.get("node_type") != "metabolite":
                    continue
                coefficient = coefficients.get(node.get("bigg_id"), 0)
                length = _displacement(
                    coefficient, reversible,
                    layout.node_radius(node.get("node_is_primary", True)))
                if which == "start":
                    toward = b1 or end
                    start = _displaced(start, toward, length, "start")
                    tip = start
                else:
                    toward = b2 or start
                    end = _displaced(toward, end, length, "end")
                    tip = end
                if reversible or coefficient > 0:
                    heads.append((tip, toward))
            if b1 and b2:
                d = (f"M{start[0]},{start[1]} C{b1[0]},{b1[1]} "
                     f"{b2[0]},{b2[1]} {end[0]},{end[1]}")
            else:
                d = f"M{start[0]},{start[1]} L{end[0]},{end[1]}"
            parts.append(f'<g class="segment-group" id="s{segment_key}">'
                         f'<path class="segment" d="{d}"/>')
            for (tip_x, tip_y), (from_x, from_y) in heads:
                # the triangle points along +y; turn it to face the node
                rotation = math.degrees(math.atan2(from_y - tip_y, from_x - tip_x)) + 90
                parts.append(
                    f'<g class="arrowheads"><path class="arrowhead" '
                    f'd="M{-ARROW_W / 2},0 L0,{ARROW_H} L{ARROW_W / 2},0 Z" '
                    f'transform="translate({tip_x},{tip_y})rotate({rotation})"/></g>')
            parts.append("</g>")
        parts.append("</g>")
    parts.append('</g><g id="nodes">')
    built = any("member_box" in reaction for reaction in body["reactions"].values())
    for key, node in nodes.items():
        kind = node.get("node_type")
        if kind == "metabolite":
            circle = "metabolite-circle"
            radius = layout.node_radius(node.get("node_is_primary", True))
        else:
            circle = f"{kind}-circle"
            radius = MARKER_R
        parts.append(f'<g class="node" id="n{key}">'
                     f'<circle class="node-circle {circle}" '
                     f'transform="translate({node["x"]},{node["y"]})" r="{radius}"/>')
        if kind == "metabolite" and "label_x" in node:
            label_x, anchor = node["label_x"], ""
            if built and label_x < node["x"]:
                # Escher starts every label at label_x, which for a label the
                # builder put on the left is its estimated width back from the
                # node; ending it at the node's own clearance instead keeps a
                # long or bold name from running into the circle whatever its
                # true width (a hand-placed label is left where it was put)
                label_x = node["x"] - layout.label_offset(node.get("node_is_primary", True))
                anchor = ' text-anchor="end"'
            parts.append(f'<text class="node-label label"{anchor} transform="translate('
                         f'{label_x},{node["label_y"]})">'
                         f'{escape(str(node.get("bigg_id", "")))}</text>')
        parts.append("</g>")
    parts.append('</g><g id="text-labels">')
    for key, label in body.get("text_labels", {}).items():
        parts.append(f'<g class="text-label" id="l{key}">'
                     f'<text class="label" transform="translate('
                     f'{label["x"]},{label["y"]})">{escape(str(label["text"]))}</text></g>')
    parts.append("</g></svg>")
    out_path = Path(out_path)
    out_path.write_text("\n".join(parts), encoding="utf-8")
    return out_path


def render_map_svg(escher_map, out_path=None, dashed=True, layout=None,
                   html=True, **processing):
    """Render a map to a finished SVG figure: member box nodes, member-coloured
    edges, cross-feeding dashed — and, beside it, the interactive HTML figure.

    ``escher_map`` is a map from :mod:`escher_edit.build_map` (or any Escher
    map) as ``[header, body]`` or a path; ``out_path`` defaults to the map's
    path with an ``.svg`` suffix. ``dashed`` dashes the cross-feeding edges
    (:func:`escher_edit.build_map.cross_feeding_segments`). Everything else
    goes to :func:`escher_edit.svg_editor.EscherSVG_processing` — ``style``,
    ``member_colors``, ``member_groups`` (colouring members by group, by
    phylum say, with a legend), ``largeEdgeLabels`` and the rest — whose defaults
    already draw a coloured box for every member, sized and placed as the
    layout reserved it; pass ``member_boxes=False`` to leave them out.
    ``layout`` is the :class:`~escher_edit.build_map.MapStyle` the map was
    built with (a default one if not given): nodes are drawn at its radii and,
    without a ``style``, labels at the sizes it reserved room for. A map laid
    out with a style of its own should pass that style here.

    ``html`` also writes the figure as an interactive page
    (:func:`escher_edit.interactive.write_interactive_html`): hovering an
    edge, a compound or a member box highlights what it connects, and the
    page zooms and pans. ``True`` puts it beside the SVG with an ``.html``
    suffix — replacing any file already there — a path puts it there, and
    ``False`` writes no page.

    Returns the path of the SVG written.
    """
    loaded = _load_map(escher_map)
    if out_path is None:
        if isinstance(escher_map, (list, tuple)):
            raise ValueError("out_path is required when the map is passed "
                             "already loaded")
        out_path = Path(escher_map).with_suffix(".svg")
    out_path = Path(out_path)
    layout = layout or MapStyle()
    if processing.get("style") is None:
        processing["style"] = EscherStyle(
            node_label_px=layout.label_font_px,
            rxn_label_px=layout.reaction_label_font_px)
    if dashed:
        processing.setdefault("dashedEdges", cross_feeding_segments(loaded))
    processing.setdefault("json_path", loaded)
    with tempfile.TemporaryDirectory() as scratch:
        raw = write_escher_svg(loaded, Path(scratch) / "map.svg", layout)
        edited = EscherSVG_processing(raw, **processing)
        shutil.move(str(edited), out_path)
    log.info("rendered %s", out_path)
    if html:
        write_interactive_html(
            out_path, loaded,
            out_path.with_suffix(".html") if html is True else html)
    return out_path


def load_member_groups(path, rank=None):
    """``{member: group}`` from the JSON object at ``path``.

    The object maps each member either to its group or, with ``rank``, to
    its lineage — ``{rank: name}`` or a ``"d__...;p__..."`` string — which
    :func:`escher_edit.palette.taxon_groups` reads the group off.
    """
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a JSON object of member -> group")
    lineages = any(isinstance(value, dict)
                   or (isinstance(value, str) and "__" in value)
                   for value in data.values())
    if rank is not None:
        if not lineages:
            raise ValueError(f"{path} maps members to groups, not lineages, "
                             f"so there is no rank {rank!r} to group by")
        return taxon_groups(data, rank)
    if lineages:
        raise ValueError(f"{path} maps members to lineages; give the rank "
                         f"to group them by")
    return data


def add_member_color_arguments(parser):
    """Add the options for colouring members by group to ``parser``; read
    them back with :func:`member_color_options`."""
    options = parser.add_argument_group("colouring members by group")
    options.add_argument(
        "--member-groups", metavar="JSON", default=None,
        help="Colour members by group rather than one colour each, with a "
             "legend naming the groups: a JSON object of member -> group, or "
             "of member -> lineage with --group-rank ({rank: name} or "
             "'d__...;p__...' strings); members without a group are grey")
    options.add_argument(
        "--group-rank", default=None,
        help="The rank to group members by when --member-groups holds "
             "lineages (e.g. phylum, family, genus)")
    options.add_argument(
        "--group-colors", metavar="JSON", default=None,
        help="A JSON object of group -> #rrggbb fixing those groups' "
             "colours, so a series of maps colours each group alike")
    options.add_argument(
        "--legend-title", default=None,
        help="Title for the group legend (default: the rank, if grouping by "
             "one)")
    options.add_argument(
        "--no-legend", action="store_true",
        help="Leave the group legend out")


def member_color_options(args, parser):
    """The :func:`render_map_svg` keywords for the options
    :func:`add_member_color_arguments` added."""
    if args.member_groups is None:
        for flag, given in (("--group-rank", args.group_rank),
                            ("--group-colors", args.group_colors),
                            ("--legend-title", args.legend_title)):
            if given is not None:
                parser.error(f"{flag} needs --member-groups")
        return {}
    try:
        groups = load_member_groups(args.member_groups, args.group_rank)
    except ValueError as error:
        parser.error(str(error))
    group_colors = (json.loads(Path(args.group_colors).read_text())
                    if args.group_colors else None)
    title = args.legend_title
    if title is None and args.group_rank:
        title = args.group_rank.capitalize()
    return {"member_groups": groups, "group_colors": group_colors,
            "member_legend": False if args.no_legend else (title or True)}


def main():
    parser = argparse.ArgumentParser(
        description="Render editEscher community maps to SVG figures: a "
                    "coloured box per member, its edges in its colour, "
                    "cross-feeding dashed.")
    parser.add_argument("maps", nargs="+", help="Escher map JSON file(s)")
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output SVG path (one map only; default: <map>.svg)")
    parser.add_argument(
        "--no-member-boxes", action="store_true",
        help="Draw members as Escher does, as the vertex their edges meet at")
    parser.add_argument(
        "--no-member-colors", action="store_true",
        help="Draw every edge in one colour")
    parser.add_argument(
        "--no-dashes", action="store_true",
        help="Draw cross-feeding edges solid")
    parser.add_argument(
        "--no-html", action="store_true",
        help="Write only the SVG; by default an interactive HTML figure goes "
             "beside it, where hovering an edge, compound or member "
             "highlights what it connects")
    add_member_color_arguments(parser)
    args = parser.parse_args()
    if args.output and len(args.maps) > 1:
        parser.error("--output takes a single map")

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    coloring = member_color_options(args, parser)
    for path in args.maps:
        render_map_svg(path, args.output, dashed=not args.no_dashes,
                       member_boxes=not args.no_member_boxes,
                       member_colors=not args.no_member_colors,
                       html=not args.no_html, **coloring)


if __name__ == "__main__":
    main()
