"""Edit a rendered Escher SVG: restyle labels/edges, fade non-highlighted
elements, draw ASV node rectangles, and recolor selected segments.

Extracted from ABX_mouse_gut/Notebooks/escher_API_mapping.ipynb, section
"editing the SVG Escher Map". The ``shapely``-based label-layout helpers are
preserved from the notebook but were never enabled in the final call.
"""
import logging
import random
import re
from dataclasses import dataclass
from pathlib import Path

# shapely is only required if the label-layout helpers below are used.
# from shapely.geometry import Point, box


log = logging.getLogger(__name__)


@dataclass
class EscherStyle:
    """Style parameters for :func:`EscherSVG_processing`.

    Hex colors are specified without a leading ``#``. Pixel fields accept
    ``int`` or ``float``; fractional values are rendered as-is into CSS.
    """
    node_label_px: float = 10
    rxn_label_px: float = 30
    node_label_hex: str = "e0865b"
    rxn_label_hex: str = "202078"
    rxn_edge_hex: str = "5aa6df"
    rxn_edge_px: float = 10
    tint_factor: float = 0.5
    stroke_reduction: float = 100


def tint_color(hex_color, tint_factor=0.5):
    """Return a lighter tint of ``hex_color`` (no leading '#') by mixing
    toward white by ``tint_factor`` (0 = unchanged, 1 = white)."""
    rgb_color = tuple(int(hex_color.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4))
    tinted_rgb = [((1 - tint_factor) * c + tint_factor * 255) for c in rgb_color]
    return "{:02x}{:02x}{:02x}".format(
        int(tinted_rgb[0]), int(tinted_rgb[1]), int(tinted_rgb[2]))


def get_text_bbox(text_element, font_size=10):
    """Approximate the bounding box of an SVG text element. Requires shapely."""
    from shapely.geometry import box
    x = float(text_element["x"])
    y = float(text_element["y"])
    width = len(text_element.get_text()) * font_size * 0.6
    height = font_size
    return box(x, y - height, x + width, y)


def get_node_bbox(node_element):   # TODO need to be updated
    from shapely.geometry import box
    cx = float(node_element["cx"])
    cy = float(node_element["cy"])
    r = float(node_element["r"])
    return box(cx - r, cy - r, cx + r, cy + r)


def get_edge_bbox(edge_element):   # TODO need to be updated
    from shapely.geometry import box
    x1 = float(edge_element["x1"])
    y1 = float(edge_element["y1"])
    x2 = float(edge_element["x2"])
    y2 = float(edge_element["y2"])
    return box(min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))


def adjust_labels(labels, nodes, edges, font_size=10, textToKeep=("l27", "l13")):
    """Randomly jitter text labels until they no longer overlap other
    labels/nodes/edges. Disabled in the final notebook call, but preserved
    here for reference. Requires shapely."""
    label_positions = [get_text_bbox(label, font_size)
                       for label in labels if label["id"] not in textToKeep]
    node_positions = [get_node_bbox(node) for node in nodes]
    edge_positions = [get_edge_bbox(edge) for edge in edges]

    all_positions = label_positions + node_positions + edge_positions
    for i, label in enumerate(labels):
        bbox = label_positions[i]
        adjusted = False
        while not adjusted:
            overlap = False
            for j, other_bbox in enumerate(all_positions):
                if i != j and bbox.intersects(other_bbox):
                    overlap = True
                    break
            if overlap:
                match = re.search(r"translate\((-?\d+\.?\d*),\s*(-?\d+\.?\d*)\)",
                                  label["transform"])
                x, y = float(match.group(1)), float(match.group(2))
                x = str(x + random.uniform(-1, 1))
                y = str(y + random.uniform(-1, 1))
                label["transform"] = f"translate({x},{y})"
                bbox = get_text_bbox(label, font_size)
                label_positions[i] = bbox
            else:
                adjusted = True


def _metadata(soup):
    return soup.find("style", attrs={"type": "text/css"})


def remove_label_groups(soup, labels_to_remove):
    """Drop any element whose ``class`` is in ``labels_to_remove`` (e.g.
    stoichiometry labels)."""
    for label in labels_to_remove:
        for lbl in soup.find_all(class_=label):
            lbl.decompose()


def restyle_reaction_labels(soup, style, largeEdgeLabels, abbrevIDs):
    """Rewrite the ``.reaction-label`` CSS rule, add a faded variant when
    ``largeEdgeLabels`` is set, and re-tag non-highlighted reactions with the
    faded class. If ``abbrevIDs`` is set, shorten each reaction label string."""
    metadata = _metadata(soup)
    metadata.string = metadata.string.replace(
        ".reaction-label{font-size:30px;fill:#202078;text-rendering:optimizelegibility}",
        ".reaction-label{" + f"font-size:{style.rxn_label_px}px;fill:#{style.rxn_label_hex};stroke:#000000;stroke-width:{style.rxn_label_px / style.stroke_reduction}px;" + "text-rendering:optimizelegibility}")
    faded_label = None
    if largeEdgeLabels is not None:
        faded_label = tint_color(style.rxn_label_hex, style.tint_factor)
        metadata.string = metadata.string + f" .m{faded_label}-label" + "{" + f"font-size:{style.rxn_label_px / 1.5}px;fill:#{faded_label};stroke:#000000;stroke-width:{style.rxn_label_px / style.stroke_reduction}px;" + "text-rendering:optimizelegibility}"
    for reaction in soup.find_all(class_="reaction"):
        label = reaction.find("text")
        if largeEdgeLabels is not None and reaction["id"] not in largeEdgeLabels:
            label["class"] = f"m{faded_label}-label"
        if abbrevIDs is not None:
            cleanedStr = label.string.split("_")[0].replace("-ABX", "").replace("RC", "").replace("WD", "")
            label.string = cleanedStr.split(".")[0][:6] + "." + cleanedStr.split(".")[1]


def restyle_nodes(soup, style, largeNodeLabels, abbrevIDs, perNodeLabelSizes=None):
    """Rewrite ``.node-label`` and ``.metabolite-circle`` CSS. When
    ``largeNodeLabels`` is set, add a faded class and apply it to every node
    outside that list; when ``abbrevIDs`` is also set, rename faded node
    labels via that mapping.

    If ``perNodeLabelSizes`` is provided (``{node_id: size_px}``), emit
    inline ``font-size`` styles on the matching node ``<text>`` elements.
    The node_id matches the SVG ``id`` attribute of each ``.node`` group
    (Escher prefixes JSON keys with ``n`` — e.g. JSON key ``"33"`` becomes
    SVG id ``"n33"`` — so both forms are checked).
    """
    metadata = _metadata(soup)
    metadata.string = metadata.string.replace(
        ".node-label{font-size:20px}",
        ".node-label{" + f"font-size:{style.node_label_px}px;fill:#{style.node_label_hex};stroke:#000000;stroke-width:{style.node_label_px / style.stroke_reduction}px;" + "}")
    metadata.string = metadata.string.replace(
        ".metabolite-circle{stroke:#a24510;fill:#e0865b}",
        ".metabolite-circle{stroke:#" + f"{style.node_label_hex};fill:#{style.node_label_hex}" + "}")

    faded = None
    if largeNodeLabels is not None:
        faded = tint_color(style.node_label_hex, style.tint_factor)
        for node in soup.find_all(class_="node"):
            circle = node.find("circle")
            label = node.find("text")
            if label is None or circle is None:
                continue
            if node["id"] in largeNodeLabels:
                continue
            circle["class"] = f"m{faded}"
            label["class"] = f"m{faded}-label"
            if abbrevIDs is not None:
                label.string = abbrevIDs[label.string]
                log.debug("faded and re-abbreviated node %s", node["id"])
        metadata.string = metadata.string + f" .m{faded}" + "{" + f"stroke:#{faded};fill:#{faded}" + "}"
        metadata.string = metadata.string + f" .m{faded}-label" + "{" + f"font-size:{style.node_label_px / 1.2}px;fill:#{faded};stroke:#000000;stroke-width:{style.node_label_px / style.stroke_reduction}px;" + "}"

    if perNodeLabelSizes:
        # Build a matcher accepting both raw keys (e.g. "33") and the
        # ``n``-prefixed form Escher uses in SVG ids ("n33").
        size_by_id = {}
        for key, size in perNodeLabelSizes.items():
            size_by_id[str(key)] = size
            size_by_id[f"n{key}"] = size
        stroke = style.node_label_px / style.stroke_reduction
        for node in soup.find_all(class_="node"):
            size = size_by_id.get(node.get("id"))
            if size is None:
                continue
            label = node.find("text")
            if label is None:
                continue
            existing = label.get("style", "")
            label["style"] = (
                existing.rstrip(";") + ";" if existing else ""
            ) + f"font-size:{size}px;stroke-width:{stroke}px"
            log.debug("enlarged label on node %s to %gpx", node.get("id"), size)


def restyle_segments(soup, style):
    """Rewrite the ``.segment`` CSS rule and the miscellaneous
    ``.text-label-input`` font-size."""
    metadata = _metadata(soup)
    metadata.string = metadata.string.replace(
        ".segment{" + "stroke:#334E75;stroke-width:10px;fill:none}",
        ".segment{" + f"stroke:#{style.rxn_edge_hex};" + f"stroke-width:{style.rxn_edge_px}px;" + "fill:none}")
    metadata.string = metadata.string.replace(
        ".text-label-input{font-size:50px}",
        ".text-label-input{font-size:" + f"{style.rxn_label_px * 1.8}" + "px}")


def dash_segments(soup, dashed_edges, dash_pattern="8,5"):
    """Dash (or dot) the listed segment groups.

    ``dashed_edges`` holds segment-group ids — ``s<id>``, as
    ``model_mapping.build_direction_tracking`` returns for consumption edges
    and ``build_map.cross_feeding_segments`` for cross-feeding edges. The
    pattern is written as an SVG ``stroke-dasharray`` attribute on the
    ``<g class="segment-group">``, which the ``<path class="segment">`` inside
    inherits; pass e.g. ``"2,6"`` for a dotted line.

    Returns the number of segments actually matched, and logs it against the
    number requested so a stale id list does not fail silently.
    """
    if not dashed_edges:
        return 0
    targets = set(dashed_edges)
    matched = 0
    for edge in soup.find_all(class_="segment-group"):
        if edge.get("id") in targets:
            edge["stroke-dasharray"] = dash_pattern
            matched += 1
    log.info("dashed %d of %d requested segment(s) with pattern %r",
             matched, len(targets), dash_pattern)
    return matched


def mark_asv_rectangles(soup, style):
    """Draw a small rectangle over the start of every segment whose path is a
    short, nearly-horizontal straight line, marking it as an ASV anchor.

    The offsets are tuned to the ABX map's viewBox, so this is only meaningful
    for that figure; pass ``mark_asv_nodes=False`` to
    :func:`EscherSVG_processing` for any other map."""
    for edge in soup.find_all(class_="segment-group"):
        segment = edge.find(class_="segment")
        coordinates = list(map(float, re.split(",| ", segment["d"]
                                               .replace("M", "").replace("M", "")
                                               .replace("C", "")
                                               .replace("undefined", "100000000"))))
        if abs(coordinates[0] - coordinates[2]) < 30 and abs(coordinates[1] - coordinates[3]) < 3:
            rect = soup.new_tag("rect")
            rect["x"] = str(coordinates[0] + float("10734.555053710938") / 2.03)
            rect["y"] = str(coordinates[1] + float("5315.0328369140625") / 7.66)
            rect["width"] = "40"
            rect["height"] = "15"
            rect["fill"] = f"#{style.rxn_label_hex}"
            soup.svg.append(rect)


def apply_color_highlights(soup, colorElements, largeEdgeLabels, style):
    """Apply per-color recoloring of segments, reactions, and nodes per
    ``colorElements`` (``{hex_color: [element_ids]}``). When
    ``largeEdgeLabels`` is set, non-highlighted segments are additionally
    tagged with the faded class."""
    if colorElements is None:
        return
    metadata = _metadata(soup)
    faded_edges = (tint_color(style.rxn_edge_hex, style.tint_factor)
                   if largeEdgeLabels is not None else None)

    edge_color = {}
    for color, lst in colorElements.items():
        edge_color.update({x: color for x in lst})
    for edge in soup.find_all(class_="segment-group"):
        if edge["id"] in edge_color:
            color = edge_color[edge["id"]]
            for child in edge.find_all(class_="segment"):
                child["class"] = f"m{color}"
        elif largeEdgeLabels is not None:
            for child in edge.find_all(class_="segment"):
                child["class"] = f"m{faded_edges}"
    for color, elements in colorElements.items():
        metadata.string = metadata.string + f" .m{color}" + "{" + f"stroke:#{color};" + f"stroke-width:{style.rxn_edge_px * 1.2}px;" + "fill:none}"
        if any("r" in x for x in elements):
            for ele in elements:
                if "r" in ele:
                    reaction = soup.find(attrs={"id": ele})
                    label = reaction.find("text")
                    label["class"] = f"m{color}-label"
            metadata.string = metadata.string + f" .m{color}-label" + "{" + f"font-size:{style.rxn_label_px}px;fill:#{color};stroke:#000000;stroke-width:{style.rxn_label_px / style.stroke_reduction}px;" + "text-rendering:optimizelegibility;font-family:sans-serif;font-style:italic;font-weight:700;}"
        if any("n" in x for x in elements):
            for ele in elements:
                if "n" in ele:
                    node = soup.find(attrs={"id": ele})
                    circle = node.find("circle")
                    circle["class"] = f"m{color}-circle"
                    label = node.find("text")
                    label["class"] = f"m{color}-label"
            metadata.string = metadata.string + f" .m{color}-circle" + "{" + f"stroke:#{color};fill:#{color}" + "}"

    if largeEdgeLabels is not None:
        metadata.string = metadata.string + f" .m{faded_edges}" + "{" + f"stroke:#{faded_edges};" + f"stroke-width:{style.rxn_edge_px / 1.2}px;" + "fill:none}"


def extract_per_node_label_sizes_from_json(json_path):
    """Read an Escher JSON map and return ``{node_id: label_size}`` for every
    metabolite node that carries the (non-standard) ``label_size`` attribute
    (written by ``layout.set_node_label_sizes``).
    """
    from json import load
    with open(json_path) as fh:
        escher_map = load(fh)
    return {
        node_id: node["label_size"]
        for node_id, node in escher_map[1]["nodes"].items()
        if node.get("node_type") == "metabolite" and "label_size" in node
    }


def EscherSVG_processing(svg_path="metabolite_focused_map.svg",
                         style=None,
                         labels_to_remove=("stoichiometry-labels",),
                         largeNodeLabels=None,
                         largeEdgeLabels=None,
                         colorElements=None,
                         abbrevIDs=None,
                         dashedEdges=None,
                         dash_pattern="8,5",
                         mark_asv_nodes=True,
                         perNodeLabelSizes=None,
                         json_path=None):
    """Post-process a rendered Escher SVG.

    Drops specified label groups, rewrites fonts/colors for nodes, reactions,
    and segments, fades everything outside ``largeNodeLabels`` /
    ``largeEdgeLabels``, draws small rectangles to mark short straight
    segments as ASV anchor nodes, applies ``colorElements`` (a dict of
    hex color -> [element ids]) to highlight specific segments/nodes/
    reactions, and dashes the segments in ``dashedEdges`` with
    ``dash_pattern`` (``"2,6"`` or similar gives dots).

    ``mark_asv_nodes`` gates the ASV anchor rectangles, whose offsets only
    make sense for the ABX map; turn it off for any other figure.

    When ``perNodeLabelSizes`` is given (or ``json_path`` points to a map
    written by ``layout.set_node_label_sizes``), each listed node has its
    ``<text>`` label enlarged to the requested pixel size via an inline
    ``style`` attribute. ``perNodeLabelSizes`` overrides ``json_path``.

    Style parameters (sizes, colors, tint/stroke scaling) are bundled in
    :class:`EscherStyle`; pass an instance via ``style=`` to override
    defaults. The edited SVG is written to ``<stem>_edited<suffix>``.
    """
    from bs4 import BeautifulSoup
    if style is None:
        style = EscherStyle()
    svg_path = Path(svg_path)
    soup = BeautifulSoup(svg_path.read_text(), "lxml-xml")

    if perNodeLabelSizes is None and json_path is not None:
        perNodeLabelSizes = extract_per_node_label_sizes_from_json(json_path)

    remove_label_groups(soup, labels_to_remove)
    restyle_reaction_labels(soup, style, largeEdgeLabels, abbrevIDs)
    restyle_nodes(soup, style, largeNodeLabels, abbrevIDs,
                  perNodeLabelSizes=perNodeLabelSizes)
    restyle_segments(soup, style)
    if mark_asv_nodes:
        mark_asv_rectangles(soup, style)
    dash_segments(soup, dashedEdges, dash_pattern)
    apply_color_highlights(soup, colorElements, largeEdgeLabels, style)

    out_path = svg_path.with_name(svg_path.stem + "_edited" + svg_path.suffix)
    out_path.write_text(soup.prettify(), encoding="utf-8")
    log.info("wrote edited Escher SVG to %s", out_path)
    return out_path


def main():
    # The original invocation from the notebook. ``consumptionEdges`` comes
    # from escher_edit.model_mapping.build_direction_tracking.
    from .model_mapping import build_direction_tracking

    logging.basicConfig(level=logging.INFO)
    _, consumptionEdges = build_direction_tracking()

    style = EscherStyle(
        node_label_px=90, rxn_label_px=90, rxn_edge_px=6,
        node_label_hex="E36400", rxn_label_hex="cc0000", rxn_edge_hex="2C79BD",
        tint_factor=0.8, stroke_reduction=100,
    )

    EscherSVG_processing(
        "metabolite_focused_map_IDs_cleaned.svg",
        style=style,
        largeNodeLabels=[
            "n202", "n311", "n900", "n963", "n33", "n32", "n830", "n71",
            "n1322", "n1443", "n1672", "n1673",
            "n2066", "n2102", "n206", "n404", "n873", "n1090", "n1182",
            "n1746", "n2154", "n2154", "n1338",  # RC
            "n539", "n1063", "n1282", "n1511", "n1576",  # cellibiose
            "n1540", "n249", "n1515", "n1540",  # Melitose
            "n85", "n104",  # WD
        ],
        largeEdgeLabels=[
            "r10", "r23", "r48", "r50", "r49", "r41", "r3", "r60", "r64",
            "r65", "r73", "r68", "r72", "r88", "r86", "r83",
            "r95", "r108", "r106", "r107",
            "r31", "r53", "r63", "r75", "r76", "r99", "r77", "r79", "r121",
            "s1435", "s1009", "s1218",  # WD
        ],
        colorElements={
            "8E0FFF": [
                "s903", "s295", "s444", "s906", "s928", "s900", "s913", "s67",
                "s1246", "s1150", "s1254", "s1238", "s1390",
                "s1370", "s1666", "s1637", "s1588", "s788", "s66", "s1589",
                "s1787", "s1777", "s387", "s915",
                "s1170", "s1186", "s1152", "s2024", "s2072", "s2042", "s2024",
                "s2085", "s1993", "s1190", "s1669",
                "s1959", "s1975",  # RC
                "s602", "s593", "s624", "s597", "s511", "s2195", "s1479", "s98"  # WD
            ],  # purple
            "11C911": [
                "s191", "s195", "r10", "r3",       # "r68", Enterococcus.3 in the 4th interval
                "r43", "n206", "s387",             # "r20", Allio in the first interval
                "s383", "s823", "s829", "n873", "r80", "r122",
            ],  # green
            "F15757": ["s1304"],  # red
        },
        dashedEdges=consumptionEdges,
    )


if __name__ == "__main__":
    main()
