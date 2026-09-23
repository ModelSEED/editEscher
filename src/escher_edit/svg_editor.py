"""Edit a rendered Escher SVG: restyle labels/edges, fade non-highlighted
elements, draw member boxes, and recolor selected segments.

Extracted from ABX_mouse_gut/Notebooks/escher_API_mapping.ipynb, section
"editing the SVG Escher Map". The ``shapely``-based label-layout helpers are
preserved from the notebook but were never enabled in the final call.
"""
import logging
import random
import re
from dataclasses import dataclass
from pathlib import Path

from .build_map import member_box_rects
from .palette import (INK, group_legend, member_colors as _member_colors,
                      text_color_on)

# shapely is only required if the label-layout helpers below are used.
# from shapely.geometry import Point, box


log = logging.getLogger(__name__)


@dataclass
class EscherStyle:
    """Style parameters for :func:`EscherSVG_processing`.

    Hex colors are specified without a leading ``#``. Pixel fields accept
    ``int`` or ``float``; fractional values are rendered as-is into CSS.

    Metabolite nodes default to ink: member colours (see
    :mod:`escher_edit.palette`) are what distinguish the members, and Escher's
    own orange nodes would read as belonging to the second member. The label
    and edge sizes are what ``build_map.MapStyle`` reserves room for, so change
    them together.
    """
    node_label_px: float = 9.8
    rxn_label_px: float = 29.4
    node_label_hex: str = "222222"
    rxn_label_hex: str = "202078"
    rxn_edge_hex: str = "5aa6df"
    rxn_edge_px: float = 14
    tint_factor: float = 0.5
    stroke_reduction: float = 100
    #: Member-box fill, border, corner rounding and the room left around the
    #: label inside it, for a reaction without a member colour (a member
    #: colour fills its box and the border is left white). An empty
    #: ``member_box_edge_hex`` takes the border from the reaction-label
    #: colour, so box, border and text read as one node.
    member_box_fill_hex: str = "ffffff"
    member_box_edge_hex: str = ""
    member_box_edge_px: float = 3
    member_box_pad: float = 12
    member_box_radius: float = 8


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


def _load_map(json_path):
    """The map at ``json_path``, or ``json_path`` itself if it is already a
    loaded ``[header, body]`` map."""
    if isinstance(json_path, (list, tuple)):
        return json_path
    from json import load
    with open(json_path) as fh:
        return load(fh)


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


def dash_segments(soup, dashed_edges, dash_pattern="11.2,7"):
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


def _path_points(d):
    """The coordinate pairs in an SVG path's ``d``, command letters dropped."""
    numbers = re.findall(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", d)
    return list(zip((float(n) for n in numbers[::2]),
                    (float(n) for n in numbers[1::2])))


def _reaction_anchor(reaction, max_length=30.0, max_rise=3.0):
    """Where a reaction sits: the point its marker segments meet.

    Escher draws ``multimarker -> midmarker`` as a short, flat, two-point
    line, one on each side of the reaction, and both of them touch the
    midmarker. Which end of a segment that is cannot be read off the path, so
    the midmarker is identified as the point those short segments have in
    common. Only one of them qualifies when the other is drawn long — a
    reaction can be left connected straight to a metabolite by filtering, and
    a metabolite drawn close in makes that segment short and flat too — so
    two fallbacks follow: the end no other segment of the reaction touches
    (the multimarker at the other end carries the metabolite edges, the
    midmarker carries none), and failing that the segment's own endpoint,
    since Escher and :mod:`escher_edit.build_map` both write these segments
    as ``multimarker -> midmarker``.

    This is what the notebook's version keyed on as well; it just added the
    answer to offsets tuned to one figure's viewBox instead of using it.

    Returns ``(x, y)``, or None when no rule leaves a single candidate.
    """
    marker_ends, other_ends = [], set()
    # by group, not by the "segment" class: apply_color_highlights replaces
    # that class on every segment it recolours or fades
    for group in reaction.find_all(class_="segment-group"):
        segment = group.find("path", recursive=False)
        if segment is None:
            continue
        d = segment.get("d", "")
        points = [(round(x, 3), round(y, 3)) for x, y in _path_points(d)]
        flat = (len(points) == 2 and "C" not in d and "undefined" not in d
                and abs(points[1][0] - points[0][0]) <= max_length
                and abs(points[1][1] - points[0][1]) <= max_rise)
        if flat:
            marker_ends += points
        else:
            other_ends.update(points)

    shared = {point for point in marker_ends if marker_ends.count(point) > 1}
    if len(shared) == 1:
        return shared.pop()
    untouched = set(marker_ends) - other_ends
    if len(untouched) == 1:
        return untouched.pop()
    return marker_ends[1] if len(marker_ends) == 2 else None


def draw_member_boxes(soup, style, largeEdgeLabels=None, box_colors=None,
                      rects=None, names=None):
    """Draw a labelled box on every reaction anchor.

    Each community member is one reaction, and without this it is drawn as
    nothing but the vertex where its edges meet, with its name floating above
    them. The box makes it a node: a filled rectangle on the anchor (see
    :func:`_reaction_anchor`), sized to its own label, with that label moved
    inside it. Escher's JSON schema has no box node — marker circles and a
    reaction label are all it can express — so, like the edge dashing, this
    happens on the rendered SVG.

    Both carry ``data-reaction`` (the reaction's id), which is how the
    interactive figure (:mod:`escher_edit.interactive`) finds a member's box.

    The box and the label go into a group appended last, so they sit above
    the edges rather than being crossed by them, and the box hides the marker
    circles it covers. That takes each label out of its ``<g class="reaction">``,
    so this must run after everything that reaches a reaction's label through
    that group — see :func:`EscherSVG_processing`.

    ``largeEdgeLabels`` is the same list :func:`restyle_reaction_labels`
    takes: the reactions outside it are drawn at a reduced font size, so
    their boxes are sized to match rather than to the full one.

    ``rects`` (``{reaction id: (x, y, width, height)}``, from
    :func:`escher_edit.build_map.member_box_rects`) replaces both guesses for
    the reactions it lists: the box is the one the layout fitted the edges
    around, exactly where the layout put it, rather than one sized to the
    label on an anchor read off the segments. ``names`` (``{reaction id:
    member}``) then labels it with the member, which is what the box stands
    for — Escher's own label is the ``bigg_id``, which in a stacked map
    repeats the condition the block caption already gives. A label too wide
    for its box widens it, with a warning, since the edges were fitted to the
    narrower one.

    Every edge of a member stops at the border of its box (see
    :func:`_trim_to_box`) rather than running on underneath it to the markers
    at its centre, and the short marker segments wholly inside it are hidden.

    ``box_colors`` (``{reaction id: "#rrggbb"}``) fills each listed member's
    box with its colour — tinted like its edges when it is outside
    ``largeEdgeLabels`` — rings it in white so the edges arriving at it stop
    short of the fill, and writes the label in white or ink, whichever reads
    better on that fill. Reactions not listed keep the ``member_box_*`` style.

    Returns the number of boxes drawn.
    """
    boxes = soup.new_tag("g")
    boxes["id"] = "member-boxes"
    edge_hex = style.member_box_edge_hex or style.rxn_label_hex

    rects, names = rects or {}, names or {}
    widened = trimmed = 0
    for reaction in soup.find_all(class_="reaction"):
        rid = reaction.get("id")
        label_group = reaction.find(class_="reaction-label-group")
        label = label_group.find("text") if label_group else None
        if label is not None and rid in rects and names.get(rid):
            label.string = names[rid]
        name = label.get_text().strip() if label is not None else ""

        font_px = style.rxn_label_px
        faded = largeEdgeLabels is not None and rid not in largeEdgeLabels
        if faded:
            font_px /= 1.5   # the size restyle_reaction_labels fades them to
        text_width = 0.6 * font_px * len(name)
        if rid in rects:
            left, top, width, height = rects[rid]
            x, y = left + width / 2, top + height / 2
            if text_width + 2 * style.member_box_pad > width:
                width = text_width + 2 * style.member_box_pad
                widened += 1
        else:
            anchor = _reaction_anchor(reaction)
            if anchor is None:
                continue
            x, y = anchor
            height = font_px + 2 * style.member_box_pad
            # a nameless reaction would otherwise get a sliver of a box
            width = max(text_width + 2 * style.member_box_pad, height)
        box = soup.new_tag("rect")
        box["class"] = "member-box"
        box["x"] = f"{x - width / 2:.3f}"
        box["y"] = f"{y - height / 2:.3f}"
        box["width"] = f"{width:.3f}"
        box["height"] = f"{height:.3f}"
        box["rx"] = str(style.member_box_radius)
        color = (box_colors or {}).get(reaction.get("id"))
        if color is None:
            box["fill"] = f"#{style.member_box_fill_hex}"
            box["stroke"] = f"#{edge_hex}"
        else:
            fill = "#" + tint_color(color, style.tint_factor) if faded else color
            box["fill"] = fill
            box["stroke"] = "#ffffff"
            if label is not None:
                _add_style(label, fill=text_color_on(fill), stroke="none")
        box["stroke-width"] = str(style.member_box_edge_px)
        box["data-reaction"] = rid
        boxes.append(box)
        trimmed += _trim_to_box(reaction, (x, y, width / 2, height / 2))

        if label_group is not None:
            label_group["data-reaction"] = rid
            # Escher anchors label text on the left, on an alphabetic
            # baseline, so centring it in the box is done by hand.
            label_group["transform"] = (
                f"translate({x - text_width / 2:.3f},"
                f"{y + 0.35 * font_px:.3f})")
            boxes.append(label_group)

    drawn = boxes.find_all("rect")
    if drawn:
        soup.svg.append(boxes)
    log.info("cut %d member edge(s) back to the border of their box", trimmed)
    if widened:
        log.warning("%d member label(s) are wider than the box the layout "
                    "reserved, so their boxes were widened past the edges "
                    "fitted to them — keep EscherStyle.rxn_label_px at or "
                    "below MapStyle.reaction_label_font_px", widened)
    overlapping = _overlapping_boxes(drawn)
    if overlapping:
        log.warning("%d of %d member boxes overlap a neighbour — the "
                    "reactions are closer together than a %gpx label needs; "
                    "lower the label size or box padding, or space the "
                    "reactions further apart (MapStyle.member_pitch for a "
                    "built map)", overlapping, len(drawn), style.rxn_label_px)
    log.info("drew %d member box(es)", len(drawn))
    return len(drawn)


def _cubic_point(points, t):
    (x0, y0), (x1, y1), (x2, y2), (x3, y3) = points
    u = 1.0 - t
    return (u ** 3 * x0 + 3 * u * u * t * x1 + 3 * u * t * t * x2 + t ** 3 * x3,
            u ** 3 * y0 + 3 * u * u * t * y1 + 3 * u * t * t * y2 + t ** 3 * y3)


def _cubic_tail(points, t):
    """The part of a cubic Bezier from ``t`` to its end, as four points."""
    def lerp(a, b):
        return a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t
    p0, p1, p2, p3 = points
    q0, q1, q2 = lerp(p0, p1), lerp(p1, p2), lerp(p2, p3)
    r0, r1 = lerp(q0, q1), lerp(q1, q2)
    return [lerp(r0, r1), r1, q2, p3]


def _trim_to_box(reaction, box, steps=200):
    """Cut a member's edges back to where they cross the border of its box.

    Escher draws every edge of a reaction from a marker node, and the markers
    sit inside the member's box, so without this each edge carries on under
    the box to its centre — hidden by the fill, but still there in the SVG,
    and visible wherever the box is not opaque or is moved. ``box`` is
    ``(centre x, centre y, half width, half height)``. An edge with one end
    inside is cut where it first leaves the box, walking out from that end;
    a segment lying wholly inside — the short marker-to-midmarker links — is
    hidden. Paths that are neither a line nor a single cubic are left alone.

    Returns how many edges were cut.
    """
    cx, cy, half_w, half_h = box

    def inside(point):
        return abs(point[0] - cx) <= half_w and abs(point[1] - cy) <= half_h

    cut = 0
    for group in reaction.find_all(class_="segment-group"):
        path = group.find("path", recursive=False)
        d = path.get("d", "") if path is not None else ""
        points = _path_points(d)
        if "C" in d and len(points) == 4:
            curve = points
        elif "C" not in d and len(points) == 2:
            curve = [points[0], points[0], points[1], points[1]]
        else:
            continue
        start_in, end_in = inside(curve[0]), inside(curve[3])
        if not (start_in or end_in):
            continue
        if start_in and end_in:
            if all(inside(_cubic_point(curve, i / steps)) for i in range(steps + 1)):
                _add_style(path, display="none")
            continue
        if end_in:
            curve = curve[::-1]          # always walk out from the inside end
        # first sample outside, then bisect back to the border
        for i in range(1, steps + 1):
            if not inside(_cubic_point(curve, i / steps)):
                lo, hi = (i - 1) / steps, i / steps
                break
        for _ in range(40):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if inside(_cubic_point(curve, mid)) else (lo, mid)
        tail = _cubic_tail(curve, hi)
        if end_in:
            tail = tail[::-1]
        if "C" in d:
            (x0, y0), (x1, y1), (x2, y2), (x3, y3) = tail
            path["d"] = f"M{x0:.3f},{y0:.3f} C{x1:.3f},{y1:.3f} {x2:.3f},{y2:.3f} {x3:.3f},{y3:.3f}"
        else:
            (x0, y0), (x3, y3) = tail[0], tail[3]
            path["d"] = f"M{x0:.3f},{y0:.3f} L{x3:.3f},{y3:.3f}"
        cut += 1
    return cut


def _overlapping_boxes(boxes):
    """How many boxes run into another one — a label too big for the layout."""
    rects = sorted((float(b["y"]), float(b["y"]) + float(b["height"]),
                    float(b["x"]), float(b["x"]) + float(b["width"]))
                   for b in boxes)
    hit = set()
    for index, (top, bottom, left, right) in enumerate(rects):
        for other in rects[index + 1:]:
            if other[0] >= bottom:
                break
            if other[2] < right and left < other[3]:
                hit.add((top, left))
                hit.add((other[0], other[2]))
    return len(hit)


def _add_style(tag, **properties):
    """Set CSS properties in ``tag``'s inline ``style``, which outranks the
    stylesheet — Escher's rules are scoped as ``svg.escher-svg .segment``, so
    a plain class rule added after them would lose."""
    declarations = {}
    for item in tag.get("style", "").split(";"):
        if ":" in item:
            key, value = item.split(":", 1)
            declarations[key.strip()] = value.strip()
    declarations.update({key.replace("_", "-"): str(value)
                         for key, value in properties.items()})
    tag["style"] = ";".join(f"{key}:{value}"
                            for key, value in declarations.items())


def color_member_edges(soup, reaction_colors, style, largeEdgeLabels=None,
                       colorElements=None):
    """Draw every member's edges, arrowheads included, in its colour.

    ``reaction_colors`` is ``{reaction id: "#rrggbb"}`` (see
    :func:`member_reaction_colors`). Each segment belongs to exactly one
    reaction, so a cross-fed compound's two edges come out in the producer's
    and the consumer's colour respectively.

    Fading follows :func:`apply_color_highlights`: when ``largeEdgeLabels``
    is given, a segment outside it is tinted towards white — from its own
    member's colour, so faded edges still say whose they are — and drawn
    thinner. Segments listed in ``colorElements`` are left for
    :func:`apply_color_highlights` to recolour.

    Returns the ids of the segment groups coloured, which
    :func:`apply_color_highlights` then leaves alone.
    """
    highlighted = {element for elements in (colorElements or {}).values()
                   for element in elements}
    colored = set()
    for reaction in soup.find_all(class_="reaction"):
        color = reaction_colors.get(reaction.get("id"))
        if color is None:
            continue
        for group in reaction.find_all(class_="segment-group"):
            group_id = group.get("id")
            if group_id in highlighted:
                continue
            faded = largeEdgeLabels is not None and group_id not in largeEdgeLabels
            hex_color = "#" + tint_color(color, style.tint_factor) if faded else color
            width = style.rxn_edge_px / 1.2 if faded else style.rxn_edge_px
            for segment in group.find_all(class_="segment"):
                _add_style(segment, stroke=hex_color, stroke_width=f"{width}px")
            for head in group.find_all(class_="arrowhead"):
                _add_style(head, fill=hex_color)
            colored.add(group_id)
    log.info("coloured %d segment(s) across %d member(s)", len(colored),
             len(set(reaction_colors.values())))
    return colored


def _reaction_members(soup, escher_map=None):
    """``{reaction id: member}``: named by the map JSON when given — a
    reaction's ``name`` — and otherwise by the label Escher drew."""
    if escher_map is not None:
        return {f"r{key}": rxn.get("name") or rxn.get("bigg_id", "")
                for key, rxn in escher_map[1]["reactions"].items()}
    names = {}
    for reaction in soup.find_all(class_="reaction"):
        label = reaction.find(class_="reaction-label")
        if label is not None and label.get_text().strip():
            names[reaction.get("id")] = label.get_text().strip()
    return names


def member_reaction_colors(soup, colors=True, json_path=None, groups=None,
                           group_colors=None):
    """``{reaction id: "#rrggbb"}`` for the member reactions of a rendered map.

    ``colors`` is either ``True`` — colour the map's own members from
    :mod:`escher_edit.palette`, in sorted order — or a ``{member: colour}``
    mapping (from :func:`escher_edit.palette.member_colors`) to share one
    assignment across several maps. With ``True``, ``groups`` (``{member:
    group}`` or a function of the member) colours members by group instead,
    ``group_colors`` fixing any group's colour — see
    :func:`escher_edit.palette.member_colors`. Members are named by the map
    JSON at ``json_path`` when given: a reaction's ``name`` is the member in
    every condition block, where the label Escher draws is its ``bigg_id``,
    which carries the condition as well. Without the JSON the label is all
    there is, so a stacked map colours a member differently in each block;
    pass ``json_path`` for those. Call this before
    :func:`restyle_reaction_labels`, which may abbreviate the labels.
    """
    names = _reaction_members(
        soup, _load_map(json_path) if json_path is not None else None)
    if colors is True:
        colors = _member_colors(sorted(set(names.values())), groups,
                                group_colors)
    elif group_colors is not None:
        log.warning("member colours were given as a mapping, so group_colors "
                    "has no effect")
    return {rid: colors[name] for rid, name in names.items() if name in colors}


def _translate(transform):
    """The ``(x, y)`` of a ``translate(x,y)`` transform, or None."""
    found = re.search(r"translate\(\s*([-\d.eE+]+)[\s,]+([-\d.eE+]+)",
                      transform or "")
    return (float(found.group(1)), float(found.group(2))) if found else None


def _drawing_bottom(soup, style):
    """The lowest point anything on the map is drawn at — circles, labels,
    member boxes and edges (their control points, so perhaps a little below
    the curve) — or None for an empty map."""
    lowest = []
    for circle in soup.find_all("circle"):
        at = _translate(circle.get("transform"))
        if at is not None:
            lowest.append(at[1] + float(circle.get("r", 0)))
    descent = 0.3 * max(style.node_label_px, style.rxn_label_px)
    labels = soup.find_all("text") + soup.find_all("g", class_="reaction-label-group")
    for label in labels:
        at = _translate(label.get("transform"))
        if at is not None:
            lowest.append(at[1] + descent)
    for box in soup.find_all("rect", class_="member-box"):
        lowest.append(float(box["y"]) + float(box["height"]))
    for segment in soup.find_all("path", class_="segment"):
        points = _path_points(segment.get("d", ""))
        lowest.extend(y for _, y in points)
    return max(lowest) if lowest else None


def draw_member_legend(soup, entries, style, title=None):
    """Draw a legend of ``entries`` (``[(label, "#rrggbb")]``, from
    :func:`escher_edit.palette.group_legend`) below the map.

    When members are coloured by group, the colour of a box says which group
    its member is in, and nothing on the map says which group a colour is:
    this does. A swatch and the group's name per entry, at the member-label
    size, in a grid of as many equal columns as fit across the canvas,
    centred under the lowest thing drawn, with ``title`` (a rank, say) above
    it. The canvas is extended by the room the legend takes, so the margin
    below the map is kept below the legend.

    Returns the number of entries drawn.
    """
    if not entries:
        return 0
    font_px = style.rxn_label_px
    canvas = soup.find(id="canvas")
    root = soup.svg
    if canvas is not None:
        left, width = float(canvas["x"]), float(canvas["width"])
        top, height = float(canvas["y"]), float(canvas["height"])
    elif root.get("viewBox"):
        left, top, width, height = map(float, root["viewBox"].split())
    else:
        log.warning("no canvas or viewBox to place the member legend on")
        return 0
    bottom = top + height
    drawn_to = _drawing_bottom(soup, style)
    drawn_to = bottom if drawn_to is None else min(drawn_to, bottom)

    swatch, text_gap, entry_gap = font_px, 0.4 * font_px, 1.6 * font_px
    row_height = 1.8 * font_px
    # a grid of equal columns, as many as fit across the canvas, filled in
    # reading order, so the entries line up whatever their names' lengths
    column = (swatch + text_gap + entry_gap
              + max(0.58 * font_px * len(label) for label, _ in entries))
    columns = max(1, min(len(entries),
                         int((width - 4 * font_px + entry_gap) // column)))
    rows = -(-len(entries) // columns)
    grid_left = left + (width - (columns * column - entry_gap)) / 2

    legend = soup.new_tag("g")
    legend["id"] = "member-legend"
    first_row = drawn_to + 2 * font_px
    lines = [(grid_left, first_row, str(title), None)] if title else []
    for index, (label, color) in enumerate(entries):
        row, col = divmod(index, columns)
        lines.append((grid_left + col * column,
                      first_row + (row + bool(title)) * row_height,
                      label, color))
    for x, row_top, label, color in lines:
        middle = row_top + row_height / 2
        if color is not None:
            rect = soup.new_tag("rect")
            rect["class"] = "member-legend-swatch"
            rect["x"] = f"{x:.3f}"
            rect["y"] = f"{middle - swatch / 2:.3f}"
            rect["width"] = rect["height"] = f"{swatch:.3f}"
            rect["rx"] = f"{min(style.member_box_radius, swatch / 4):g}"
            rect["fill"] = color
            legend.append(rect)
            x += swatch + text_gap
        text = soup.new_tag("text")
        text["class"] = "member-legend-label"
        text["x"] = f"{x:.3f}"
        text["y"] = f"{middle + 0.35 * font_px:.3f}"
        text.string = label
        _add_style(text, font_family="sans-serif", font_size=f"{font_px}px",
                   fill=INK, stroke="none",
                   font_weight="700" if color is None else "400")
        legend.append(text)
    soup.svg.append(legend)

    grow = first_row + (rows + bool(title)) * row_height - drawn_to
    if canvas is not None:
        canvas["height"] = f"{height + grow:g}"
    if root.get("viewBox"):
        vx, vy, vw, vh = map(float, root["viewBox"].split())
        root["viewBox"] = f"{vx:g} {vy:g} {vw:g} {vh + grow:g}"
    try:
        root["height"] = f"{float(root['height']) + grow:g}"
    except (KeyError, ValueError):
        pass
    log.info("drew a member legend of %d group(s) in %d row(s)",
             len(entries), rows)
    return len(entries)


def apply_color_highlights(soup, colorElements, largeEdgeLabels, style,
                           member_colored=()):
    """Apply per-color recoloring of segments, reactions, and nodes per
    ``colorElements`` (``{hex_color: [element_ids]}``). When
    ``largeEdgeLabels`` is set, non-highlighted segments are additionally
    tagged with the faded class, except those in ``member_colored``, which
    :func:`color_member_edges` has already faded in their member's colour."""
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
        elif largeEdgeLabels is not None and edge["id"] not in member_colored:
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
    escher_map = _load_map(json_path)
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
                         dash_pattern="11.2,7",
                         member_boxes=True,
                         perNodeLabelSizes=None,
                         json_path=None,
                         member_colors=True,
                         member_groups=None,
                         group_colors=None,
                         member_legend=True):
    """Post-process a rendered Escher SVG.

    Drops specified label groups, rewrites fonts/colors for nodes, reactions,
    and segments, fades everything outside ``largeNodeLabels`` /
    ``largeEdgeLabels``, applies ``colorElements`` (a dict of hex color ->
    [element ids]) to highlight specific segments/nodes/reactions, dashes the
    segments in ``dashedEdges`` with ``dash_pattern`` (``"2,6"`` or similar
    gives dots), and draws a labelled box on each reaction anchor.

    ``member_colors`` gives every member its own colour, on its box and on
    all of its edges and arrowheads (:func:`color_member_edges`), from the
    palette in :mod:`escher_edit.palette`. Pass ``json_path`` with it so a
    member is recognised by name in every condition block of a stacked map,
    and a ``{member: colour}`` mapping instead of ``True`` to colour a series
    of maps alike; ``False`` leaves every edge in ``style.rxn_edge_hex`` —
    as for ``member_boxes``, the choice for a map whose reactions are not
    community members.
    Segments and reactions in ``colorElements`` keep their highlight colour.

    ``member_groups`` colours the members by group rather than one colour
    each — ``{member: group}`` or a function of the member, such as
    :func:`escher_edit.palette.taxon_groups` makes from a taxonomy to colour
    by phylum or genus — and ``group_colors`` (``{group: colour}``) fixes any
    group's colour, so a series of maps colours each group alike (see
    :func:`escher_edit.palette.member_colors`). The groups are then named in a
    legend below the map (:func:`draw_member_legend`), which ``member_legend``
    titles when it is a string and leaves out when ``False``. With a
    ``{member: colour}`` mapping, ``member_groups`` only labels the legend.

    ``member_boxes`` gates those boxes (:func:`draw_member_boxes`), which are
    what turns each community member from a bare vertex with a label floating
    over the edges into a node. It runs last because it moves the reaction
    labels into the boxes, out of the ``<g class="reaction">`` groups the
    earlier steps look them up in. With ``json_path``, a map from
    :mod:`escher_edit.build_map` supplies each box itself — the size its
    edges were laid out around, centred on the member's midmarker, labelled
    with the member (:func:`escher_edit.build_map.member_box_rects`); other
    maps get boxes sized to their labels. Turn it off for a map whose
    reactions are not members.

    ``json_path`` may also be the map already loaded, ``[header, body]``.

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

    escher_map = _load_map(json_path) if json_path is not None else None
    if perNodeLabelSizes is None and escher_map is not None:
        perNodeLabelSizes = extract_per_node_label_sizes_from_json(escher_map)

    remove_label_groups(soup, labels_to_remove)
    # before the labels are restyled, which can abbreviate them
    reaction_colors = (member_reaction_colors(soup, member_colors, escher_map,
                                              member_groups, group_colors)
                       if member_colors else {})
    members = (_reaction_members(soup, escher_map)
               if member_groups is not None else {})
    restyle_reaction_labels(soup, style, largeEdgeLabels, abbrevIDs)
    restyle_nodes(soup, style, largeNodeLabels, abbrevIDs,
                  perNodeLabelSizes=perNodeLabelSizes)
    restyle_segments(soup, style)
    dash_segments(soup, dashedEdges, dash_pattern)
    member_colored = (color_member_edges(soup, reaction_colors, style,
                                         largeEdgeLabels, colorElements)
                      if reaction_colors else set())
    apply_color_highlights(soup, colorElements, largeEdgeLabels, style,
                           member_colored)
    # last: this moves the reaction labels out of their reaction groups
    if member_boxes:
        # a reaction highlighted by colorElements keeps that colour on its box
        box_colors = dict(reaction_colors)
        for color, elements in (colorElements or {}).items():
            box_colors.update({element: f"#{color.lstrip('#')}"
                               for element in elements
                               if element.startswith("r")})
        rects = member_box_rects(escher_map) if escher_map is not None else {}
        names = ({f"r{key}": reaction.get("name")
                  for key, reaction in escher_map[1]["reactions"].items()}
                 if escher_map is not None else {})
        draw_member_boxes(soup, style, largeEdgeLabels, box_colors,
                          rects, names)
    if member_groups is not None and member_legend and reaction_colors:
        drawn = {members[rid]: color for rid, color in reaction_colors.items()}
        draw_member_legend(
            soup, group_legend(drawn, member_groups), style,
            title=member_legend if isinstance(member_legend, str) else None)

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
        # the notebook's highlight colours would collide with member hues
        member_colors=False,
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
