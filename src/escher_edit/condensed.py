"""Condensed interaction map: a slide-sized, bipartite view of one community.

The full community map (:mod:`escher_edit.build_map`) draws every member's net
exchange with the shared pool, collapsing each compound to one node and
routing the cross-fed ones through lanes between the members. That is the map
to read on screen; on a slide its labels are a few points high. The condensed
map is the companion for a slide or a figure panel. It was first made, ad hoc,
for a reactor-community presentation, and is offered here as an optional
aesthetic of the same data:

* **condensed** — only the compounds asked for (or all), and for each compound
  and side only the members with the largest ``|flux|``, taken in order until
  they cover ``cover`` of that side's total (:func:`condense_fluxes`);
* **bipartite** — what members take up sits in a column on the left, what
  they excrete in a column on the right, so a compound both made and used
  appears once on each side (:func:`split_sides`); no lanes;
* **large type** — member labels about 10 pt and compound labels about 9 pt
  on a figure about 8.5 inches wide, on a landscape canvas
  (:class:`CondensedStyle`);
* **edges from the box outline** — each edge leaves the face of its member's
  box nearest its compound and runs to it as an S-curve that is horizontal at
  both ends; an excreted compound gets an arrowhead, an uptake none;
* **edge width = |flux|** on one log scale (:attr:`CondensedStyle.width_scale`),
  and the edges of compounds that some members make and others use dashed
  (member-to-member), the rest solid (medium / effluent);
* **one key on the right** — the member (or group) colours, the width scale,
  what dashed and solid mean, and any short notes; a header over each column.

    from escher_edit.condensed import render_condensed_map
    render_condensed_map({"memberA": {"cpd00029": -0.4, "cpd00011": 0.8}, ...},
                         "panel.svg", compound_names=names,
                         compounds=["cpd00029", "cpd00011"], title="Day 560")

    escher-edit-condensed ASVMetaboliteInteractions.csv --condition RC-ABX_12.5 \\
        --id-csv ids.csv --names-csv names.csv --compounds cpd00029 cpd00011 \\
        --cover 0.6 -o panel.svg

The map itself is an ordinary Escher map (:func:`build_condensed_map`), so it
also opens in Escher; the condensed look is in the SVG this module renders.
"""
import argparse
import logging
import math
import re
import textwrap
from collections import defaultdict
from pathlib import Path

from .build_map import (MapStyle, build_escher_map, build_member_reactions,
                        load_compound_names, parse_interaction_matrix)
from .palette import group_legend, member_colors as _member_colors, text_color_on
from .render import add_member_color_arguments, member_color_options, render_map_svg
from .svg_editor import EscherStyle

log = logging.getLogger(__name__)

INK = "#222222"
UPTAKE, EXCRETION = "taken up", "excreted"
#: Escher collapses nodes by ``bigg_id`` and labels them with it, so each
#: side's node is keyed by the compound's display name, the excretion side's
#: with this mark after it (a no-break space: invisible, but a different id).
OUT_MARK = "\u00a0"
_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
_NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
_TRANSLATE_Y = re.compile(r"translate\(\s*[-\d.eE+]+[\s,]+([-\d.eE+]+)")


class CondensedStyle:
    """Sizes of the condensed map, in map units (the SVG's px).

    ``member_label_px`` / ``node_label_px`` are the member-box and compound
    label sizes; on a figure 8.5 in wide drawn from a canvas about 3,000 units
    wide, the defaults print at about 10 pt and 9 pt. The layout reserves
    room for labels a little larger than drawn (``reserve_member`` /
    ``reserve_node``) so bold or italic glyphs never touch. ``width_scale``
    maps ``|flux|`` to edge width on a log scale, ``((flux0, px0), (flux1,
    px1))``, clamped at both ends; ``key_fluxes`` are the widths shown in the
    key and ``flux_units`` its unit line.
    """

    def __init__(self, member_label_px=56.0, node_label_px=50.0, member_pitch=140.0,
                 reserve_member=1.2, reserve_node=1.15,
                 width_scale=((1e-5, 2.5), (1.0, 20.0)),
                 key_fluxes=(1e-4, 1e-3, 1e-2, 1e-1), flux_units=""):
        F, N = float(member_label_px), float(node_label_px)
        self.F, self.node_px = F, N
        self.box_pad = 0.4 * F
        self.layout = MapStyle(
            member_pitch=member_pitch, mixed_lane_dx=250.0, lane_dx_fraction=0.04,
            column_clearance=4 * N, min_node_spacing=130.0, lane_node_spacing=170.0,
            lane_span_fraction=1.0, edge_curve="s", metabolite_radius=34.0,
            crossfed_radius=26.0, label_font_px=N * reserve_node,
            reaction_label_font_px=F * reserve_member, block_label_font_px=2.0 * F,
            block_label_gap=4.0 * F, member_node_pad=self.box_pad,
            member_label_gap=1.5 * F, max_aspect=None, canvas_margin=F)
        self.escher = EscherStyle(node_label_px=N, rxn_label_px=F, rxn_edge_px=6.0,
                                  member_box_pad=self.box_pad, member_box_edge_px=3.0,
                                  member_box_radius=10.0)
        self.width_scale = width_scale
        self.key_fluxes = tuple(key_fluxes)
        self.flux_units = flux_units
        self.legend_px = 0.8 * F
        self.margin = F

    def edge_width(self, flux):
        """Edge width (px) for ``|flux|`` on the log scale, clamped."""
        (f0, w0), (f1, w1) = self.width_scale
        t = ((math.log10(max(abs(flux), 1e-300)) - math.log10(f0))
             / (math.log10(f1) - math.log10(f0)))
        return w0 + min(max(t, 0.0), 1.0) * (w1 - w0)


# ------------------------------------------------------------------ the data


def condense_fluxes(fluxes_by_member, compounds=None, cover=0.6):
    """Keep, per compound and side, the largest exchanges covering ``cover``.

    ``fluxes_by_member`` is ``{member: {compound: flux}}`` (negative = taken
    up, positive = excreted). For each compound in ``compounds`` (default:
    every compound) and each side, members are taken in order of ``|flux|``
    until they account for at least ``cover`` of that side's total; ``cover``
    = 1 keeps every exchange of the chosen compounds.

    Returns ``(condensed, coverage)``: the kept ``{member: {compound: flux}}``
    and ``{(compound, side): {"kept", "of", "share"}}`` with ``side`` one of
    :data:`UPTAKE` / :data:`EXCRETION` — say what the panel leaves out in its
    caption.
    """
    if not 0 < cover <= 1:
        raise ValueError(f"cover must be in (0, 1]: {cover!r}")
    every = sorted({c for fluxes in fluxes_by_member.values() for c in fluxes})
    chosen = list(compounds) if compounds else every
    missing = [c for c in chosen if c not in every]
    if missing:
        log.warning("compounds no member exchanges: %s", ", ".join(missing))
    kept, coverage = defaultdict(dict), {}
    for cpd in chosen:
        for sign, side in ((-1, UPTAKE), (1, EXCRETION)):
            entries = sorted(((abs(f[cpd]), m) for m, f in fluxes_by_member.items()
                              if f.get(cpd, 0.0) * sign > 0), reverse=True)
            if not entries:
                continue
            total = sum(a for a, _ in entries)
            done, n = 0.0, 0
            for a, member in entries:
                if done >= cover * total:
                    break
                kept[member][cpd] = fluxes_by_member[member][cpd]
                done += a
                n += 1
            coverage[(cpd, side)] = {"kept": n, "of": len(entries), "share": done / total}
    return dict(kept), coverage


def crossfed_compounds(fluxes_by_member):
    """Compounds some member takes up and another excretes (drawn dashed)."""
    up, out = set(), set()
    for fluxes in fluxes_by_member.values():
        for cpd, flux in fluxes.items():
            (up if flux < 0 else out).add(cpd)
    return up & out


def split_sides(fluxes_by_member, compound_names=None):
    """Give each compound one node per side, keyed by its display name.

    The uptake node is keyed by the name and the excretion node by the name
    plus :data:`OUT_MARK`; two compounds sharing a name keep their ids in
    brackets. :func:`escher_edit.build_map.classify_compounds` then puts every
    uptake node in the input column and every excretion node in the output
    column, which is the bipartite layout.

    Returns ``(split_fluxes, sides)``: ``{member: {key: flux}}`` and ``{key:
    (compound, side)}`` with ``side`` :data:`UPTAKE` or :data:`EXCRETION`.
    """
    compound_names = compound_names or {}
    used = defaultdict(set)
    for fluxes in fluxes_by_member.values():
        for cpd in fluxes:
            used[compound_names.get(cpd, cpd)].add(cpd)
    label = {cpd: (name if len(cpds) == 1 else f"{name} ({cpd})")
             for name, cpds in used.items() for cpd in cpds}
    split, sides = {}, {}
    for member, fluxes in fluxes_by_member.items():
        row = {}
        for cpd, flux in fluxes.items():
            if not flux:
                continue
            side = UPTAKE if flux < 0 else EXCRETION
            key = label[cpd] + ("" if side == UPTAKE else OUT_MARK)
            row[key] = flux
            sides[key] = (cpd, side)
        if row:
            split[member] = row
    return split, sides


def build_condensed_map(fluxes_by_member, compound_names=None, compounds=None, cover=0.6,
                        title="", style=None, min_abs_flux=0.0, crossfed=None,
                        map_name="condensed_interaction_map"):
    """Build the condensed, bipartite Escher map of one community.

    ``fluxes_by_member`` is ``{member: {compound: flux}}``; exchanges below
    ``min_abs_flux`` are dropped first. ``compounds`` and ``cover`` choose what
    is drawn (:func:`condense_fluxes`). ``crossfed`` is the set of compounds
    whose edges are dashed; by default the compounds some member takes up and
    another excretes in the whole (uncondensed) input. ``title`` labels the
    block.

    Returns ``(escher_map, info)``; ``info`` holds the ``segments`` (``{"s<key>":
    {"member", "compound", "flux", "rid", "dashed"}}``), the ``coverage`` and
    the ``style`` — pass it to :func:`render_condensed_svg`.
    """
    style = style or CondensedStyle()
    full = {m: {c: f for c, f in fluxes.items() if f and abs(f) >= min_abs_flux}
            for m, fluxes in fluxes_by_member.items()}
    crossfed = crossfed_compounds(full) if crossfed is None else set(crossfed)
    condensed, coverage = condense_fluxes(full, compounds, cover)
    if not condensed:
        raise ValueError("nothing to draw: no member exchanges the chosen compounds")
    split, sides = split_sides(condensed, compound_names)
    names = {key: key for key in sides}
    members = build_member_reactions(split, model_id="", compound_names=names)
    escher_map = build_escher_map([(title, members)], compound_names=names,
                                  map_name=map_name, style=style.layout)
    body = escher_map[1]
    nodes = body["nodes"]
    segments = {}
    for rkey, rxn in body["reactions"].items():
        coef = {met["bigg_id"]: met["coefficient"] for met in rxn["metabolites"]}
        for skey, seg in rxn["segments"].items():
            for nid in (seg["from_node_id"], seg["to_node_id"]):
                node = nodes.get(nid, {})
                if node.get("node_type") == "metabolite" and node.get("bigg_id") in coef:
                    cpd, side = sides[node["bigg_id"]]
                    segments[f"s{skey}"] = {"member": rxn["name"], "compound": cpd, "side": side,
                                            "flux": coef[node["bigg_id"]],
                                            "rid": f"r{rkey}", "dashed": cpd in crossfed}
    return escher_map, {"segments": segments, "coverage": coverage, "sides": sides,
                        "crossfed": sorted(crossfed), "style": style}


# ------------------------------------------------------------------ the SVG


def _set_style(tag, **props):
    decl = {}
    for item in tag.get("style", "").split(";"):
        if ":" in item:
            k, v = item.split(":", 1)
            decl[k.strip()] = v.strip()
    decl.update({k.replace("_", "-"): str(v) for k, v in props.items()})
    tag["style"] = ";".join(f"{k}:{v}" for k, v in decl.items())


def _text(soup, x, y, string, px, weight="400", anchor=None, fill=INK):
    t = soup.new_tag("text")
    t["x"], t["y"] = f"{x:.1f}", f"{y:.1f}"
    if anchor:
        t["text-anchor"] = anchor
    t.string = string
    _set_style(t, font_family="sans-serif", font_size=f"{px:g}px", fill=fill, stroke="none",
               font_weight=weight)
    return t


def _pt(p):
    return f"{p[0]:.2f},{p[1]:.2f}"


def dash_pattern(width):
    """(dash, gap) for a dashed edge of this width, scaled with the width."""
    d = max(1.0, width / 6.0)
    return round(11.2 * d, 1), round(7 * d, 1)


def anchor_edges_to_faces(soup, segments, F):
    """Move each edge's member end onto the outline of its member box.

    The member end goes to the box face nearest the compound (left face for a
    compound left of the box centre, else the right face), at a height spread
    over the middle 60 % of the face in the order of the compound ends (so the
    edges on one face do not cross beside the box); the curve is rebuilt as an
    S with horizontal tangents at both ends, and an arrowhead at the compound
    is turned to point along it. Only single-cubic segments are moved. Returns
    the number of edges moved.
    """
    rects = {r.get("data-reaction"): tuple(float(r[k]) for k in ("x", "y", "width", "height"))
             for r in soup.find_all("rect", class_="member-box")}
    groups = defaultdict(list)
    for group in soup.find_all("g", class_="segment-group"):
        info = segments.get(group.get("id"))
        if not info or info["rid"] not in rects:
            continue
        rx, ry, rw, rh = rects[info["rid"]]
        cx, cy = rx + rw / 2, ry + rh / 2
        for path in group.find_all("path", class_="segment"):
            if "display:none" in path.get("style", ""):
                continue
            nums = [float(v) for v in _NUM.findall(path.get("d", ""))]
            if len(nums) != 8 or not re.match(r"\s*M[^MLCSQZ]*C[^MLCSQZ]*$", path["d"]):
                continue
            pts = list(zip(nums[0::2], nums[1::2]))
            member_first = ((pts[0][0] - cx) ** 2 + (pts[0][1] - cy) ** 2
                            <= (pts[3][0] - cx) ** 2 + (pts[3][1] - cy) ** 2)
            comp = pts[3] if member_first else pts[0]
            face = "L" if comp[0] < cx else "R"
            groups[(info["rid"], face)].append((path, member_first, comp, group))
    moved = 0
    for (rid, face), items in groups.items():
        rx, ry, rw, rh = rects[rid]
        fx = rx if face == "L" else rx + rw
        out = -1.0 if face == "L" else 1.0
        items.sort(key=lambda it: it[2][1])
        k = len(items)
        for i, (path, member_first, comp, group) in enumerate(items):
            y = ry + rh / 2 if k == 1 else ry + rh * (0.2 + 0.6 * i / (k - 1))
            p0 = (fx, y)
            reach = max(abs(comp[0] - fx) * 0.5, 1.2 * F)
            c1, c2 = (fx + out * reach, y), (comp[0] - out * reach, comp[1])
            path["d"] = (f"M{_pt(p0)} C{_pt(c1)} {_pt(c2)} {_pt(comp)}" if member_first
                         else f"M{_pt(comp)} C{_pt(c2)} {_pt(c1)} {_pt(p0)}")
            for head in group.find_all("path", class_="arrowhead"):
                m = re.match(r"\s*translate\(\s*([-\d.eE+]+)[\s,]+([-\d.eE+]+)\)\s*rotate\(\s*[-\d.eE+]+\)(.*)$",
                             head.get("transform", ""))
                if m and abs(float(m.group(1)) - comp[0]) < 1.0 and abs(float(m.group(2)) - comp[1]) < 1.0:
                    # the arrow shape points along +y: rotate 270 to point right, 90 to point left
                    head["transform"] = (f"translate({m.group(1)},{m.group(2)})"
                                         f"rotate({270.0 if out > 0 else 90.0}){m.group(3)}")
            moved += 1
    return moved


def _canvas_resize(soup, canvas, left=0.0, top=0.0, right=0.0, bottom=0.0):
    """Grow the canvas and the root viewBox / size by the given margins."""
    canvas["x"] = f"{float(canvas['x']) - left:g}"
    canvas["y"] = f"{float(canvas['y']) - top:g}"
    canvas["width"] = f"{float(canvas['width']) + left + right:g}"
    canvas["height"] = f"{float(canvas['height']) + top + bottom:g}"
    root = soup.svg
    if root.get("viewBox"):
        vx, vy, vw, vh = map(float, root["viewBox"].split())
        root["viewBox"] = f"{vx - left:g} {vy - top:g} {vw + left + right:g} {vh + top + bottom:g}"
    for attr, extra in (("width", left + right), ("height", top + bottom)):
        if root.get(attr):
            root[attr] = f"{float(root[attr]) + extra:g}"


def _extent(soup):
    """(top, bottom) of the drawn boxes, nodes and edges."""
    ys = []
    for r in soup.find_all("rect", class_="member-box"):
        ys += [float(r["y"]), float(r["y"]) + float(r["height"])]
    for c in soup.find_all("circle"):
        m = _TRANSLATE_Y.search(c.get("transform", ""))
        if m:
            ys += [float(m.group(1)) - float(c.get("r", 0)), float(m.group(1)) + float(c.get("r", 0))]
    for p in soup.find_all("path", class_="segment"):
        ys += [float(y) for y in re.findall(r"[-\d.eE+]+,([-\d.eE+]+)", p.get("d", ""))]
    return min(ys), max(ys)


def _columns(escher_map, sides):
    """x of the uptake and excretion columns (mean x of each side's nodes)."""
    xs = defaultdict(list)
    for node in escher_map[1]["nodes"].values():
        if node.get("node_type") == "metabolite" and node.get("bigg_id") in sides:
            xs[sides[node["bigg_id"]][1]].append(node["x"])
    return {side: sum(v) / len(v) for side, v in xs.items() if v}


def render_condensed_svg(escher_map, info, out_path, member_colors=None, legend=None,
                         legend_title="Members", headers=(UPTAKE, EXCRETION), notes=(),
                         dash_labels=("member-to-member", "medium / effluent"), html=False):
    """Render a map from :func:`build_condensed_map` in the condensed style.

    ``member_colors`` is ``{member: "#rrggbb"}`` (default: one distinct colour
    each, :func:`escher_edit.palette.member_colors`); ``legend`` the key's
    ``[(label, "#rrggbb")]`` entries (default: one per member, as coloured).
    ``headers`` are the two column headers (``None`` for none), ``notes`` short
    lines set under the key, ``dash_labels`` the key's words for dashed and
    solid edges. ``html`` also writes the interactive page
    (:func:`escher_edit.interactive.write_interactive_html`).

    Returns the path of the SVG written.
    """
    from bs4 import BeautifulSoup
    from .interactive import write_interactive_html

    style = info["style"]
    out_path = Path(out_path)
    names = sorted({r["name"] for r in escher_map[1]["reactions"].values()})
    if member_colors is None:
        member_colors = _member_colors(names)
    if legend is None:
        legend = [(n, member_colors[n]) for n in names if n in member_colors]
    segments = info["segments"]
    dashed = [sid for sid, s in segments.items() if s["dashed"]]
    render_map_svg(escher_map, out_path, dashed=False, dashedEdges=dashed, layout=style.layout,
                   html=False, member_colors=member_colors, member_legend=False, style=style.escher)

    F, N, L = style.F, style.node_px, style.legend_px
    soup = BeautifulSoup(out_path.read_text(encoding="utf-8"), "lxml-xml")
    anchor_edges_to_faces(soup, segments, F)
    # edge width = |flux|; dashes scaled with the width; arrowheads with the edge
    for group in soup.find_all("g", class_="segment-group"):
        seg = segments.get(group.get("id"))
        if not seg:
            continue
        width = style.edge_width(seg["flux"])
        for path in group.find_all("path", class_="segment"):
            _set_style(path, stroke_width=f"{width:.2f}px")
        for head in group.find_all("path", class_="arrowhead"):
            head["transform"] = head.get("transform", "") + f" scale({min(max(width / 10.0, 1.0), 1.6):.2f})"
        if group.get("stroke-dasharray"):
            group["stroke-dasharray"] = "{:.1f},{:.1f}".format(*dash_pattern(width))
    # compound labels on a white halo; member labels centred in their boxes
    css = soup.find("style", attrs={"type": "text/css"})
    if css is not None and css.string:
        css.string = css.string + (" svg.escher-svg .node-label{paint-order:stroke;stroke:#ffffff;"
                                   f"stroke-width:{0.22 * N:.1f}px;stroke-linejoin:round}}")
    boxes = soup.find(id="member-boxes")
    if boxes is not None:
        rects = {r.get("data-reaction"): r for r in boxes.find_all("rect", class_="member-box")}
        for label in boxes.find_all("g", class_="reaction-label-group"):
            rect = rects.get(label.get("data-reaction"))
            if rect is None:
                continue
            cx = float(rect["x"]) + float(rect["width"]) / 2
            cy = float(rect["y"]) + float(rect["height"]) / 2
            label["transform"] = f"translate({cx:.3f},{cy + 0.35 * F:.3f})"
            text = label.find("text")
            if text is not None:
                text["text-anchor"] = "middle"
                fill = rect.get("fill") or ""
                if re.fullmatch(r"#?[0-9a-fA-F]{6}", fill):
                    _set_style(text, fill=text_color_on(fill.lstrip("#")))
    canvas = soup.find(id="canvas")
    top, bottom = _extent(soup)
    # column headers, a row above the drawing (the canvas grows up if need be)
    if headers:
        cols = _columns(escher_map, info["sides"])
        base = top - 0.9 * N
        title = soup.select_one("#text-labels text")
        if title is not None:
            m = re.search(r"translate\(\s*([-\d.eE+]+)[\s,]+([-\d.eE+]+)", title.get("transform", ""))
            if m and float(m.group(2)) + 0.6 * N > base - 0.9 * N:
                shift = float(m.group(2)) + 0.6 * N - (base - 0.9 * N)
                title["transform"] = f"translate({m.group(1)},{float(m.group(2)) - shift:.1f})"
        g = soup.new_tag("g")
        g["id"] = "column-headers"
        for side, text in zip((UPTAKE, EXCRETION), headers):
            if side in cols and text:
                g.append(_text(soup, cols[side], base, text, N, weight="700", anchor="middle",
                               fill="#4a4a4a"))
        soup.svg.append(g)
    lows = [float(canvas["y"])]
    for t in soup.find_all("text"):
        m = _TRANSLATE_Y.search(t.get("transform", "") or "")
        y = float(t["y"]) if t.get("y") is not None else (float(m.group(1)) if m else None)
        if y is not None:
            lows.append(y - 1.2 * F)
    room_up = max(0.0, float(canvas["y"]) - (min(lows) - style.margin / 2))
    if room_up:
        _canvas_resize(soup, canvas, top=room_up)
    # the key, in a column right of the map
    x0 = float(canvas["x"]) + float(canvas["width"])
    y = float(canvas["y"]) + style.margin + 1.2 * L
    key = soup.new_tag("g")
    key["id"] = "condensed-key"
    widest = 0.0
    if legend:
        key.append(_text(soup, x0, y, legend_title, L, weight="700"))
        widest = 0.58 * L * len(legend_title)
    for label, color in legend:
        y += 1.45 * L
        sw = soup.new_tag("rect")
        sw["class"] = "legend-swatch"
        sw["x"], sw["y"] = f"{x0:.1f}", f"{y - 0.82 * L:.1f}"
        sw["width"] = sw["height"] = f"{L:.1f}"
        sw["rx"] = "8"
        sw["fill"] = color
        key.append(sw)
        key.append(_text(soup, x0 + 1.4 * L, y, label, L))
        widest = max(widest, 1.4 * L + 0.58 * L * len(label))
    y += 1.9 * L
    key.append(_text(soup, x0, y, "Edge width = |flux|", L, weight="700"))
    if style.flux_units:
        y += 1.0 * L
        key.append(_text(soup, x0, y, style.flux_units, 0.8 * L))
    run = 2.4 * L
    for flux in style.key_fluxes:
        y += 1.4 * L
        p = soup.new_tag("path")
        p["d"] = f"M{x0:.1f},{y - 0.35 * L:.1f} L{x0 + run:.1f},{y - 0.35 * L:.1f}"
        _set_style(p, stroke=INK, stroke_width=f"{style.edge_width(flux):.2f}px", fill="none")
        key.append(p)
        exp = math.log10(flux)
        label = ("10" + str(round(exp)).translate(_SUP)) if abs(exp - round(exp)) < 1e-9 else f"{flux:g}"
        key.append(_text(soup, x0 + run + 0.4 * L, y, label, L))
    for is_dashed, words in zip((True, False), dash_labels):
        if not words:
            continue
        y += 1.4 * L
        width = style.edge_width(1e-2)
        p = soup.new_tag("path")
        r = run
        if is_dashed:
            dash, gap = dash_pattern(width)
            n = max(1, round((run + gap) / (dash + gap)))
            r = n * (dash + gap) - gap
            p["stroke-dasharray"] = f"{dash:.1f},{gap:.1f}"
        p["d"] = f"M{x0:.1f},{y - 0.35 * L:.1f} L{x0 + r:.1f},{y - 0.35 * L:.1f}"
        _set_style(p, stroke=INK, stroke_width=f"{width:.2f}px", fill="none")
        key.append(p)
        line = ("dashed: " if is_dashed else "solid: ") + words
        key.append(_text(soup, x0 + run + 0.4 * L, y, line, L))
        widest = max(widest, run + 0.4 * L + 0.58 * L * len(line))
    for note in notes:
        for i, part in enumerate(textwrap.wrap(note, 30)):
            y += (1.5 if i == 0 else 1.15) * L
            key.append(_text(soup, x0, y, part, 0.85 * L))
            widest = max(widest, 0.55 * 0.85 * L * len(part))
    soup.svg.append(key)
    _canvas_resize(soup, canvas, right=widest + style.margin)
    need = max(y + L, bottom) + style.margin - (float(canvas["y"]) + float(canvas["height"]))
    if need > 0:
        _canvas_resize(soup, canvas, bottom=need)
    out_path.write_text(str(soup), encoding="utf-8")
    log.info("condensed map %s: %d edges, %d dashed", out_path, len(segments), len(dashed))
    if html:
        write_interactive_html(out_path, escher_map,
                               out_path.with_suffix(".html") if html is True else html)
    return out_path


def render_condensed_map(fluxes_by_member, out_path, compound_names=None, compounds=None,
                         cover=0.6, title="", style=None, min_abs_flux=0.0, crossfed=None,
                         member_colors=None, member_groups=None, group_colors=None,
                         legend_title=None, **render):
    """:func:`build_condensed_map` then :func:`render_condensed_svg`, in one call.

    With ``member_groups`` (``{member: group}``) members are coloured by group
    and the key lists the groups (:func:`escher_edit.palette.group_legend`);
    otherwise each member gets its own colour unless ``member_colors`` is
    given. Returns ``(svg_path, escher_map, info)``.
    """
    escher_map, info = build_condensed_map(fluxes_by_member, compound_names, compounds, cover,
                                           title, style, min_abs_flux, crossfed)
    names = sorted({r["name"] for r in escher_map[1]["reactions"].values()})
    legend = render.pop("legend", None)
    if render.pop("member_legend", True) is False:
        legend = []                                       # --no-legend: the key keeps only the edge scale
    if member_groups is not None:
        member_colors = _member_colors(names, groups=member_groups, group_colors=group_colors)
        if legend is None:
            legend = group_legend(member_colors, member_groups)
    svg = render_condensed_svg(escher_map, info, out_path, member_colors=member_colors,
                               legend=legend, legend_title=legend_title or (
                                   "Member group" if member_groups is not None else "Members"),
                               **render)
    return svg, escher_map, info


def main():
    parser = argparse.ArgumentParser(
        description="Draw a condensed, bipartite interaction map of one condition: the chosen "
                    "compounds, each side's largest exchanges, edges from the box outline, "
                    "edge width = |flux|, one key on the right.")
    parser.add_argument("csv", help="Wide interaction matrix (ASVMetaboliteInteractions.csv)")
    parser.add_argument("--condition", default=None,
                        help="The condition (model id) to draw; default the only one in the file")
    parser.add_argument("--id-csv", default=None, help="Compound id CSV (with --names-csv)")
    parser.add_argument("--names-csv", default=None, help="Compound names CSV (with --id-csv)")
    parser.add_argument("--compounds", nargs="+", default=None,
                        help="Compound ids to draw (default: all)")
    parser.add_argument("--cover", type=float, default=0.6,
                        help="Per compound and side, draw the largest exchanges until they cover "
                             "this share of the side's total (default 0.6; 1 draws all)")
    parser.add_argument("--min-abs-flux", type=float, default=0.0,
                        help="Ignore exchanges below this |flux| first")
    parser.add_argument("--title", default=None, help="Title (default: the condition)")
    parser.add_argument("--flux-units", default="", help="Unit line under the edge-width key")
    parser.add_argument("--note", action="append", default=[], help="A note under the key (repeatable)")
    parser.add_argument("--html", action="store_true", help="Also write the interactive page")
    parser.add_argument("-o", "--output", default=None, help="Output SVG (default: <condition>_condensed.svg)")
    add_member_color_arguments(parser)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    data = parse_interaction_matrix(args.csv)
    if args.condition is None:
        if len(data) != 1:
            parser.error(f"the file has {len(data)} conditions; choose one with --condition: "
                         + ", ".join(sorted(data)))
        args.condition = next(iter(data))
    if args.condition not in data:
        parser.error(f"no condition {args.condition!r}; the file has: " + ", ".join(sorted(data)))
    if (args.id_csv is None) != (args.names_csv is None):
        parser.error("--id-csv and --names-csv go together")
    names = load_compound_names(args.id_csv, args.names_csv) if args.id_csv else None
    coloring = member_color_options(args, parser)
    out = args.output or f"{args.condition}_condensed.svg"
    render_condensed_map(data[args.condition], out, compound_names=names, compounds=args.compounds,
                         cover=args.cover, title=args.title if args.title is not None else args.condition,
                         style=CondensedStyle(flux_units=args.flux_units),
                         min_abs_flux=args.min_abs_flux, member_groups=coloring.get("member_groups"),
                         group_colors=coloring.get("group_colors"),
                         legend_title=coloring.get("member_legend") if isinstance(
                             coloring.get("member_legend"), str) else None,
                         notes=args.note, html=args.html,
                         member_legend=coloring.get("member_legend", True))
    print(out)


if __name__ == "__main__":
    main()
