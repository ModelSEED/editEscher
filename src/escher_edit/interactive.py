"""Interactive HTML figures: point at part of a map to pick out what it touches.

A community map of any size is a tangle of edges, and in print the only way to
follow one is by eye. The HTML figure written beside every rendered SVG
(:func:`escher_edit.render.render_map_svg`) is the same drawing, made to
answer "what is this connected to?" directly:

- hovering an edge thickens it, lifts it above the others and enlarges the
  compound and the member box at its two ends, with their labels;
- hovering a compound, or its label, does that for every edge it has and
  every member at their other ends;
- hovering a member box does it for every edge of that member and every
  compound they reach;

and fades everything else back. A tooltip names what is under the pointer
and gives the flux. Clicking (or tapping, where there is no hover) pins the
highlight until the next click, and the address then links to it
(``figure.html#n36``). Scrolling zooms and dragging pans; enlarged nodes and
labels are brought up to a readable size on screen while only a few are lit,
so a highlight reads from the whole-map view without zooming in first.

An edge highlights its own two ends, not the path on through a cross-fed
compound — hover the compound for that. Nodes are per condition block, so a
highlight stays within its block.

Everything interactive is added here, to the HTML: the SVG figure is left as
it is for print, apart from the ``data-reaction`` its member boxes carry. The
page is self-contained — no scripts or fonts are fetched — so it opens
offline and can go into supplementary material as a single file.
"""
import html
import json
import logging
from pathlib import Path

from .svg_editor import _load_map

log = logging.getLogger(__name__)

_CSS = """
html, body { margin: 0; height: 100%; background: #ffffff; overflow: hidden; }
svg.escher-svg { display: block; width: 100vw; height: 100vh;
  touch-action: none; cursor: grab; user-select: none; -webkit-user-select: none; }
svg.escher-svg.panning { cursor: grabbing; }
svg.escher-svg * { pointer-events: none; }
svg.escher-svg .segment-hit { fill: none; stroke: transparent; stroke-width: 12px;
  vector-effect: non-scaling-stroke; pointer-events: stroke; cursor: pointer; }
svg.escher-svg .hover-target { fill: transparent; stroke: transparent; stroke-width: 8px;
  vector-effect: non-scaling-stroke; pointer-events: all; cursor: pointer; }
svg.escher-svg.dimmed :is(.segment-group, .node, #member-boxes > *,
  .reaction > .reaction-label-group):not(.hl) { opacity: 0.12; }
#escher-tip { position: fixed; z-index: 2; pointer-events: none; max-width: 340px;
  padding: 6px 9px; background: #ffffff; color: #222222; border: 1px solid #d4d4d4;
  border-radius: 6px; box-shadow: 0 2px 10px rgba(0, 0, 0, 0.12);
  font: 13px/1.4 system-ui, -apple-system, "Segoe UI", sans-serif; }
#escher-tip b { font-weight: 650; }
#escher-tip .muted { color: #666666; }
#escher-controls { position: fixed; z-index: 2; top: 10px; right: 10px; display: flex;
  gap: 4px; align-items: center; font: 12px/1.3 system-ui, -apple-system, "Segoe UI", sans-serif;
  color: #555555; }
#escher-controls span { margin-right: 6px; background: rgba(255, 255, 255, 0.85);
  padding: 3px 6px; border-radius: 4px; }
#escher-controls button { min-width: 30px; height: 30px; padding: 0 8px; font: inherit;
  font-size: 15px; color: #222222; background: #ffffff; border: 1px solid #cfcfcf;
  border-radius: 5px; cursor: pointer; }
#escher-controls button:hover { background: #f2f2f2; }
@media (max-width: 600px) { #escher-controls span { display: none; } }
"""

_JS = r"""
(function () {
  "use strict";
  var svg = document.querySelector("svg.escher-svg");
  var names = JSON.parse(document.getElementById("escher-names").textContent);
  var tip = document.getElementById("escher-tip");
  var reactionsLayer = svg.querySelector("#reactions");
  var nodesLayer = svg.querySelector("#nodes");
  var boxesLayer = svg.querySelector("#member-boxes");

  // On-screen floors, in CSS pixels, for what a highlight enlarges. Labels
  // are brought up to a readable size only while few are lit — a compound
  // fifty members eat would otherwise pile fifty enlarged boxes on each other;
  // the tooltip names them instead.
  var MIN_EDGE_PX = 4, MIN_RADIUS_PX = 6, MIN_LABEL_PX = 14;
  var READABLE_NODES = 12, READABLE_MEMBERS = 1;

  var edges = {}, byNode = {}, byMember = {};
  svg.querySelectorAll(".segment-group[data-member]").forEach(function (g) {
    var edge = {
      id: g.id, g: g, member: g.dataset.member,
      nodes: g.dataset.node ? g.dataset.node.split(" ") : [],
      flux: g.dataset.flux === undefined ? null : parseFloat(g.dataset.flux)
    };
    edges[edge.id] = edge;
    (byMember[edge.member] = byMember[edge.member] || []).push(edge);
    edge.nodes.forEach(function (n) { (byNode[n] = byNode[n] || []).push(edge); });
  });

  // What is pointed at is a fixed, invisible area per compound (circle and
  // label) and per member (box and label), laid over the drawing: the drawn
  // things grow when lit, and pointing at them directly would let a label
  // grow out from under the pointer and flicker. Compounds go under members,
  // as they are drawn; edges' own hit paths lie under both.
  (function () {
    var areas = {}, toUser = svg.getScreenCTM().inverse();
    svg.querySelectorAll("[data-key]:not(.segment-hit)").forEach(function (el) {
      var key = el.getAttribute("data-key"), matrix = el.getScreenCTM(), bbox;
      if (!matrix) return;
      try { bbox = el.getBBox(); } catch (err) { return; }
      matrix = toUser.multiply(matrix);
      [[bbox.x, bbox.y], [bbox.x + bbox.width, bbox.y],
       [bbox.x, bbox.y + bbox.height], [bbox.x + bbox.width, bbox.y + bbox.height]].forEach(function (corner) {
        var point = svg.createSVGPoint();
        point.x = corner[0]; point.y = corner[1];
        point = point.matrixTransform(matrix);
        var area = areas[key] || (areas[key] = [point.x, point.y, point.x, point.y]);
        area[0] = Math.min(area[0], point.x); area[1] = Math.min(area[1], point.y);
        area[2] = Math.max(area[2], point.x); area[3] = Math.max(area[3], point.y);
      });
    });
    var layer = document.createElementNS("http://www.w3.org/2000/svg", "g");
    layer.setAttribute("id", "hover-targets");
    Object.keys(areas).sort(function (a, b) {
      return (a.charAt(0) === "r") - (b.charAt(0) === "r");
    }).forEach(function (key) {
      var area = areas[key], rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      rect.setAttribute("class", "hover-target");
      rect.setAttribute("data-key", key);
      rect.setAttribute("x", area[0]); rect.setAttribute("y", area[1]);
      rect.setAttribute("width", area[2] - area[0]); rect.setAttribute("height", area[3] - area[1]);
      layer.appendChild(rect);
    });
    svg.appendChild(layer);
  })();

  function unique(list) { return list.filter(function (x, i) { return list.indexOf(x) === i; }); }

  function targets(key) {
    if (!key) return null;
    var kind = key.charAt(0), list;
    if (kind === "s" && edges[key]) {
      var edge = edges[key];
      return { edges: [edge], nodes: edge.nodes, members: [edge.member] };
    }
    if (kind === "n" && names[key] !== undefined) {
      list = byNode[key] || [];
      return { edges: list, nodes: [key], members: unique(list.map(function (e) { return e.member; })) };
    }
    if (kind === "r" && names[key] !== undefined) {
      list = byMember[key] || [];
      return { edges: list, members: [key],
               nodes: unique([].concat.apply([], list.map(function (e) { return e.nodes; }))) };
    }
    return null;
  }

  // Every change a highlight makes is recorded with how to undo it.
  var undo = [];
  function setStyle(el, prop, value) {
    // the attribute is restored verbatim, not re-serialised by the browser
    var had = el.hasAttribute("style"), old = el.getAttribute("style");
    el.style.setProperty(prop, value);
    undo.push(function () { if (had) el.setAttribute("style", old); else el.removeAttribute("style"); });
  }
  function setAttr(el, name, value) {
    var had = el.hasAttribute(name), old = el.getAttribute(name);
    el.setAttribute(name, value);
    undo.push(function () { if (had) el.setAttribute(name, old); else el.removeAttribute(name); });
  }
  function mark(el) {
    el.classList.add("hl");
    undo.push(function () { el.classList.remove("hl"); });
  }
  function lift(el, layer) {
    if (!layer || !el.parentNode) return;
    var placeholder = document.createComment("");
    el.parentNode.insertBefore(placeholder, el);
    layer.appendChild(el);
    undo.push(function () { placeholder.parentNode.insertBefore(el, placeholder); placeholder.remove(); });
  }
  function translateOf(el) {
    var found = /translate\(\s*([-\d.eE+]+)[\s,]+([-\d.eE+]+)/.exec(el.getAttribute("transform") || "");
    return found ? [parseFloat(found[1]), parseFloat(found[2])] : null;
  }

  function emphasizeEdge(edge, perUnit) {
    mark(edge.g);
    lift(edge.g, reactionsLayer);
    var factor = 1;
    edge.g.querySelectorAll("path").forEach(function (path) {
      if (path.closest(".arrowheads")) return;
      var computed = getComputedStyle(path);
      var width = parseFloat(computed.strokeWidth) || 1;
      var wider = Math.max(width * 2, MIN_EDGE_PX / perUnit);
      factor = wider / width;
      setStyle(path, "stroke-width", wider + "px");
      if (computed.strokeDasharray && computed.strokeDasharray !== "none") {
        setStyle(path, "stroke-dasharray", computed.strokeDasharray.split(/[\s,]+/)
          .map(function (d) { return parseFloat(d) * factor; }).join(" "));
      }
    });
    edge.g.querySelectorAll(".arrowhead").forEach(function (head) {
      setAttr(head, "transform", (head.getAttribute("transform") || "") + " scale(" + factor + ")");
    });
  }

  function emphasizeNode(id, perUnit, readable) {
    var g = document.getElementById(id);
    if (!g) return;
    mark(g);
    lift(g, nodesLayer);
    var circle = g.querySelector("circle"), label = g.querySelector("text");
    var radius = circle ? parseFloat(circle.getAttribute("r")) : 0, grown = radius;
    if (circle) {
      grown = Math.max(radius * 1.6, MIN_RADIUS_PX / perUnit);
      setAttr(circle, "r", grown);
    }
    if (label) {
      var size = parseFloat(getComputedStyle(label).fontSize) || 10;
      var larger = readable ? Math.max(size * 1.4, MIN_LABEL_PX / perUnit) : size * 1.4;
      setStyle(label, "font-size", larger + "px");
      var at = translateOf(label), centre = circle && translateOf(circle);
      if (at) {
        // out past the grown circle, and still centred on the same line
        var dx = centre ? Math.sign(at[0] - centre[0]) * (grown - radius) : 0;
        setAttr(label, "transform", "translate(" + (at[0] + dx) + "," + (at[1] + 0.35 * (larger - size)) + ")");
      }
    }
  }

  function emphasizeMember(id, perUnit, readable) {
    var parts = svg.querySelectorAll('[data-key="' + id + '"]');
    var rect = null, label = null;
    parts.forEach(function (el) {
      mark(el);
      if (el.tagName.toLowerCase() === "rect") rect = el;
      else if (!label) label = el.querySelector("text") || el;
    });
    if (!rect) return;              // a member without a box: its label stays lit
    var x = parseFloat(rect.getAttribute("x")), y = parseFloat(rect.getAttribute("y"));
    var cx = x + parseFloat(rect.getAttribute("width")) / 2;
    var cy = y + parseFloat(rect.getAttribute("height")) / 2;
    var size = label ? parseFloat(getComputedStyle(label).fontSize) || 10 : 10;
    var scale = readable ? Math.min(Math.max(1.2, MIN_LABEL_PX / (size * perUnit)), 5) : 1.15;
    var about = "translate(" + cx + "," + cy + ") scale(" + scale + ") translate(" + (-cx) + "," + (-cy) + ") ";
    parts.forEach(function (el) {
      setAttr(el, "transform", about + (el.getAttribute("transform") || ""));
      lift(el, boxesLayer);
    });
    setStyle(rect, "stroke", "#222222");
    setStyle(rect, "stroke-width", "2px");
    setStyle(rect, "vector-effect", "non-scaling-stroke");
  }

  function fmt(value) {
    var size = Math.abs(value);
    return size !== 0 && (size >= 1000 || size < 0.001) ? value.toExponential(2) : String(+value.toPrecision(3));
  }
  function listing(members) {
    var shown = unique(members).map(function (m) { return names[m]; });
    return shown.length > 6 ? shown.slice(0, 6).join(", ") + " and " + (shown.length - 6) + " more" : shown.join(", ");
  }
  function line(parts) {
    var div = document.createElement("div");
    parts.forEach(function (part) {
      var el = document.createElement(part[1] === "b" ? "b" : "span");
      if (part[1] === "muted") el.className = "muted";
      el.textContent = part[0];
      div.appendChild(el);
    });
    tip.appendChild(div);
  }
  function describe(key, found) {
    tip.textContent = "";
    var kind = key.charAt(0);
    if (kind === "s") {
      var edge = found.edges[0];
      var verb = edge.flux === null ? "exchanges" : edge.flux < 0 ? "consumes" : "produces";
      line([[names[edge.member], "b"], [" " + verb + " "], [edge.nodes.length ? names[edge.nodes[0]] : "", "b"]]);
      if (edge.flux !== null) line([["flux ", "muted"], [fmt(edge.flux)]]);
    } else if (kind === "n") {
      line([[names[key], "b"]]);
      var eaten = found.edges.filter(function (e) { return e.flux < 0; }).map(function (e) { return e.member; });
      var made = found.edges.filter(function (e) { return e.flux > 0; }).map(function (e) { return e.member; });
      if (made.length) line([["produced by ", "muted"], [listing(made)]]);
      if (eaten.length) line([["consumed by ", "muted"], [listing(eaten)]]);
    } else {
      line([[names[key], "b"]]);
      var consumes = found.edges.filter(function (e) { return e.flux < 0; }).length;
      var produces = found.edges.filter(function (e) { return e.flux > 0; }).length;
      line([["consumes ", "muted"], [String(consumes)], [" \u00b7 produces ", "muted"], [String(produces)]]);
    }
  }
  function placeTip(clientX, clientY) {
    if (tip.hidden) return;
    var pad = 14, w = tip.offsetWidth, h = tip.offsetHeight;
    var left = clientX + pad, top = clientY + pad;
    if (left + w > window.innerWidth - 4) left = Math.max(4, clientX - pad - w);
    if (top + h > window.innerHeight - 4) top = Math.max(4, clientY - pad - h);
    tip.style.left = left + "px";
    tip.style.top = top + "px";
  }

  var current = null, pinned = false, lastPointer = null;
  function show(key) {
    while (undo.length) undo.pop()();
    var found = targets(key);
    current = found ? key : null;
    svg.classList.toggle("dimmed", !!current);
    if (!current) { tip.hidden = true; return; }
    var perUnit = svg.getScreenCTM().a;
    found.edges.forEach(function (e) { emphasizeEdge(e, perUnit); });
    var readableNodes = found.nodes.length <= READABLE_NODES;
    var readableMembers = found.members.length <= READABLE_MEMBERS;
    found.nodes.forEach(function (n) { emphasizeNode(n, perUnit, readableNodes); });
    found.members.forEach(function (m) { emphasizeMember(m, perUnit, readableMembers); });
    describe(key, found);
    tip.hidden = false;
    if (lastPointer) placeTip(lastPointer[0], lastPointer[1]);
    else { tip.style.left = "10px"; tip.style.top = "10px"; }
  }
  function keyOf(target) {
    var el = target && target.closest ? target.closest(".segment-hit, .hover-target") : null;
    return el ? el.getAttribute("data-key") : null;
  }
  function pin(key) {
    pinned = !!targets(key);
    show(pinned ? key : null);
    try {
      history.replaceState(null, "", pinned ? "#" + key : location.pathname + location.search);
    } catch (err) { /* some viewers refuse history changes; the pin still holds */ }
  }

  // hover
  svg.addEventListener("mouseover", function (ev) {
    if (pinned || (drag && drag.moved)) return;
    var key = keyOf(ev.target);
    if (key !== current) show(key);
  });
  svg.addEventListener("mouseleave", function () { if (!pinned) show(null); });

  // pan, zoom, click to pin
  var initial = svg.getAttribute("viewBox").trim().split(/[\s,]+/).map(parseFloat);
  var home = { x: initial[0], y: initial[1], w: initial[2], h: initial[3] };
  var view = { x: home.x, y: home.y, w: home.w, h: home.h };
  var drag = null, refresh = 0;
  function applyView(rescale) {
    svg.setAttribute("viewBox", [view.x, view.y, view.w, view.h].join(" "));
    if (rescale && current && !refresh) {
      refresh = requestAnimationFrame(function () { refresh = 0; show(current); });
    }
  }
  function zoomAt(clientX, clientY, factor) {
    var point = svg.createSVGPoint();
    point.x = clientX; point.y = clientY;
    point = point.matrixTransform(svg.getScreenCTM().inverse());
    var width = Math.min(Math.max(view.w * factor, home.w / 400), home.w * 4);
    factor = width / view.w;
    view.x = point.x - (point.x - view.x) * factor;
    view.y = point.y - (point.y - view.y) * factor;
    view.w *= factor; view.h *= factor;
    applyView(true);
  }
  svg.addEventListener("wheel", function (ev) {
    ev.preventDefault();
    var step = ev.deltaMode === 1 ? 0.05 : ev.deltaMode === 2 ? 1 : 0.0015;
    zoomAt(ev.clientX, ev.clientY, Math.exp(ev.deltaY * step));
  }, { passive: false });
  svg.addEventListener("pointerdown", function (ev) {
    if (ev.button !== 0) return;
    drag = { id: ev.pointerId, x: ev.clientX, y: ev.clientY, moved: false,
             key: keyOf(ev.target), perUnit: svg.getScreenCTM().a,
             from: { x: view.x, y: view.y } };
  });
  window.addEventListener("pointermove", function (ev) {
    lastPointer = [ev.clientX, ev.clientY];
    if (drag && ev.pointerId === drag.id) {
      var dx = ev.clientX - drag.x, dy = ev.clientY - drag.y;
      if (!drag.moved && Math.sqrt(dx * dx + dy * dy) > 4) {
        drag.moved = true;
        svg.classList.add("panning");
      }
      if (drag.moved) {
        view.x = drag.from.x - dx / drag.perUnit;
        view.y = drag.from.y - dy / drag.perUnit;
        applyView(false);
        return;
      }
    }
    if (!pinned) placeTip(ev.clientX, ev.clientY);
  });
  window.addEventListener("pointerup", function (ev) {
    if (!drag || ev.pointerId !== drag.id) return;
    var ended = drag;
    drag = null;
    svg.classList.remove("panning");
    if (ended.moved) return;
    lastPointer = [ev.clientX, ev.clientY];
    if (pinned && ended.key === current) {
      pin(null);
      // a mouse is still over it, so it stays lit as a hover
      if (ev.pointerType === "mouse") show(ended.key);
    } else {
      pin(ended.key);
    }
  });
  window.addEventListener("keydown", function (ev) { if (ev.key === "Escape") pin(null); });

  document.getElementById("escher-zoom-in").addEventListener("click", function () {
    zoomAt(window.innerWidth / 2, window.innerHeight / 2, 1 / 1.5);
  });
  document.getElementById("escher-zoom-out").addEventListener("click", function () {
    zoomAt(window.innerWidth / 2, window.innerHeight / 2, 1.5);
  });
  document.getElementById("escher-fit").addEventListener("click", function () {
    view = { x: home.x, y: home.y, w: home.w, h: home.h };
    applyView(true);
  });

  function followHash() {
    var linked = decodeURIComponent(location.hash.slice(1));
    if (linked !== (pinned ? current : null)) {
      pinned = !!targets(linked);
      show(pinned ? linked : null);
    }
  }
  window.addEventListener("hashchange", followHash);
  if (location.hash.length > 1) followHash();
  svg.setAttribute("data-ready", "true");
})();
"""

_PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{css}</style>
</head>
<body>
{svg}
<div id="escher-controls">
<span>Hover or tap to highlight &middot; scroll to zoom &middot; drag to pan</span>
<button type="button" id="escher-zoom-in" aria-label="Zoom in">+</button>
<button type="button" id="escher-zoom-out" aria-label="Zoom out">&minus;</button>
<button type="button" id="escher-fit">Fit</button>
</div>
<div id="escher-tip" role="status" hidden></div>
<script type="application/json" id="escher-names">{names}</script>
<script>{js}</script>
</body>
</html>
"""


def _visible(path):
    return "display:none" not in path.get("style", "").replace(" ", "")


def annotate_map_svg(soup, escher_map):
    """Tie the drawing in ``soup`` to the map's topology for
    :func:`write_interactive_html`.

    Each segment group gets ``data-member`` (its reaction's id), ``data-node``
    (the metabolite node it reaches) and ``data-flux``; a transparent copy of
    each drawn edge, wide on screen at any zoom, goes into a
    ``#segment-hits`` layer under everything else — so where an edge passes
    under a compound, a label or a member box, pointing there picks the
    thing on top — to be pointed at; and every
    hoverable thing — hit edge, metabolite node, member box and its label —
    gets the ``data-key`` the page highlights by (the page lays a fixed
    hover area over each compound and member from these). Returns ``{key: name}`` for
    the metabolite nodes and members, for the tooltips.
    """
    body = escher_map[1]
    nodes = body["nodes"]
    names, segments = {}, {}
    for key, reaction in body["reactions"].items():
        rid = f"r{key}"
        names[rid] = reaction.get("name") or reaction.get("bigg_id", "")
        fluxes = {met["bigg_id"]: met["coefficient"]
                  for met in reaction.get("metabolites", [])}
        for segment_key, segment in reaction["segments"].items():
            ends = [node_id for node_id in (segment["from_node_id"], segment["to_node_id"])
                    if nodes.get(node_id, {}).get("node_type") == "metabolite"]
            flux = fluxes.get(nodes[ends[0]].get("bigg_id")) if ends else None
            segments[f"s{segment_key}"] = (rid, [f"n{node_id}" for node_id in ends], flux)
    for key, node in nodes.items():
        if node.get("node_type") == "metabolite":
            name, bigg_id = node.get("name"), node.get("bigg_id", "")
            names[f"n{key}"] = (f"{name} ({bigg_id})" if name and bigg_id and name != bigg_id
                                else name or bigg_id)

    hits = soup.new_tag("g")
    hits["id"] = "segment-hits"
    for group in soup.find_all("g", class_="segment-group"):
        sid = group.get("id")
        if sid not in segments:
            continue
        rid, ends, flux = segments[sid]
        group["data-member"] = rid
        if ends:
            group["data-node"] = " ".join(ends)
        if flux is not None:
            group["data-flux"] = f"{flux:.6g}"
        for path in group.find_all("path", recursive=False):
            if _visible(path) and path.get("d"):
                hit = soup.new_tag("path")
                hit["class"] = "segment-hit"
                hit["d"] = path["d"]
                hit["data-key"] = sid
                hits.append(hit)
    layer = soup.find(id="reactions")
    if layer is not None:
        layer.insert_before(hits)
    else:
        soup.svg.append(hits)

    for group in soup.find_all("g", class_="node"):
        if group.get("id") in names:
            group["data-key"] = group["id"]
    for element in soup.find_all(attrs={"data-reaction": True}):
        element["data-key"] = element["data-reaction"]
    for reaction in soup.find_all("g", class_="reaction"):
        label = reaction.find("g", class_="reaction-label-group")
        if label is not None and reaction.get("id") in names:
            label["data-key"] = reaction["id"]
    return names


def write_interactive_html(svg_path, escher_map, out_path=None, title=None):
    """Write the interactive HTML figure for a finished map SVG.

    ``svg_path`` is a figure from :func:`escher_edit.render.render_map_svg`
    (or :func:`escher_edit.svg_editor.EscherSVG_processing` given the map),
    ``escher_map`` the map it was drawn from, as ``[header, body]`` or a path.
    ``out_path`` defaults to the SVG's path with an ``.html`` suffix and
    ``title`` to the map's caption when it has exactly one, else the SVG's
    stem. See the module docstring for what the page does.

    Returns the path written.
    """
    from bs4 import BeautifulSoup
    svg_path = Path(svg_path)
    out_path = Path(out_path) if out_path is not None else svg_path.with_suffix(".html")
    soup = BeautifulSoup(svg_path.read_text(encoding="utf-8"), "lxml-xml")
    if soup.svg is None:
        raise ValueError(f"{svg_path} holds no <svg>")
    escher_map = _load_map(escher_map)
    names = annotate_map_svg(soup, escher_map)
    if title is None:
        captions = [label.get("text") for label in escher_map[1].get("text_labels", {}).values()]
        title = captions[0] if len(captions) == 1 and captions[0] else svg_path.stem

    root = soup.svg
    for size in ("width", "height"):
        if size in root.attrs:
            del root[size]
    if not root.get("viewBox"):
        canvas = soup.find(id="canvas")
        if canvas is None:
            raise ValueError(f"{svg_path} has neither a viewBox nor a canvas "
                             f"to size the page by")
        root["viewBox"] = " ".join(canvas.get(k, "0") for k in ("x", "y", "width", "height"))
    root["preserveAspectRatio"] = "xMidYMid meet"
    # prettify() wraps text in newlines, which HTML would keep as a space
    for text in soup.find_all("text"):
        if text.string is not None:
            text.string = text.string.strip()

    page = _PAGE.format(
        title=html.escape(str(title)),
        css=_CSS,
        svg=str(root),
        names=json.dumps(names, ensure_ascii=False).replace("</", "<\\/"),
        js=_JS)
    out_path.write_text(page, encoding="utf-8")
    log.info("wrote interactive figure %s", out_path)
    return out_path
