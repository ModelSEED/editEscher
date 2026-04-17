"""Edit a rendered Escher SVG: restyle labels/edges, fade non-highlighted
elements, draw ASV node rectangles, and recolor selected segments.

Extracted from ABX_mouse_gut/Notebooks/escher_API_mapping.ipynb, section
"editing the SVG Escher Map". The ``shapely``-based label-layout helpers are
preserved from the notebook but were never enabled in the final call.
"""
import random
import re

# shapely is only required if the label-layout helpers below are used.
# from shapely.geometry import Point, box


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


def EscherSVG_processing(svg_path="metabolite_focused_map.svg",
                         labels_to_remove=("stoichiometry-labels",),
                         node_label_px=10, rxn_label_px=30,
                         node_label_hex="e0865b", rxn_label_hex="202078",
                         rxn_edge_hex="5aa6df", rxn_edge_px=10,
                         colorElements=None, largeNodeLabels=None,
                         largeEdgeLabels=None, tint_factor=0.5,
                         stroke_reduction=100, abbrevIDs=None, msdb=None,
                         large_node_offset=100, dashedEdges=None):
    """Post-process a rendered Escher SVG.

    Drops specified label groups, rewrites fonts/colors for nodes, reactions,
    and segments, fades everything outside ``largeNodeLabels`` /
    ``largeEdgeLabels``, draws small rectangles to mark short straight
    segments as ASV anchor nodes, applies ``colorElements`` (a dict of
    hex color -> [element ids]) to highlight specific segments/nodes/
    reactions, and dashes the segments in ``dashedEdges``.

    The edited SVG is written to ``*_edited.svg``.
    """
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(open(svg_path, "r").read(), "lxml-xml")

    # remove all stoichiometric or other undesirable labels
    for label in labels_to_remove:
        stoich_labels = soup.find_all(class_=label)
        for lbl in stoich_labels:
            lbl.decompose()

    # label font-size and color change
    metadata = soup.find("style", attrs={"type": "text/css"})

    # reactions
    metadata.string = metadata.string.replace(
        ".reaction-label{font-size:30px;fill:#202078;text-rendering:optimizelegibility}",
        ".reaction-label{" + f"font-size:{rxn_label_px}px;fill:#{rxn_label_hex};stroke:#000000;stroke-width:{rxn_label_px / stroke_reduction}px;" + "text-rendering:optimizelegibility}")
    if largeEdgeLabels is not None:
        fadedEdgesColor = tint_color(rxn_edge_hex, tint_factor)
        fadedEdgesLabelsColor = tint_color(rxn_label_hex, tint_factor)
        metadata.string = metadata.string + f" .m{fadedEdgesLabelsColor}-label" + "{" + f"font-size:{rxn_label_px / 1.5}px;fill:#{fadedEdgesLabelsColor};stroke:#000000;stroke-width:{rxn_label_px / stroke_reduction}px;" + "text-rendering:optimizelegibility}"
    for reaction in soup.find_all(class_="reaction"):
        label = reaction.find("text")
        if largeEdgeLabels is not None and reaction["id"] not in largeEdgeLabels:
            label["class"] = f"m{fadedEdgesLabelsColor}-label"
        if abbrevIDs is not None:
            cleanedStr = label.string.split("_")[0].replace("-ABX", "").replace("RC", "").replace("WD", "")
            label.string = cleanedStr.split(".")[0][:6] + "." + cleanedStr.split(".")[1]

    # nodes
    metadata.string = metadata.string.replace(
        ".node-label{font-size:20px}",
        ".node-label{" + f"font-size:{node_label_px}px;fill:#{node_label_hex};stroke:#000000;stroke-width:{node_label_px / stroke_reduction}px;" + "}")
    metadata.string = metadata.string.replace(
        ".metabolite-circle{stroke:#a24510;fill:#e0865b}",
        ".metabolite-circle{stroke:#" + f"{node_label_hex};fill:#{node_label_hex}" + "}")
    if largeNodeLabels is not None:
        fadedNodesColor = tint_color(node_label_hex, tint_factor)
        for node in soup.find_all(class_="node"):
            circle = node.find("circle")
            label = node.find("text")
            if label is None:
                continue
            if node["id"] not in largeNodeLabels:
                if any(x is None for x in [circle, label]):
                    continue
                circle["class"] = f"m{fadedNodesColor}"
                label["class"] = f"m{fadedNodesColor}-label"
                if abbrevIDs is not None:
                    cpdID = label.string
                    label.string = abbrevIDs[cpdID]
                    print(node["id"])

        metadata.string = metadata.string + f" .m{fadedNodesColor}" + "{" + f"stroke:#{fadedNodesColor};fill:#{fadedNodesColor}" + "}"
        metadata.string = metadata.string + f" .m{fadedNodesColor}-label" + "{" + f"font-size:{node_label_px / 1.2}px;fill:#{fadedNodesColor};stroke:#000000;stroke-width:{node_label_px / stroke_reduction}px;" + "}"

    # segments
    metadata.string = metadata.string.replace(
        ".segment{" + "stroke:#334E75;stroke-width:10px;fill:none}",
        ".segment{" + f"stroke:#{rxn_edge_hex};" + f"stroke-width:{rxn_edge_px}px;" + "fill:none}")

    # miscellaneous text
    metadata.string = metadata.string.replace(
        ".text-label-input{font-size:50px}",
        ".text-label-input{font-size:" + f"{rxn_label_px * 1.8}" + "px}")

    # add ASV nodes (rectangles at the start of short straight segments)
    for edge in soup.find_all(class_="segment-group"):
        if dashedEdges is not None:
            if edge["id"] in dashedEdges:
                edge["stroke-dasharray"] = "8,5"
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
            rect["fill"] = f"#{rxn_label_hex}"
            soup.svg.append(rect)

    # recolor select edges/nodes/reactions per colorElements
    if colorElements is not None:
        edgeColor = {}
        for color, lst in colorElements.items():
            edgeColor.update({x: color for x in lst})
        for edge in soup.find_all(class_="segment-group"):
            if edge["id"] in edgeColor:
                color = edgeColor[edge["id"]]
                for child in edge.find_all(class_="segment"):
                    child["class"] = f"m{color}"
            elif largeEdgeLabels is not None:
                for child in edge.find_all(class_="segment"):
                    child["class"] = f"m{fadedEdgesColor}"
        for color, elements in colorElements.items():
            metadata.string = metadata.string + f" .m{color}" + "{" + f"stroke:#{color};" + f"stroke-width:{rxn_edge_px * 1.2}px;" + "fill:none}"
            if any("r" in x for x in elements):
                for ele in elements:
                    if "r" in ele:
                        reaction = soup.find(attrs={"id": ele})
                        label = reaction.find("text")
                        label["class"] = f"m{color}-label"
                metadata.string = metadata.string + f" .m{color}-label" + "{" + f"font-size:{rxn_label_px}px;fill:#{color};stroke:#000000;stroke-width:{rxn_label_px / stroke_reduction}px;" + "text-rendering:optimizelegibility;font-family:sans-serif;font-style:italic;font-weight:700;}"
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
            metadata.string = metadata.string + f" .m{fadedEdgesColor}" + "{" + f"stroke:#{fadedEdgesColor};" + f"stroke-width:{rxn_edge_px / 1.2}px;" + "fill:none}"

    out_path = svg_path.replace(".svg", "_edited.svg")
    with open(out_path, "w", encoding="utf-8") as file:
        file.write(soup.prettify())
    return out_path


if __name__ == "__main__":
    # The original invocation from the notebook. ``consumptionEdges`` comes
    # from escher_model_mapping.build_direction_tracking, and ``msdb`` from
    # ``modelseedpy.biochem.from_local("../../ModelSEEDDatabase")``.
    from escher_model_mapping import build_direction_tracking

    _, consumptionEdges = build_direction_tracking()

    msdb = None  # optionally: from modelseedpy.biochem import from_local; msdb = from_local("../../ModelSEEDDatabase")

    EscherSVG_processing(
        "metabolite_focused_map_IDs_cleaned.svg",
        node_label_px=90, rxn_label_px=90, rxn_edge_px=6,
        node_label_hex="E36400", rxn_label_hex="cc0000", rxn_edge_hex="2C79BD",
        largeNodeLabels=[
            "n202", "n311", "n900", "n963", "n33", "n32", "n830", "n71",
            "n1322", "n1443", "n1672", "n1673",
            "n2066", "n2102", "n206", "n404", "n873", "n1090", "n1182",
            "n1746", "n2154", "n2154", "n1338",  # RC
            "n539", "n1063", "n1282", "n1511", "n1576"  # cellibiose
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
        msdb=msdb, tint_factor=0.8, stroke_reduction=100,
        dashedEdges=consumptionEdges,
    )
