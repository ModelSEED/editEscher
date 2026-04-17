"""Build the model<->SVG direction-tracking map and categorize Escher segments.

Extracted from ABX_mouse_gut/Notebooks/escher_API_mapping.ipynb, the block
that builds ``modelSVG_mapping.json`` and classifies segments as reactant/
product and nodes as consumed/produced/intermediate.
"""
from json import load, dump

from escher_clean_json import CPD_ID_ABBRV


DEFAULT_DIETS = ["RC", "WD"]
DEFAULT_DAYS = ["-1.5", "1.5", "4.0", "6.0", "9.0", "12.5"]


def build_direction_tracking(model_path="ASVInteractionModel.json",
                             escher_path="metabolite_focused_map_IDs_cleaned.json",
                             output_path="modelSVG_mapping.json",
                             diets=DEFAULT_DIETS,
                             days=DEFAULT_DAYS,
                             abbrev_map=CPD_ID_ABBRV):
    """Build per-(diet, day, reaction) mapping of reactants/products and the
    SVG segment IDs that feed each. Returns ``(directionTracking,
    consumptionEdges)`` and persists ``directionTracking`` to JSON.
    """
    with open(model_path, "r") as jsonIn:
        model = load(jsonIn)

    directionTracking = {diet: {day: {} for day in days} for diet in diets}

    for content in model["reactions"]:
        nameParts = content["id"].split("-ABX")
        diet = nameParts[0][-2:]
        nameParts[0] = nameParts[0][:-2]
        day = nameParts[1].split("_")[1]
        simpleName = (nameParts[0].split(".")[0][:6]
                      + "." + nameParts[0].split(".")[1])
        directionTracking[diet][day][simpleName] = {
            "ID": "",
            "segments": {"consumption": {}, "production": {}},
            "reactants": [],
            "products": [],
        }
        for ID, flux in content["metabolites"].items():
            flux = float(flux)
            if flux < 0:
                directionTracking[diet][day][simpleName]["reactants"].append(ID)
            elif flux > 0:
                directionTracking[diet][day][simpleName]["products"].append(ID)

    with open(escher_path, "r") as jsonIn:
        escherMap = load(jsonIn)

    consumptionEdges = []
    for rxnNum, content in escherMap[1]["reactions"].items():
        nameParts = content["bigg_id"].split("-ABX")
        diet = nameParts[0][-2:]
        day = nameParts[1].split("_")[1]
        rxnName = content["name"]
        directionTracking[diet][day][rxnName]["ID"] = f"r{rxnNum}"
        for segNum, seg in content["segments"].items():
            fromNode, toNode = seg["from_node_id"], seg["to_node_id"]
            fromID = escherMap[1]["nodes"][fromNode].get("bigg_id", "")
            toID = escherMap[1]["nodes"][toNode].get("bigg_id", "")
            if fromID != "":
                print(fromID)
            if toID in directionTracking[diet][day][rxnName]["reactants"]:
                directionTracking[diet][day][rxnName]["segments"]["consumption"][f"s{segNum}"] = (
                    f"n{toNode}", abbrev_map[toID])
                consumptionEdges.append(f"s{segNum}")
            if toID in directionTracking[diet][day][rxnName]["products"]:
                directionTracking[diet][day][rxnName]["segments"]["production"][f"s{segNum}"] = (
                    f"n{toNode}", abbrev_map[toID])

    if output_path:
        with open(output_path, "w") as jsonOut:
            dump(directionTracking, jsonOut, indent=3)

    return directionTracking, consumptionEdges


def categorize_segments_and_nodes(escher_path="metabolite_focused_map_IDs_cleaned.json"):
    """Walk the Escher JSON and classify each segment as reactant/product
    edge and each metabolite node as consumed/produced/intermediate."""
    with open(escher_path, "r") as jsonIn:
        escherMap = load(jsonIn)

    nodeCategories = {"producedNodes": set(), "consumedNodes": set(),
                      "intermediateNodes": set()}
    segment_categorization = {"reactant": [], "product": []}

    for rxnNum, content in escherMap[1]["reactions"].items():
        for segNum, seg in content["segments"].items():
            fromNode, toNode = seg["from_node_id"], seg["to_node_id"]
            from_type = escherMap[1]["nodes"][fromNode]["node_type"]
            to_type = escherMap[1]["nodes"][toNode]["node_type"]

            if "marker" in to_type and from_type == "metabolite":
                segment_categorization["reactant"].append(segNum)
                if fromNode in nodeCategories["intermediateNodes"]:
                    pass
                elif fromNode not in nodeCategories["producedNodes"]:
                    nodeCategories["consumedNodes"].add(fromNode)
                elif fromNode in nodeCategories["producedNodes"]:
                    nodeCategories["producedNodes"].remove(fromNode)
                    nodeCategories["intermediateNodes"].add(fromNode)
            elif "marker" in from_type and to_type == "metabolite":
                segment_categorization["product"].append(segNum)
                if toNode in nodeCategories["intermediateNodes"]:
                    pass
                elif toNode not in nodeCategories["consumedNodes"]:
                    nodeCategories["producedNodes"].add(toNode)
                elif toNode in nodeCategories["consumedNodes"]:
                    nodeCategories["consumedNodes"].remove(toNode)
                    nodeCategories["intermediateNodes"].add(toNode)

    return segment_categorization, nodeCategories


if __name__ == "__main__":
    build_direction_tracking()
    categorize_segments_and_nodes()
