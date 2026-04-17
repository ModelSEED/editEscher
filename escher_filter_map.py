"""Filter an Escher JSON map to drop low-coefficient and excluded metabolites.

Extracted from ABX_mouse_gut/Notebooks/escher_API_mapping.ipynb, section
"filtering the Escher map for fewer reactions based on the metabolite fluxes".

The original notebook block loads ``metabolite_focused_map.json``, writes two
lookup files (``Escher_nodeNum_cpdIDs.json`` and ``EscherNodeMapping.json``),
then emits a filtered map at ``metabolite_focused_map0.json``. The logic is
preserved here as reusable functions.
"""
import logging
from json import load, dump
from pathlib import Path


log = logging.getLogger(__name__)

CONFIG_DIR = Path(__file__).parent / "config"

with (CONFIG_DIR / "filter.json").open() as _fh:
    _filter_cfg = load(_fh)

AA_NAMES = list(_filter_cfg["aa_names"])
DEFAULT_SKIP_NAMES = AA_NAMES + list(_filter_cfg["skip_extras"])


def build_node_lookups(escher_map,
                       node_num_path="Escher_nodeNum_cpdIDs.json",
                       name_path="EscherNodeMapping.json"):
    """Build and persist the nodeNum->biggID and biggID->name lookups."""
    nodes = escher_map[1]["nodes"]
    node_id = {num: content["bigg_id"]
               for num, content in nodes.items() if "bigg_id" in content}
    id_name = {content["bigg_id"]: content["name"]
               for content in nodes.values() if "bigg_id" in content}

    if node_num_path:
        with open(node_num_path, "w") as jsonOut:
            dump(node_id, jsonOut, indent=3)
    if name_path:
        with open(name_path, "w") as jsonOut:
            dump(id_name, jsonOut, indent=3)
    return node_id, id_name


def filter_escher_map(input_path="metabolite_focused_map.json",
                      output_path="metabolite_focused_map0.json",
                      skip_names=DEFAULT_SKIP_NAMES,
                      coefficient_threshold=0.05):
    """Filter the Escher JSON map, dropping nodes that fall below the
    coefficient threshold or that belong to ``skip_names``.

    Metabolite nodes whose absolute flux is below ``coefficient_threshold`` are
    removed, as are any nodes whose display name is in ``skip_names`` (amino
    acids plus a few default metabolites like Niacin and Acetate). All other
    reactions and map elements are passed through unchanged.
    """
    with open(input_path, "r") as jsonIn:
        escherMap = load(jsonIn)

    og_escherMap = [ele for ele in escherMap]
    nodeID, idName = build_node_lookups(og_escherMap)

    newEscher = []
    for index, ele in enumerate(og_escherMap):
        if index == 0:
            newEscher = [ele]
        newEscher.append({})
        nodeNumCoef = {}

        # pass reactions through and accumulate per-node coefficients
        for k, val in ele.items():
            if k != "reactions":
                continue
            newEscher[index][k] = val
            for rxnNum, content in val.items():
                # map metabolite bigg_ids to coefficients
                metCoefs = {met["bigg_id"]: met["coefficient"]
                            for met in content["metabolites"]}
                # collect nodeNum -> bigg_id referenced by this reaction's segments
                nodeNumIDs = {seg["to_node_id"]: nodeID[seg["to_node_id"]]
                              for seg in content["segments"].values()
                              if seg["to_node_id"] in nodeID}
                nodeNumIDs.update({seg["from_node_id"]: nodeID[seg["from_node_id"]]
                                   for seg in content["segments"].values()
                                   if seg["from_node_id"] in nodeID})
                nodeNumCoef.update({num: metCoefs[bid]
                                    for num, bid in nodeNumIDs.items()
                                    if bid in metCoefs})
                newEscher[index][k][rxnNum] = content

        # filter nodes
        for k, val in ele.items():
            if k == "nodes":
                new_nodes = {}
                for nodeNum, content in val.items():
                    new_nodes[nodeNum] = {}
                    if nodeNum in nodeNumCoef:
                        if abs(nodeNumCoef[nodeNum]) < coefficient_threshold:
                            continue
                        if idName.get(content.get("bigg_id", "")) in skip_names:
                            continue
                    for k2, v in content.items():
                        if k2 == "node_is_primary":
                            new_nodes[nodeNum][k2] = True
                        else:
                            new_nodes[nodeNum][k2] = v
                newEscher[index][k] = new_nodes
            elif k != "reactions":
                newEscher[index][k] = val

    with open(output_path, "w") as jsonOut:
        dump(newEscher, jsonOut, indent=3)
    log.info("wrote filtered Escher map to %s", output_path)
    return newEscher


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    filter_escher_map()
