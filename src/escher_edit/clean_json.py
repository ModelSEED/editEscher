"""Clean up the Escher JSON map: swap metabolite node names to short
abbreviations and shorten reaction names for display.

Extracted from ABX_mouse_gut/Notebooks/escher_API_mapping.ipynb, section
"editing the SVG Escher Map" > "updating the JSON file".
"""
import logging
from json import load, dump
from pathlib import Path


log = logging.getLogger(__name__)

CONFIG_DIR = Path(__file__).parent / "config"

with (CONFIG_DIR / "cpd_abbrev.json").open() as _fh:
    CPD_ID_ABBRV = load(_fh)


def build_name_abbrev_table(msdb, output_csv="nameAbbrev.csv"):
    """Build a DataFrame of (name, abbreviation, compound ID) using a
    ModelSEED database instance (``modelseedpy.biochem.from_local``)."""
    from pandas import DataFrame

    nameAbbrev = DataFrame(columns=["name", "abbrev", "ID"])
    for cpdID, abbrev in CPD_ID_ABBRV.items():
        cpd = msdb.compounds.get_by_id(cpdID)
        nameAbbrev.loc[len(nameAbbrev)] = [cpd.name, abbrev, cpdID]
    if output_csv:
        nameAbbrev.to_csv(output_csv)
    return nameAbbrev


def cleanEscherJSON(escherPath, abbrev_map=CPD_ID_ABBRV):
    """Rewrite metabolite node names to short abbreviations and shorten
    reaction names, writing the result to ``*_cleaned0<suffix>``."""
    escherPath = Path(escherPath)
    with escherPath.open("r") as jsonIn:
        escherMap = load(jsonIn)

    # replace metabolite node names with their abbreviations; drop empties
    nodes_to_delete = []
    for nodeNum, content in escherMap[1]["nodes"].items():
        if content == {}:
            nodes_to_delete.append(nodeNum)
            continue
        if content["node_type"] == "metabolite":
            content["name"] = abbrev_map[content["bigg_id"]]

    for nodeNum in nodes_to_delete:
        escherMap[1]["nodes"].pop(nodeNum)

    # shorten reaction names (first six chars of the first dotted part,
    # plus the second dotted part)
    for content in escherMap[1]["reactions"].values():
        content["name"] = (content["name"].split(".")[0][:6]
                           + "." + content["name"].split(".")[1])

    out_path = escherPath.with_name(escherPath.stem + "_cleaned0" + escherPath.suffix)
    with out_path.open("w") as jsonOut:
        dump(escherMap, jsonOut, indent=3)
    log.info("wrote cleaned Escher JSON to %s", out_path)
    return out_path


def main():
    logging.basicConfig(level=logging.INFO)
    cleanEscherJSON("metabolite_focused_map_IDs.json")


if __name__ == "__main__":
    main()
