"""Clean up the Escher JSON map: swap metabolite node names to short
abbreviations and shorten reaction names for display.

Extracted from ABX_mouse_gut/Notebooks/escher_API_mapping.ipynb, section
"editing the SVG Escher Map" > "updating the JSON file".
"""
from json import load, dump


# Mapping from ModelSEED compound IDs to short abbreviations used in the
# Escher map labels.
CPD_ID_ABBRV = {
    "cpd00076": "sucr", "cpd00141": "prpa", "cpd00211": "butr",
    "cpd00130": "malat", "cpd00137": "citr", "cpd00024": "akg",
    "cpd00036": "succ", "cpd00064": "ornth", "cpd03847": "myrst",
    "cpd00106": "fumr", "cpd00108": "galct", "cpd00382": "melit",
    "cpd00105": "ribs", "cpd00121": "inost", "cpd00751": "fucos",
    "cpd00122": "acglum", "cpd00158": "cellb",
    "cpd00214": "palm", "cpd03198": "melib", "cpd01171": "dulco",
    "cpd00020": "pyr", "cpd00224": "arbns", "cpd00082": "fru",
    "cpd00851": "4hpro", "cpd01107": "decac", "cpd00396": "rhmn",
    "cpd01055": "allos", "cpd00027": "glu",
    "cpd01242": "2drib",
}


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
    reaction names, writing the result to ``*_cleaned0.json``."""
    with open(escherPath, "r") as jsonIn:
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

    out_path = escherPath.replace(".json", "_cleaned0.json")
    with open(out_path, "w") as jsonOut:
        dump(escherMap, jsonOut, indent=3)
    return out_path


if __name__ == "__main__":
    cleanEscherJSON("metabolite_focused_map_IDs.json")
