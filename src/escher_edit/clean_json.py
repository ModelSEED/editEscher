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


def _shorten_reaction_name(name):
    """Shorten a reaction name using the original string-slice convention:
    ``name.split(".")[0][:6] + "." + name.split(".")[1]``. Returns the name
    unchanged when it does not contain a "." (otherwise the slice would
    raise IndexError)."""
    if "." not in name:
        return name
    parts = name.split(".")
    return parts[0][:6] + "." + parts[1]


def cleanEscherJSON(
    escherPath,
    abbrev_map=CPD_ID_ABBRV,
    rxn_abbrev_map=None,
    rewrite_bigg_id=True,
):
    """Rewrite metabolite node names to short abbreviations and shorten
    reaction names, writing the result to ``*_cleaned0<suffix>``.

    Parameters
    ----------
    escherPath : str or Path
        Path to the Escher JSON map.
    abbrev_map : dict[str, str]
        Mapping from metabolite ``bigg_id`` to the short abbreviation to use
        as the new ``name`` (and ``bigg_id`` when ``rewrite_bigg_id`` is
        True). Metabolites whose ``bigg_id`` is absent from the map are
        left unchanged.
    rxn_abbrev_map : dict[str, str] or None
        Optional mapping from reaction ``bigg_id`` to a short display label.
        When provided, matching reactions have their ``name`` (and
        ``bigg_id`` when ``rewrite_bigg_id`` is True) replaced. When None,
        reaction names fall back to the original string-slice shortener
        (safe for names without a "." — those are left unchanged).
    rewrite_bigg_id : bool
        Escher displays a node's ``bigg_id`` as the on-map label (``name``
        is only tooltip metadata), so rewriting ``bigg_id`` is required for
        the new labels to appear. Disable to preserve original IDs (e.g.,
        when the map must remain matchable to a specific model).
    """
    escherPath = Path(escherPath)
    with escherPath.open("r") as jsonIn:
        escherMap = load(jsonIn)

    # Replace metabolite node names (and optionally bigg_ids); drop empties
    nodes_to_delete = []
    for nodeNum, content in escherMap[1]["nodes"].items():
        if content == {}:
            nodes_to_delete.append(nodeNum)
            continue
        if content["node_type"] == "metabolite":
            original_bigg = content.get("bigg_id")
            if original_bigg in abbrev_map:
                abbrev = abbrev_map[original_bigg]
                content["name"] = abbrev
                if rewrite_bigg_id:
                    content["bigg_id"] = abbrev

    for nodeNum in nodes_to_delete:
        escherMap[1]["nodes"].pop(nodeNum)

    # Replace reaction names (and optionally bigg_ids). Update the
    # metabolite bigg_ids inside each reaction's metabolites list so the
    # reaction remains internally consistent with the rewritten nodes.
    for content in escherMap[1]["reactions"].values():
        original_rxn_bigg = content.get("bigg_id")
        if rxn_abbrev_map and original_rxn_bigg in rxn_abbrev_map:
            content["name"] = rxn_abbrev_map[original_rxn_bigg]
            if rewrite_bigg_id:
                content["bigg_id"] = rxn_abbrev_map[original_rxn_bigg]
        else:
            content["name"] = _shorten_reaction_name(content["name"])
        if rewrite_bigg_id:
            for met in content.get("metabolites", []):
                if met.get("bigg_id") in abbrev_map:
                    met["bigg_id"] = abbrev_map[met["bigg_id"]]

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
