"""Reverse the direction of reactions in an Escher JSON map.

Escher renders arrow direction based on two pieces of state per reaction:

* the sign of each entry in ``reaction["metabolites"][i]["coefficient"]``
  (negative = reactant, positive = product);
* the ``from_node_id`` / ``to_node_id`` on every segment in
  ``reaction["segments"]``.

Reversing a reaction flips the signs of all coefficients *and* swaps
``from_node_id`` with ``to_node_id`` on every segment. This mirrors the
behaviour of the "reverse reaction" control in Escher's interactive editor.
"""
import logging
from json import load, dump
from pathlib import Path


log = logging.getLogger(__name__)


def reverse_reaction(rxn_content):
    """Reverse a single reaction dict in place.

    Parameters
    ----------
    rxn_content : dict
        A reaction entry from ``escherMap[1]["reactions"]``.
    """
    for met in rxn_content.get("metabolites", []):
        met["coefficient"] = -met["coefficient"]
    for seg in rxn_content.get("segments", {}).values():
        seg["from_node_id"], seg["to_node_id"] = seg["to_node_id"], seg["from_node_id"]


def reverse_reactions_in_map(escher_path, rxn_bigg_ids, output_path=None):
    """Reverse every reaction whose ``bigg_id`` appears in ``rxn_bigg_ids``.

    Parameters
    ----------
    escher_path : str or Path
        Path to the Escher JSON map to read.
    rxn_bigg_ids : iterable[str]
        Reaction ``bigg_id`` values to reverse.
    output_path : str or Path or None
        Where to write the modified map. Defaults to
        ``<stem>_reversed<suffix>`` alongside the input.

    Returns
    -------
    pathlib.Path
        The output path written.
    """
    escher_path = Path(escher_path)
    with escher_path.open("r") as f:
        escher_map = load(f)

    targets = set(rxn_bigg_ids)
    reversed_count = 0
    for rxn_content in escher_map[1]["reactions"].values():
        if rxn_content.get("bigg_id") in targets:
            reverse_reaction(rxn_content)
            reversed_count += 1

    if output_path is None:
        output_path = escher_path.with_name(
            escher_path.stem + "_reversed" + escher_path.suffix
        )
    output_path = Path(output_path)
    with output_path.open("w") as f:
        dump(escher_map, f, indent=3)
    log.info("reversed %d reaction(s); wrote %s", reversed_count, output_path)
    return output_path


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("escher_path", help="Input Escher JSON map")
    parser.add_argument(
        "rxn_bigg_ids", nargs="+",
        help="One or more reaction bigg_id values to reverse",
    )
    parser.add_argument("-o", "--output", default=None, help="Output path")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    reverse_reactions_in_map(args.escher_path, args.rxn_bigg_ids, args.output)


if __name__ == "__main__":
    main()
