"""Layout utilities for Escher JSON maps.

Functions here operate on the ``escher_map`` list structure
(``[metadata, content]``) and mutate ``content["nodes"]``,
``content["reactions"]``, and ``content["text_labels"]`` in place. All
helpers accept an already-loaded map so callers can chain them without
re-reading the file.

Typical use alongside ``cleanEscherJSON``: after abbreviations replace long
metabolite and reaction names, the original (verbose) label positions are
usually too far from their anchors. ``normalize_metabolite_label_offsets``
restores a consistent dx/dy offset from each node, and
``set_reaction_label_positions`` moves each reaction label to an explicit
(x, y). Custom free-form annotations left over from the verbose design can
be cleared with ``remove_text_labels``.
"""
import logging
from json import load, dump
from pathlib import Path


log = logging.getLogger(__name__)


def normalize_metabolite_label_offsets(escher_map, dx=15, dy=0):
    """Set ``label_x = node.x + dx`` and ``label_y = node.y + dy`` for every
    metabolite node. Mutates ``escher_map`` in place.
    """
    for node in escher_map[1]["nodes"].values():
        if node.get("node_type") == "metabolite":
            if "x" in node and "y" in node:
                node["label_x"] = node["x"] + dx
                node["label_y"] = node["y"] + dy


def set_node_label_positions(escher_map, positions_by_node_id):
    """Move metabolite labels to explicit absolute coordinates.

    Use this when multiple nodes share the same ``bigg_id`` (e.g., a
    cofactor drawn in several places) and per-node positioning is required.

    Parameters
    ----------
    escher_map : list
        Escher map structure.
    positions_by_node_id : dict[str, tuple[float, float]]
        Mapping from node id (the key in ``escherMap[1]["nodes"]``) to
        ``(label_x, label_y)``.
    """
    nodes = escher_map[1]["nodes"]
    for node_id, (x, y) in positions_by_node_id.items():
        if node_id in nodes and nodes[node_id].get("node_type") == "metabolite":
            nodes[node_id]["label_x"] = x
            nodes[node_id]["label_y"] = y


def set_reaction_label_positions(escher_map, positions):
    """Move reaction labels to explicit coordinates, keyed by reaction
    ``bigg_id``.

    Parameters
    ----------
    escher_map : list
        Escher map structure.
    positions : dict[str, tuple[float, float]]
        Mapping from reaction ``bigg_id`` to ``(label_x, label_y)``.
    """
    for rxn in escher_map[1]["reactions"].values():
        bigg = rxn.get("bigg_id")
        if bigg in positions:
            x, y = positions[bigg]
            rxn["label_x"] = x
            rxn["label_y"] = y


def set_reaction_label_positions_by_id(escher_map, positions_by_reaction_id):
    """Move reaction labels to explicit coordinates, keyed by the reaction's
    key in ``escherMap[1]["reactions"]`` (unique even when ``bigg_id``
    duplicates exist).
    """
    rxns = escher_map[1]["reactions"]
    for rxn_id, (x, y) in positions_by_reaction_id.items():
        if rxn_id in rxns:
            rxns[rxn_id]["label_x"] = x
            rxns[rxn_id]["label_y"] = y


def set_node_label_sizes(escher_map, label_size, bigg_ids=None, node_ids=None):
    """Write a (non-standard) ``label_size`` attribute onto selected
    metabolite nodes.

    Escher's stock web renderer ignores unknown fields, but downstream SVG
    post-processors (e.g., ``svg_editor.restyle_nodes``) can read this
    attribute and emit an inline ``font-size`` style per node. Useful for
    emphasizing specific metabolites beyond what ``node_is_primary`` alone
    provides (which uses Escher's global primary font size).

    Parameters
    ----------
    escher_map : list
        Escher map structure.
    label_size : float
        Pixel size to write.
    bigg_ids : iterable[str] or None
        Metabolite ``bigg_id`` values to target.
    node_ids : iterable[str] or None
        Node keys in ``escherMap[1]["nodes"]`` to target.
    """
    bigg_set = set(bigg_ids or ())
    node_set = set(node_ids or ())
    for nid, node in escher_map[1]["nodes"].items():
        if node.get("node_type") != "metabolite":
            continue
        if nid in node_set or node.get("bigg_id") in bigg_set:
            node["label_size"] = label_size


def scale_node_label_offsets(escher_map, scale, bigg_ids=None, node_ids=None):
    """Multiply the vector from node center to label position by ``scale``
    for every matching metabolite node. Used to push labels further away
    when their font size has been enlarged, so the larger text does not
    overlap the node.

    Parameters
    ----------
    escher_map : list
        Escher map structure.
    scale : float
        Multiplier applied to ``(label_x - x, label_y - y)``. A value of
        2.0 doubles the distance; 1.0 leaves offsets unchanged.
    bigg_ids : iterable[str] or None
        Metabolite ``bigg_id`` values to target.
    node_ids : iterable[str] or None
        Node keys in ``escherMap[1]["nodes"]`` to target.
    """
    bigg_set = set(bigg_ids or ())
    node_set = set(node_ids or ())
    for nid, node in escher_map[1]["nodes"].items():
        if node.get("node_type") != "metabolite":
            continue
        if not (nid in node_set or node.get("bigg_id") in bigg_set):
            continue
        if "x" in node and "label_x" in node:
            node["label_x"] = node["x"] + (node["label_x"] - node["x"]) * scale
        if "y" in node and "label_y" in node:
            node["label_y"] = node["y"] + (node["label_y"] - node["y"]) * scale


def set_primary_nodes(escher_map, bigg_ids=None, node_ids=None, demote_others=False):
    """Mark specific metabolite nodes as primary so Escher renders them
    larger (a common form of highlighting).

    Parameters
    ----------
    escher_map : list
        Escher map structure.
    bigg_ids : iterable[str] or None
        Metabolite ``bigg_id`` values to promote. Every matching node is
        promoted (useful when the same compound appears in multiple
        places).
    node_ids : iterable[str] or None
        Specific node keys in ``escherMap[1]["nodes"]`` to promote.
    demote_others : bool
        When True, every metabolite not listed is set to
        ``node_is_primary=False``, so the listed nodes are the *only*
        primary nodes. This produces a sharper visual highlight.
    """
    bigg_set = set(bigg_ids or ())
    node_set = set(node_ids or ())
    for nid, node in escher_map[1]["nodes"].items():
        if node.get("node_type") != "metabolite":
            continue
        is_target = nid in node_set or node.get("bigg_id") in bigg_set
        if is_target:
            node["node_is_primary"] = True
        elif demote_others:
            node["node_is_primary"] = False


def remove_text_labels(escher_map, label_ids=None):
    """Remove ``text_labels`` entries. When ``label_ids`` is None, removes
    all text labels; otherwise removes only the given ids.
    """
    text_labels = escher_map[1].get("text_labels", {})
    if label_ids is None:
        escher_map[1]["text_labels"] = {}
        return
    for lid in list(label_ids):
        text_labels.pop(lid, None)


def apply_layout(
    escher_path,
    output_path=None,
    metabolite_label_dx=15,
    metabolite_label_dy=0,
    node_label_positions=None,
    reaction_label_positions=None,
    reaction_label_positions_by_id=None,
    primary_bigg_ids=None,
    primary_node_ids=None,
    demote_other_primaries=False,
    highlight_label_size=None,
    highlight_offset_scale=None,
    remove_all_text_labels=False,
    remove_text_label_ids=None,
):
    """Read an Escher map, apply layout adjustments, and write the result.

    Parameters
    ----------
    escher_path : str or Path
        Input Escher JSON map.
    output_path : str or Path or None
        Output path; defaults to ``<stem>_layout<suffix>``.
    metabolite_label_dx, metabolite_label_dy : float
        Offsets applied by ``normalize_metabolite_label_offsets``. Pass
        ``None`` for either to skip metabolite label normalization.
    node_label_positions : dict[str, tuple[float, float]] or None
        Per-node-id overrides passed to ``set_node_label_positions``. These
        run after normalization so explicit positions override offsets.
    reaction_label_positions : dict[str, tuple[float, float]] or None
        Passed to ``set_reaction_label_positions`` (keyed by reaction
        ``bigg_id``).
    reaction_label_positions_by_id : dict[str, tuple[float, float]] or None
        Passed to ``set_reaction_label_positions_by_id`` (keyed by the
        reaction's key in ``escherMap[1]["reactions"]``).
    primary_bigg_ids, primary_node_ids : iterable[str] or None
        Passed to ``set_primary_nodes`` to highlight key metabolites by
        promoting them to primary (larger) nodes.
    demote_other_primaries : bool
        When True, every metabolite not listed in ``primary_bigg_ids``
        or ``primary_node_ids`` is demoted so only the targets are
        primary. Makes the highlight visually distinct.
    highlight_label_size : float or None
        If set, writes ``label_size=highlight_label_size`` on the listed
        primary metabolites (via ``set_node_label_sizes``).
    highlight_offset_scale : float or None
        If set, multiplies the label-to-node offset vector of the listed
        primary metabolites by this factor (via
        ``scale_node_label_offsets``) so the enlarged labels do not
        overlap their now-larger node circles.
    remove_all_text_labels : bool
        When True, all ``text_labels`` are cleared.
    remove_text_label_ids : iterable[str] or None
        Specific text-label ids to remove (ignored when
        ``remove_all_text_labels`` is True).

    Returns
    -------
    pathlib.Path
        Output path written.
    """
    escher_path = Path(escher_path)
    with escher_path.open("r") as f:
        escher_map = load(f)

    if metabolite_label_dx is not None and metabolite_label_dy is not None:
        normalize_metabolite_label_offsets(
            escher_map, metabolite_label_dx, metabolite_label_dy,
        )
    if node_label_positions:
        set_node_label_positions(escher_map, node_label_positions)
    if reaction_label_positions:
        set_reaction_label_positions(escher_map, reaction_label_positions)
    if reaction_label_positions_by_id:
        set_reaction_label_positions_by_id(escher_map, reaction_label_positions_by_id)
    if primary_bigg_ids or primary_node_ids or demote_other_primaries:
        set_primary_nodes(
            escher_map,
            bigg_ids=primary_bigg_ids,
            node_ids=primary_node_ids,
            demote_others=demote_other_primaries,
        )
    if highlight_label_size is not None and (primary_bigg_ids or primary_node_ids):
        set_node_label_sizes(
            escher_map, highlight_label_size,
            bigg_ids=primary_bigg_ids, node_ids=primary_node_ids,
        )
    if highlight_offset_scale is not None and (primary_bigg_ids or primary_node_ids):
        scale_node_label_offsets(
            escher_map, highlight_offset_scale,
            bigg_ids=primary_bigg_ids, node_ids=primary_node_ids,
        )
    if remove_all_text_labels:
        remove_text_labels(escher_map, None)
    elif remove_text_label_ids:
        remove_text_labels(escher_map, remove_text_label_ids)

    if output_path is None:
        output_path = escher_path.with_name(
            escher_path.stem + "_layout" + escher_path.suffix,
        )
    output_path = Path(output_path)
    with output_path.open("w") as f:
        dump(escher_map, f, indent=3)
    log.info("wrote laid-out Escher JSON to %s", output_path)
    return output_path
