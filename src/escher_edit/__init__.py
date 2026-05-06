"""Tools for editing Escher metabolic-map JSON and rendered SVG output.

Public surface:
    clean_json.cleanEscherJSON, CPD_ID_ABBRV
    filter_map.filter_escher_map, AA_NAMES, DEFAULT_SKIP_NAMES
    layout.apply_layout, normalize_metabolite_label_offsets,
        set_reaction_label_positions, remove_text_labels
    model_mapping.build_direction_tracking, categorize_segments_and_nodes
    reverse_reactions.reverse_reaction, reverse_reactions_in_map
    svg_editor.EscherSVG_processing, EscherStyle
"""
from .clean_json import CPD_ID_ABBRV, cleanEscherJSON, build_name_abbrev_table
from .filter_map import AA_NAMES, DEFAULT_SKIP_NAMES, filter_escher_map, build_node_lookups
from .layout import (
    apply_layout,
    normalize_metabolite_label_offsets,
    remove_text_labels,
    scale_node_label_offsets,
    set_node_label_positions,
    set_node_label_sizes,
    set_primary_nodes,
    set_reaction_label_positions,
    set_reaction_label_positions_by_id,
)
from .model_mapping import (
    DEFAULT_DIETS, DEFAULT_DAYS,
    build_direction_tracking, categorize_segments_and_nodes,
)
from .reverse_reactions import reverse_reaction, reverse_reactions_in_map
from .svg_editor import (
    EscherSVG_processing,
    EscherStyle,
    extract_per_node_label_sizes_from_json,
)

__all__ = [
    "CPD_ID_ABBRV", "cleanEscherJSON", "build_name_abbrev_table",
    "AA_NAMES", "DEFAULT_SKIP_NAMES", "filter_escher_map", "build_node_lookups",
    "apply_layout", "normalize_metabolite_label_offsets",
    "remove_text_labels", "scale_node_label_offsets",
    "set_node_label_positions", "set_node_label_sizes", "set_primary_nodes",
    "set_reaction_label_positions", "set_reaction_label_positions_by_id",
    "DEFAULT_DIETS", "DEFAULT_DAYS",
    "build_direction_tracking", "categorize_segments_and_nodes",
    "reverse_reaction", "reverse_reactions_in_map",
    "EscherSVG_processing", "EscherStyle", "extract_per_node_label_sizes_from_json",
]
