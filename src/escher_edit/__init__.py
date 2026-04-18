"""Tools for editing Escher metabolic-map JSON and rendered SVG output.

Public surface:
    clean_json.cleanEscherJSON, CPD_ID_ABBRV
    filter_map.filter_escher_map, AA_NAMES, DEFAULT_SKIP_NAMES
    model_mapping.build_direction_tracking, categorize_segments_and_nodes
    svg_editor.EscherSVG_processing, EscherStyle
"""
from .clean_json import CPD_ID_ABBRV, cleanEscherJSON, build_name_abbrev_table
from .filter_map import AA_NAMES, DEFAULT_SKIP_NAMES, filter_escher_map, build_node_lookups
from .model_mapping import (
    DEFAULT_DIETS, DEFAULT_DAYS,
    build_direction_tracking, categorize_segments_and_nodes,
)
from .svg_editor import EscherSVG_processing, EscherStyle

__all__ = [
    "CPD_ID_ABBRV", "cleanEscherJSON", "build_name_abbrev_table",
    "AA_NAMES", "DEFAULT_SKIP_NAMES", "filter_escher_map", "build_node_lookups",
    "DEFAULT_DIETS", "DEFAULT_DAYS",
    "build_direction_tracking", "categorize_segments_and_nodes",
    "EscherSVG_processing", "EscherStyle",
]
