"""Member colours: one distinct colour per community member.

A community map colours everything that belongs to a member — its box and the
edges it consumes and excretes along — in that member's colour, so which
member an edge belongs to can be read off it without tracing it back to the
member column. The 3H11/R12 SynCom figure is the reference: 3H11 blue, R12
orange, metabolite nodes and their labels left neutral so no compound reads as
belonging to one member.

The first eight colours are a validated categorical palette, in its documented
order — lightness held in one band (OKLCH L 0.43-0.77), chroma high enough not
to read as grey, and neighbouring slots kept apart under simulated protanopia
and deuteranopia as well as normal vision. The SynCom's blue and orange are
its first two slots. That validation is of neighbouring slots, though, and on
a map any two members can end up side by side: up to three members every
pair clears the thresholds, but from the fourth some pairs fall short (orange
against yellow, and against red, under normal vision), so the member label
has to carry identity there as well.

Colours follow the member, not its position in the layout:
:func:`member_colors` assigns them in the order it is given, and
:func:`map_member_colors` sorts the member names first, so the same member
keeps its colour in every condition block of a map.

A community larger than eight still gets one colour per member. Past the
eighth, colours are generated rather than drawn from the palette: each is the
colour, within the same lightness band and chroma floor, furthest from every
colour already assigned under normal and simulated colour-blind vision. They
are as far apart as the band allows, but no set of forty hues can all be told
apart, so on a large community the colour is a grouping cue and the member's
label — which every member box carries — is what identifies it.

Colouring by group
------------------
The scheme can colour groups of members instead: pass ``groups`` (``{member:
group}``, or a function of the member) to :func:`member_colors` and every
member of a group shares that group's colour, the groups taking the palette
slots largest first — the same rules, one level up. That is how a map is
coloured by phylogeny: :func:`taxon_groups` reads the group off a taxonomy at
any rank, so members of one phylum, say, come out in one colour, and a legend
(:func:`group_legend`) says which. A member no group claims is drawn grey
(:data:`UNGROUPED`), so a gap in the taxonomy cannot pass for a taxon of its
own. ``group_colors`` fixes the colour of any group outright; build it once
from every group with :func:`group_palette` to colour a series of maps with
different members alike.
"""
import logging
import math
from collections import Counter
from functools import lru_cache

log = logging.getLogger(__name__)

#: The validated categorical slots, in order (light surface).
MEMBER_PALETTE = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                  "#e87ba4", "#008300", "#4a3aa7", "#e34948")

#: Ink for text and neutral marks.
INK = "#222222"

#: The colour of a member no group claims, when members are coloured by group.
UNGROUPED = "#8a8985"
#: The legend's name for those members.
UNGROUPED_LABEL = "Unassigned"

# Rank prefixes of lineage strings ("d__Bacteria;p__Bacillota;..."), and the
# names one rank goes by in different taxonomies.
_RANK_PREFIXES = {"d": "domain", "k": "domain", "p": "phylum", "c": "class",
                  "o": "order", "f": "family", "g": "genus", "s": "species"}
_RANK_ALIASES = {"kingdom": "domain", "superkingdom": "domain"}

# Generated colours past the palette stay inside the palette's own lightness
# band and above its chroma floor, with a little margin on both.
_BAND = (0.47, 0.73)
_LIGHTNESS_STEPS = 4
_HUE_STEP = 5.0
_CHROMA_TARGET, _CHROMA_FLOOR = 0.16, 0.11
# the separations the palette is validated to (OKLab distance x100)
_NORMAL_FLOOR, _CVD_TARGET = 15.0, 8.0

# Machado, Oliveira & Fernandes (2009) dichromacy at severity 1.0, applied in
# linear RGB — the simulation the palette itself was validated against.
_CVD = {
    "protan": ((0.152286, 1.052583, -0.204868),
               (0.114503, 0.786281, 0.099216),
               (-0.003882, -0.048116, 1.051998)),
    "deutan": ((0.367322, 0.860646, -0.227968),
               (0.280085, 0.672501, 0.047413),
               (-0.011820, 0.042940, 0.968881)),
}


def _hex_to_linear(hex_color):
    h = hex_color.lstrip("#")
    return tuple(((c / 255 / 12.92) if c / 255 <= 0.04045
                  else ((c / 255 + 0.055) / 1.055) ** 2.4)
                 for c in (int(h[i:i + 2], 16) for i in (0, 2, 4)))


def _linear_to_hex(rgb):
    def encode(c):
        c = min(1.0, max(0.0, c))
        c = 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055
        return round(c * 255)
    return "#{:02x}{:02x}{:02x}".format(*(encode(c) for c in rgb))


def _cbrt(x):
    return math.copysign(abs(x) ** (1 / 3), x)


def _oklab(rgb):
    r, g, b = rgb
    l = _cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b)
    m = _cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b)
    s = _cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b)
    return (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s)


def _oklch_to_linear(lightness, chroma, hue_deg):
    a = chroma * math.cos(math.radians(hue_deg))
    b = chroma * math.sin(math.radians(hue_deg))
    l = (lightness + 0.3963377774 * a + 0.2158037573 * b) ** 3
    m = (lightness - 0.1055613458 * a - 0.0638541728 * b) ** 3
    s = (lightness - 0.0894841775 * a - 1.2914855480 * b) ** 3
    return (4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
            -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
            -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s)


def _views(hex_color):
    """The colour in OKLab as seen with normal vision, protanopia and
    deuteranopia."""
    rgb = _hex_to_linear(hex_color)
    views = [_oklab(rgb)]
    for matrix in _CVD.values():
        views.append(_oklab(tuple(
            min(1.0, max(0.0, sum(row[i] * rgb[i] for i in range(3))))
            for row in matrix)))
    return views


def _separation(views_a, views_b):
    """How far apart two colours are, as a fraction of what the palette
    requires: OKLab distance (x100) of 15 under normal vision and 8 under
    protanopia or deuteranopia, whichever is closer to failing."""
    normal = 100 * math.dist(views_a[0], views_b[0]) / _NORMAL_FLOOR
    cvd = 100 * min(math.dist(a, b) for a, b in zip(views_a[1:], views_b[1:]))
    return min(normal, cvd / _CVD_TARGET)


@lru_cache(maxsize=1)
def _candidates():
    """Every in-gamut colour in the band on a lightness x hue grid, at the
    target chroma or the most the gamut allows there."""
    found = []
    for step in range(_LIGHTNESS_STEPS):
        lightness = _BAND[0] + step * (_BAND[1] - _BAND[0]) / (_LIGHTNESS_STEPS - 1)
        for index in range(int(360 / _HUE_STEP)):
            hue = index * _HUE_STEP
            lo, hi = 0.0, _CHROMA_TARGET
            if not all(0 <= c <= 1 for c in _oklch_to_linear(lightness, hi, hue)):
                for _ in range(30):           # largest chroma still in gamut
                    mid = (lo + hi) / 2
                    in_gamut = all(0 <= c <= 1
                                   for c in _oklch_to_linear(lightness, mid, hue))
                    lo, hi = (mid, hi) if in_gamut else (lo, mid)
                hi = lo
            if hi < _CHROMA_FLOOR:
                continue
            hex_color = _linear_to_hex(_oklch_to_linear(lightness, hi, hue))
            found.append((hex_color, _views(hex_color)))
    return tuple(found)


@lru_cache(maxsize=None)
def member_palette(n):
    """``n`` member colours: the validated slots first, then generated ones.

    Deterministic — the same ``n`` always gives the same list, and every
    shorter list is a prefix of a longer one, so adding members never
    repaints the ones already coloured.
    """
    if n <= len(MEMBER_PALETTE):
        return MEMBER_PALETTE[:n]
    colors = list(member_palette(n - 1))
    chosen = [_views(c) for c in colors]
    best = max(_candidates(),
               key=lambda cand: min(_separation(cand[1], v) for v in chosen))
    return tuple(colors + [best[0]])


def _group_lookup(groups):
    """``groups`` — a mapping, anything else with ``.get`` (a pandas Series)
    or a function — as a function of the member giving its group, or None for
    a member without one (missing, empty or NaN)."""
    get = groups.get if hasattr(groups, "get") else groups

    def group_of(member):
        group = get(member)
        if isinstance(group, str):
            group = group.strip()
        if group is None or group == "" or group != group:
            return None
        return group
    return group_of


def group_palette(groups, group_colors=None):
    """``{group: "#rrggbb"}``, one colour per distinct group in ``groups``.

    ``groups`` gives each member's group — ``phylum_of.values()``, say — so
    a group appears once per member. The groups take the member palette's
    slots largest first, ties in sorted order: the groups most of the map is
    drawn in get the validated colours, and any past the eighth, the
    generated ones, go to the smallest. ``group_colors`` (``{group:
    colour}``) fixes the colour of the groups it lists, and the rest take the
    slots it leaves free. None is not a group and is skipped.

    To colour a series of maps alike when not every group is in every map,
    build this once from every group and pass it to each as ``group_colors``.
    """
    sizes = Counter(group for group in groups if group is not None)
    distinct = sorted(sizes, key=lambda group: (-sizes[group], str(group)))
    fixed = {group: "#" + color.lstrip("#")
             for group, color in (group_colors or {}).items()}
    free = [group for group in distinct if group not in fixed]
    taken = {color.lower() for color in fixed.values()}
    if len(free) > len(MEMBER_PALETTE) - len(taken & set(MEMBER_PALETTE)):
        log.warning(
            "%d groups but %d validated member colours: the rest are "
            "generated to be as distinct as possible, but not all %d can be "
            "told apart by colour alone — group at a broader rank, or rely on "
            "the legend and member labels",
            len(distinct), len(MEMBER_PALETTE), len(distinct))
    slots = [color for color in member_palette(len(free) + len(taken))
             if color not in taken]
    colors = {group: fixed[group] for group in distinct if group in fixed}
    colors.update(zip(free, slots))
    return colors


def member_colors(names, groups=None, group_colors=None):
    """``{member: "#rrggbb"}``, one distinct colour per member, in the order
    ``names`` gives them (duplicates keep their first colour).

    With ``groups`` — ``{member: group}`` or a function of the member, such
    as :func:`taxon_groups` gives — members are coloured by group instead:
    every member of a group gets its colour from :func:`group_palette`
    (``group_colors`` fixing any of them), and a member without a group is
    drawn in :data:`UNGROUPED` grey, with a warning naming it.

    To colour a series of maps consistently — one per condition, say, where
    a member may be missing from some — build this once from every member and
    hand the same mapping to each.
    """
    unique = list(dict.fromkeys(names))
    if groups is None:
        if group_colors:
            log.warning("group_colors has no effect without groups")
        if len(unique) > len(MEMBER_PALETTE):
            log.warning(
                "%d members but %d validated member colours: the rest are "
                "generated to be as distinct as possible, but not all %d can "
                "be told apart by colour alone — rely on the member labels",
                len(unique), len(MEMBER_PALETTE), len(unique))
        return dict(zip(unique, member_palette(len(unique))))

    group_of = _group_lookup(groups)
    member_group = {name: group_of(name) for name in unique}
    ungrouped = [name for name, group in member_group.items() if group is None]
    if ungrouped:
        log.warning("%d of %d members have no group and are drawn grey: %s",
                    len(ungrouped), len(unique),
                    ", ".join(map(str, ungrouped[:10]))
                    + (", ..." if len(ungrouped) > 10 else ""))
    colors = group_palette(member_group.values(), group_colors)
    return {name: UNGROUPED if group is None else colors[group]
            for name, group in member_group.items()}


def map_member_colors(escher_map, colors=None, groups=None, group_colors=None):
    """``{reaction key: "#rrggbb"}`` for every member reaction in a map.

    A member reaction's ``name`` is the member, whichever condition block it
    sits in (its ``bigg_id`` carries the condition too), so every block
    colours a member alike. ``colors`` is a ``{member: colour}`` mapping from
    :func:`member_colors`; without it the map's own members are coloured in
    sorted order — by ``groups`` and ``group_colors`` when given (see
    :func:`member_colors`) — so the assignment depends only on who is in the
    community.
    """
    reactions = escher_map[1]["reactions"]
    names = {key: rxn.get("name") or rxn.get("bigg_id", "")
             for key, rxn in reactions.items()}
    if colors is None:
        colors = member_colors(sorted(set(names.values())), groups, group_colors)
    return {key: colors[name] for key, name in names.items() if name in colors}


def group_legend(colors, groups):
    """``[(group, "#rrggbb")]``: the legend for members coloured by group.

    ``colors`` is ``{member: colour}`` for the members drawn, ``groups`` what
    they were grouped by. One entry per group among them, largest first as
    :func:`group_palette` assigns the colours, and last, if any member has no
    group, :data:`UNGROUPED_LABEL`. A group whose members were drawn in different
    colours — from a ``{member: colour}`` mapping not built from these groups
    — is listed in the first, with a warning.
    """
    group_of = _group_lookup(groups)
    found = {}
    for member, color in colors.items():
        found.setdefault(group_of(member), []).append(color)
    entries = []
    order = sorted((group for group in found if group is not None),
                   key=lambda group: (-len(found[group]), str(group)))
    for group in order + ([None] if None in found else []):
        shades = list(dict.fromkeys(found[group]))
        label = UNGROUPED_LABEL if group is None else str(group)
        if len(shades) > 1:
            log.warning("the members of %s are drawn in %d colours; the "
                        "legend shows %s", label, len(shades), shades[0])
        entries.append((label, shades[0]))
    return entries


def _rank_key(rank):
    rank = str(rank).strip().lower()
    return _RANK_ALIASES.get(rank, rank)


def taxon_groups(taxonomy, rank):
    """``{member: taxon}``: each member's taxon at ``rank``, as ``groups`` for
    :func:`member_colors` — colouring a map by phylogeny.

    ``taxonomy`` maps each member to its lineage, either as ``{rank: name}``
    — rank keys in any case, ``kingdom`` standing for ``domain`` — or as a
    lineage string of prefixed ranks (``"d__Bacteria;p__Bacillota;..."``, as
    GTDB writes them). Names are stripped of surrounding whitespace; a member
    whose lineage leaves that rank empty (``"g__"``, ``""``, NaN) or out gets
    None, and is drawn as ungrouped.
    """
    wanted = _rank_key(rank)
    groups = {}
    for member, lineage in taxonomy.items():
        if isinstance(lineage, str):
            ranks = {}
            for part in lineage.split(";"):
                prefix, sep, name = part.strip().partition("__")
                if sep and prefix.lower() in _RANK_PREFIXES:
                    ranks[_RANK_PREFIXES[prefix.lower()]] = name
            lineage = ranks
        elif not hasattr(lineage, "items"):
            raise TypeError(f"{member}: a lineage is a {{rank: name}} mapping "
                            f"or a 'd__...;p__...' string, not "
                            f"{type(lineage).__name__}")
        taxa = {_rank_key(key): name for key, name in lineage.items()}
        groups[member] = _group_lookup(taxa)(wanted)
    if taxonomy and all(group is None for group in groups.values()):
        log.warning("no member has a taxon at rank %r", rank)
    return groups


def text_color_on(fill_hex):
    """White or ink, whichever contrasts more with ``fill_hex`` — for a label
    drawn on a member-coloured box."""
    def luminance(hex_color):
        r, g, b = _hex_to_linear(hex_color)
        return 0.2126 * r + 0.7152 * g + 0.0722 * b
    fill = luminance(fill_hex)
    on_white = (1.05) / (fill + 0.05)
    on_ink = (fill + 0.05) / (luminance(INK) + 0.05)
    return "#ffffff" if on_white >= on_ink else INK
