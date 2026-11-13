"""
Read parameters from the Fortran sources (controles.f90, main.f90,
system_parameters.f90, global.f90).
"""
import re
from collections import namedtuple

LOGICAL = {"yes": True, "no": False, ".true.": True, ".false.": False}
NUM = r"([-\d.dDeE+]+)"
UNITS = {"kpc": 3.085677581e21, "pc": 3.085677581e18}

MapFreq = namedtuple("MapFreq", "nu label sigma on") # sigma: (x, y) ["]


def number(s):
    return float(s.strip().lower().replace("d", "e"))


def source(path):
    # file text without comments (! outside strings) and with joined & lines
    lines = []
    with open(path) as f:
        for line in f:
            in_str = False
            for i, ch in enumerate(line):
                if ch == "'":
                    in_str = not in_str
                elif ch == "!" and not in_str:
                    line = line[:i]
                    break
            lines.append(line.rstrip())
    return re.sub(r"&\s*\n\s*&?", " ", "\n".join(lines))


def search(paths, pattern):
    # first match of pattern in a list of files (None if not found)
    for p in paths:
        try:
            m = re.search(pattern, source(p), re.I)
        except OSError:
            continue
        if m:
            return m
    return None


def read_controls(path):
    # {name: bool} for every 'logical, parameter', names in lower case
    flags = {}
    for name, val in re.findall(r"logical\s*,\s*parameter\s*::\s*(\w+)\s*=\s*([\w.]+)",
                                source(path), re.I):
        flags[name.lower()] = LOGICAL.get(val.lower(), flags.get(val.lower(), False))
    return flags


def read_parameter(paths, name, default=None):
    m = search(paths, rf"\b{re.escape(name)}\s*=\s*{NUM}")
    return number(m.group(1)) if m else default


def read_string(paths, name, default=None):
    m = search(paths, rf"\b{re.escape(name)}\s*=\s*'([^']*)'")
    return m.group(1) if m else default


def read_distance(paths, default=None):
    # distance = <number>*kpc (or pc) -> cm
    m = search(paths, rf"\bdistance\s*=\s*{NUM}\s*\*\s*(kpc|pc)\b")
    return number(m.group(1)) * UNITS[m.group(2).lower()] if m else default


def read_map_catalogue(path):
    # map_catalogue and map_avg of controles.f90
    item = rf"map_freq_t\(\s*{NUM}\s*,\s*'([^']*)'\s*,\s*\[\s*([^,\]]+),\s*([^\]]+)\]\s*,\s*([\w.]+)\s*\)"

    def freq(m):
        return MapFreq(number(m[0]), m[1].strip(), (number(m[2]), number(m[3])),
                       LOGICAL.get(m[4].lower(), False))

    text = source(path)
    block = re.search(r"map_catalogue\s*\(\s*\*\s*\)\s*=\s*\[(.*?)\]\s*\n", text, re.I | re.S)
    avg = re.search(r"map_avg\s*=\s*" + item, text, re.I)
    return ([freq(m) for m in re.findall(item, block.group(1), re.I)] if block else [],
            freq(avg.groups()) if avg else None)
