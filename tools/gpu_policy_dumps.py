#!/usr/bin/env python3
"""Add binary dumps of complete policy/value arrays for direct comparison."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


INTEGER_NEEDLE = "\t\ttopct,ntop,topind(endognum)"
INTEGER_REACH_NEEDLE = (
    "\t\ttopct,ntop,topind(endognum),reach_pos,reach_global_count,&\n"
    "\t\treach_local_min,reach_local_max,reach_global_min,reach_global_max,reach_updates"
)
EXACT_HEADER_NEEDLE = (
    "subroutine politeration()\n"
    "\tuse params\n"
    "\tuse base_lib\n"
    "\tuse omp_lib\n"
    "\timplicit none\n"
)
POLITERATION_HEADER_RE = re.compile(
    r"(subroutine politeration\(\)\n(?:(?!end subroutine politeration).)*?"
    r"^[ \t]*implicit none[ \t]*\n)",
    re.MULTILINE | re.DOTALL,
)

DUMP_BLOCK = """open(newunit=dump_unit,file='RCC_kbpolind.bin',access='stream',&
\tform='unformatted',status='replace',action='write')
write(dump_unit) kbpolind
close(dump_unit)
open(newunit=dump_unit,file='RCC_dfpolind.bin',access='stream',&
\tform='unformatted',status='replace',action='write')
write(dump_unit) dfpolind
close(dump_unit)
open(newunit=dump_unit,file='RCC_Vtheta.bin',access='stream',&
\tform='unformatted',status='replace',action='write')
write(dump_unit) Vtheta
close(dump_unit)
open(newunit=dump_unit,file='RCC_Vthetaold.bin',access='stream',&
\tform='unformatted',status='replace',action='write')
write(dump_unit) Vthetaold
close(dump_unit)
open(newunit=dump_unit,file='RCC_qmat.bin',access='stream',&
\tform='unformatted',status='replace',action='write')
write(dump_unit) qmat
close(dump_unit)

"""


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise ValueError(f"expected one {label}; found {count}")
    return text.replace(old, new, 1)


def transform(text: str) -> str:
    if text.count(INTEGER_REACH_NEEDLE) == 1:
        text = replace_once(
            text,
            INTEGER_REACH_NEEDLE,
            INTEGER_REACH_NEEDLE + ",dump_unit",
            "reach integer declaration",
        )
    elif text.count(INTEGER_NEEDLE) == 1:
        text = replace_once(
            text,
            INTEGER_NEEDLE,
            INTEGER_NEEDLE + ",dump_unit",
            "base integer declaration",
        )
    else:
        match = POLITERATION_HEADER_RE.search(text)
        if match is None:
            text = replace_once(
                text,
                EXACT_HEADER_NEEDLE,
                EXACT_HEADER_NEEDLE + "\n\tinteger :: dump_unit\n",
                "exact politeration header",
            )
        else:
            text = text[: match.end()] + "\n\tinteger :: dump_unit\n" + text[match.end() :]

    marker = "deallocate(EVtheta,EVthetaT,qmatT,currentdiv,ynetcache,accache,currentbase,&"
    position = text.find(marker)
    if position < 0:
        marker = "deallocate(EVthetaT,qmatT,currentdiv)"
        position = text.find(marker)
    if position < 0:
        start = text.find("subroutine politeration()")
        end = text.find("end subroutine politeration", start)
        if start < 0 or end < 0:
            raise ValueError("policy subroutine not found")
        match = re.search(r"^[ \t]*deallocate\(", text[start:end], re.MULTILINE)
        if match is None:
            raise ValueError("policy deallocation marker not found")
        position = start + match.start()
    return text[:position] + DUMP_BLOCK + text[position:]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    args.source.write_text(transform(args.source.read_text()))


if __name__ == "__main__":
    main()
