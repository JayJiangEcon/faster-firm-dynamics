#!/usr/bin/env python3
"""Dump a compact, complete RE or DE CUDA fixture after one ordinary CPU solve."""

from __future__ import annotations

import argparse
from pathlib import Path

from gpu_policy_dumps import transform as add_policy_dumps


PROFIT_END = """end do !endogct
end do !exogct

end subroutine politeration"""

GPU_DUMPS = """end do !endogct
end do !exogct

open(newunit=dump_unit,file='RCC_pr_mat_exog_theta.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) pr_mat_exog_theta
close(dump_unit)
open(newunit=dump_unit,file='RCC_pr_mat_exog.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) pr_mat_exog
close(dump_unit)
open(newunit=dump_unit,file='RCC_pr_mat_A.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) pr_mat_A
close(dump_unit)
open(newunit=dump_unit,file='RCC_pr_mat_z.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) pr_mat_z
close(dump_unit)
open(newunit=dump_unit,file='RCC_exogind.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) exogind
close(dump_unit)
open(newunit=dump_unit,file='RCC_endogind.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) endogind
close(dump_unit)
open(newunit=dump_unit,file='RCC_exogval.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) exogval
close(dump_unit)
open(newunit=dump_unit,file='RCC_endogval.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) endogval
close(dump_unit)
open(newunit=dump_unit,file='RCC_A0.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) A0
close(dump_unit)
open(newunit=dump_unit,file='RCC_z0.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) z0
close(dump_unit)
open(newunit=dump_unit,file='RCC_k0.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) k0
close(dump_unit)
open(newunit=dump_unit,file='RCC_b0.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) b0
close(dump_unit)
open(newunit=dump_unit,file='RCC_dist.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) dist
close(dump_unit)
open(newunit=dump_unit,file='RCC_pimat.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) pimat
close(dump_unit)
open(newunit=dump_unit,file='RCC_divposmat.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) divposmat
close(dump_unit)

end subroutine politeration"""

QMAT_DUMP_END = """write(dump_unit) qmat
close(dump_unit)

"""

CURRENT_DIVIDEND_DUMP = """write(dump_unit) qmat
close(dump_unit)
open(newunit=dump_unit,file='RCC_currentdiv.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) currentdiv
close(dump_unit)
open(newunit=dump_unit,file='RCC_ynetcache.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) ynetcache
close(dump_unit)
open(newunit=dump_unit,file='RCC_accache.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) accache
close(dump_unit)
open(newunit=dump_unit,file='RCC_actionkbase.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(dump_unit) actionkbase
close(dump_unit)

"""

SIMULATION_END = """write(*,*) " "

end subroutine simulate"""

SIMULATION_DUMPS = """write(*,*) " "

! Fixed-draw simulation inputs and independently produced CPU targets for the
! CUDA full-evaluator gate.  These dumps occur only after the released
! simulation and moment calculation have completed.
open(unit=47,file='RCC_Asimpos.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) Asimpos
close(47)
open(unit=47,file='RCC_zfirmpos.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) zfirmpos
close(47)
open(unit=47,file='RCC_pifirmshocks.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) pifirmshocks
close(47)
open(unit=47,file='RCC_sim_endogfirmpos.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) endogfirmpos
close(47)
open(unit=47,file='RCC_sim_kfirm.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) kfirm
close(47)
open(unit=47,file='RCC_sim_bfirm.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) bfirm
close(47)
open(unit=47,file='RCC_sim_ifirm.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) ifirm
close(47)
open(unit=47,file='RCC_sim_pifirm.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) pifirm
close(47)
open(unit=47,file='RCC_sim_pifcstfirm.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) pifcstfirm
close(47)
open(unit=47,file='RCC_sim_sprfirm.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) sprfirm
close(47)
open(unit=47,file='RCC_sim_Ysim.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) Ysim
close(47)
open(unit=47,file='RCC_sim_dfsim.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) dfsim
close(47)
open(unit=47,file='RCC_sim_spreadsim.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) spreadsim
close(47)
open(unit=47,file='RCC_GMMmodelmom.bin',access='stream',&
  form='unformatted',status='replace',action='write')
write(47) GMMmodelmom
close(47)

end subroutine simulate"""


SIMULATION_END += "\n"
SIMULATION_DUMPS += "\n"

def transform(source: str) -> str:
    source = add_policy_dumps(source)
    count = source.count(QMAT_DUMP_END)
    if count != 1:
        raise ValueError(f"expected one qmat dump end; found {count}")
    source = source.replace(QMAT_DUMP_END, CURRENT_DIVIDEND_DUMP, 1)
    count = source.count(PROFIT_END)
    if count != 1:
        raise ValueError(f"expected one final profit block; found {count}")
    source = source.replace(PROFIT_END, GPU_DUMPS, 1)
    count = source.count(SIMULATION_END)
    if count != 1:
        raise ValueError(f"expected one simulation end; found {count}")
    return source.replace(SIMULATION_END, SIMULATION_DUMPS, 1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    arguments = parser.parse_args()
    arguments.source.write_text(
        transform(arguments.source.read_text(encoding="utf-8")), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
