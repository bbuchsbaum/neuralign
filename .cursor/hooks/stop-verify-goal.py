#!/usr/bin/env python3
"""stop hook: continue until scripts/verify-goal succeeds.

Completion is a predicate, not the model's claim that it is done.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VERIFY = ROOT / "scripts" / "verify-goal"


def emit(payload: dict) -> None:
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def main() -> int:
    # Drain stdin so the hook host does not see a broken pipe.
    try:
        sys.stdin.read()
    except OSError:
        pass

    if not VERIFY.is_file():
        emit(
            {
                "followup_message": (
                    "scripts/verify-goal is missing. Recreate it from GOAL.md / PLAN.md "
                    "and keep executing Stage 1 only. Do not stop."
                )
            }
        )
        return 0

    try:
        proc = subprocess.run(
            [sys.executable, str(VERIFY)],
            cwd=ROOT,
            text=True,
            capture_output=True,
            timeout=960,
            check=False,
        )
    except subprocess.TimeoutExpired:
        emit(
            {
                "followup_message": (
                    "scripts/verify-goal timed out. Diagnose the hang, then continue "
                    "Stage 1. Do not weaken tests. Do not start Stage 2."
                )
            }
        )
        return 0

    if proc.returncode == 0:
        emit({})
        return 0

    detail = (proc.stderr or proc.stdout or "verify-goal failed").strip()
    if len(detail) > 4000:
        detail = detail[-4000:]

    emit(
        {
            "followup_message": (
                "Stage 1 is not complete. scripts/verify-goal failed:\n\n"
                f"{detail}\n\n"
                "Continue executing PLAN.md in order. Update checkboxes only after "
                "the relevant tests pass. Do not weaken tests or acceptance criteria. "
                "Do not start Stage 2. Run scripts/verify-goal before stopping again."
            )
        }
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
