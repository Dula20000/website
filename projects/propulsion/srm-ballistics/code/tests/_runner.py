"""Minimal runner so the tests work with plain `python3 tests/test_x.py` when pytest is absent."""
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run(module_globals):
    names = [n for n in module_globals if n.startswith("test_") and callable(module_globals[n])]
    fails = 0
    for n in names:
        t = time.time()
        try:
            module_globals[n]()
            print(f"PASS {n} ({time.time() - t:.2f}s)")
        except Exception:
            fails += 1
            print(f"FAIL {n}")
            traceback.print_exc()
    print(f"{len(names) - fails}/{len(names)} passed")
    sys.exit(1 if fails else 0)
