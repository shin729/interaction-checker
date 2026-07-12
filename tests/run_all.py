# -*- coding: utf-8 -*-
"""
オフライン単体テストを一括実行する。

  python tests/run_all.py    # test_pure と test_parsers を続けて実行し、
                             # どちらかに失敗があれば exit 1

ネットワーク・キャッシュに依存する統合テストは validate.py 側（別実行）。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import test_matrix
import test_parsers
import test_pure


def main():
    rc = 0
    for mod in (test_pure, test_parsers, test_matrix):
        rc |= mod.main()
    print("\n" + "=" * 70)
    print("[全体] " + ("すべて通過" if rc == 0 else "失敗あり（上記参照）"))
    return rc


if __name__ == "__main__":
    sys.exit(main())
