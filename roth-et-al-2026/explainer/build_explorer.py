"""Build the standalone Sliding RP explorer page.

Inlines srp_core.js (the tested numerical core) into explorer_template.html so
the result is one self-contained file with no network dependencies: it opens
from disk in any browser. Writes it next to this script and to the paper's
explainer folder.

Run:  python build_explorer.py     (after: python export_reference.py && node test_core.js)
"""
from pathlib import Path

HERE = Path(__file__).parent
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/explainer")
MARK = "/*__SRP_CORE__*/"


def main():
    tpl = (HERE / "explorer_template.html").read_text(encoding="utf-8")
    core = (HERE / "srp_core.js").read_text(encoding="utf-8")
    assert tpl.count(MARK) == 1, "template must contain the core marker exactly once"
    assert "</script" not in core.lower(), "core must not close the script element"
    page = tpl.replace(MARK, core)
    for d in (HERE, OUTDIR):
        d.mkdir(parents=True, exist_ok=True)
        (d / "sliding_rp_explorer.html").write_text(page, encoding="utf-8")
        print("wrote", d / "sliding_rp_explorer.html", f"({len(page) / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
