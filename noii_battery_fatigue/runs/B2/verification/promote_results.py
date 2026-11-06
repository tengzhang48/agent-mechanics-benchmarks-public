"""Historical ordinary-strain promotion helper, retained for provenance.

The enhanced formulation uses per-run source snapshots instead. This helper
does no promotion when that formulation is selected.
"""
from pathlib import Path
import json
import hashlib

ROOT=Path(__file__).resolve().parents[1]
FINAL=ROOT/"element/uel_contract.json"


def normalized(path):
    p=Path(path);fixed=p.suffix==".for"
    return "\n".join(line.rstrip() for line in p.read_text().splitlines()
        if line.strip() and not (line.lstrip().startswith("!") or (fixed and line[0] in "Cc*")))


def main():
    c=json.loads(FINAL.read_text())
    if c.get("element_local_variables"):
        print("Historical helper: enhanced-source runs already preserve snapshots; no promotion performed.")
        return
    prototypes=[
        ("D0_lag_lumped16_long","D0_structured_16_s0_spc40_cy20"),
        ("D0_lag_lumped32_long","D0_structured_32_s0_spc80_cy20"),
        ("D0_lag_lumped64_long","D0_structured_64_s0_spc40_cy20"),
        ("severe32_lag_lumped","V_severe_structured_32_s0_spc80_cy5")]
    for source,destination in prototypes:
        src=ROOT/"verification/exploration"/source
        if not (src/"run.json").exists():
            continue
        result=json.loads((src/"run.json").read_text())
        if result["status"]!="ok":
            continue
        old=json.loads((src/"contract.json").read_text())
        assert len(old["source"])==len(c["source"])
        for before,after in zip(old["source"],c["source"]):
            assert normalized(ROOT/before)==normalized(ROOT/after),(before,after)
        result["original_source_sha256"]=result.get("source_sha256",{})
        result["source_sha256"]={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in c["source"]}
        result["source_note"]="Prototype implementation is identical to final source after excluding comments and blank lines. Only wrapper comments differ. Original compile fingerprints are retained."
        current={**c,"props":old["props"]}
        (src/"contract.json").write_text(json.dumps(current,indent=2)+"\n")
        (src/"run.json").write_text(json.dumps(result,indent=2)+"\n")
        dst=ROOT/"verification/runs"/destination
        assert not dst.exists(),dst
        src.rename(dst)
        print("Promoted",destination)


if __name__=="__main__":
    main()
