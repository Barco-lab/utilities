#!/usr/bin/env python3
"""
旧版 snakemake_bulk_chip / snakemake_bulk_atac 的 bigwig 改为带归一化方式的新命名，不重新计算：
  results/04_bigwig/{sample}.bw         -> results/04_bigwig/{sample}.{tag}.bw
  results/04_bigwig/groups/{group}.bw   -> results/04_bigwig/groups/{group}.{tag}.bw
tag 为旧 bigwig 实际使用的 bigwig.normalize（默认读项目 config.yaml，旧版只能是一个值）。
只改样本表中的样本和组（group 列为空时同流程，去掉样本名末尾的 _repN）；目录中的其它 bw（外部数据、软链接）不动；目标已存在时跳过。
默认只打印，加 --apply 才执行。

用法：python rename_bigwig_norm.py <项目输出目录，含 config.yaml> [--tag RPGC] [--apply]
"""
import argparse
import csv
import re
import sys
from pathlib import Path

import yaml

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("outdir", type=Path)
ap.add_argument("--tag", help="旧 bigwig 的归一化方式；默认取 config.yaml 的 bigwig.normalize")
ap.add_argument("--apply", action="store_true", help="执行改名（默认只打印）")
args = ap.parse_args()

cfg_path = args.outdir / "config.yaml"
cfg = yaml.safe_load(cfg_path.read_text())
tag = args.tag or cfg.get("bigwig", {}).get("normalize", "RPGC")
if not isinstance(tag, str):
    sys.exit(f"config 的 bigwig.normalize 为 {tag}，不是旧版的单个值；用 --tag 指定旧 bigwig 的归一化方式")

with open(cfg_path.parent / cfg["samples"]) as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
samples = [r["sample"] for r in rows]
# group 列为空时按样本名推断（同流程 common.smk 的 REP_SUFFIX）：去掉末尾的 _repN
REP_SUFFIX = re.compile(r"[._-]?rep\d+$", re.I)
groups = sorted({r.get("group") or REP_SUFFIX.sub("", r["sample"]) for r in rows})

bw = args.outdir / "results" / "04_bigwig"
pairs = [(bw / f"{s}.bw", bw / f"{s}.{tag}.bw") for s in samples]
pairs += [(bw / "groups" / f"{g}.bw", bw / "groups" / f"{g}.{tag}.bw") for g in groups]

n = 0
for old, new in pairs:
    if not old.exists() or old.is_symlink():
        continue
    if new.exists():
        print(f"跳过（目标已存在）：{new}")
        continue
    print(f"{old.relative_to(args.outdir)} -> {new.name}")
    if args.apply:
        old.rename(new)  # 同一目录内改名，保留修改时间，snakemake 不会因此重跑下游
    n += 1
print(f"{'已改名' if args.apply else '将改名（加 --apply 执行）'} {n} 个文件，tag = {tag}")
