"""
生成 bulk RNA / ChIP / ATAC 流程 README 中的流程图（SVG + 2 倍 PNG）。

用法：
    python make_workflow_flowcharts.py            # 写入 ../snakemake_bulk_*/docs/images/
    python make_workflow_flowcharts.py OUT_DIR    # 全部写入 OUT_DIR

流程有改动时修改下面 rna() / peaks() 中的节点和连线后重新运行。PNG 需要 rsvg-convert。
"""
import os
import shutil
import subprocess
import sys
from html import escape

FONT = ("Helvetica, Arial, 'PingFang SC', 'Microsoft YaHei', 'Noto Sans CJK SC', "
        "'WenQuanYi Micro Hei', sans-serif")
MONO = "Menlo, Consolas, 'DejaVu Sans Mono', monospace"

STYLE = {  # fill, stroke, title color, dashed
    "input": ("#F1F3F5", "#495057", "#212529", False),
    "main": ("#E7F0FA", "#2B6CB0", "#1A4E8A", False),
    "qc": ("#E8F6EE", "#2F855A", "#22643F", False),
    "opt": ("#FFFFFF", "#B7791F", "#8A5A12", True),
    "report": ("#F1EAFB", "#6B46C1", "#4C2F94", False),
}
OUT_COLOR = "#C0501F"

GUTTER = 128      # stage label column
COLW = 212        # column pitch
BOXW = 192
PAD_X = 24
TOP = 92          # title + legend
ROW_PAD = 26
LINE_T, LINE_S, LINE_O = 19, 16, 16


class Node:
    def __init__(self, id, row, col, title, sub=(), out=None, kind="main", span=1.0):
        self.id, self.row, self.col, self.title = id, row, col, title
        self.sub, self.out, self.kind, self.span = list(sub), out, kind, span

    @property
    def w(self):
        return BOXW + (self.span - 1) * COLW

    @property
    def h(self):
        return 14 + LINE_T + LINE_S * len(self.sub) + (LINE_O + 2 if self.out else 0) + 8


def render(title, subtitle, stages, nodes, edges, ncols, path):
    nodes = {n.id: n for n in nodes}
    width = GUTTER + PAD_X + ncols * COLW + PAD_X - (COLW - BOXW)
    # row geometry: height = tallest node in that row
    row_h, row_y = {}, {}
    y = TOP
    for r in range(len(stages)):
        hs = [n.h for n in nodes.values() if n.row == r] or [40]
        row_h[r] = max(hs) + 2 * ROW_PAD
        row_y[r] = y
        y += row_h[r]
    height = y + 20

    for n in nodes.values():
        n.cx = GUTTER + PAD_X + n.col * COLW + BOXW / 2
        n.cy = row_y[n.row] + row_h[n.row] / 2
        n.x, n.y = n.cx - n.w / 2, n.cy - n.h / 2

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
         f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="{escape(FONT)}">']
    o.append(f'<rect width="100%" height="100%" fill="#FFFFFF"/>')
    o.append(f'<text x="24" y="38" font-size="22" font-weight="700" fill="#1F2933">{escape(title)}</text>')
    o.append(f'<text x="24" y="62" font-size="13" fill="#52606D">{escape(subtitle)}</text>')

    # legend (top right)
    items = [("main", "主流程"), ("opt", "可选 / 条件步骤"), ("qc", "QC"), ("report", "汇总")]
    ws = [28 + text_w(lab, 12.5) + 22 for _, lab in items]
    x = width - 24 - sum(ws) - text_w("★ 交付结果", 12.5)
    for (k, lab), w in zip(items, ws):
        f, s, _, d = STYLE[k]
        dash = ' stroke-dasharray="5 3"' if d else ""
        o.append(f'<rect x="{x:.1f}" y="30" width="22" height="14" rx="3" fill="{f}" stroke="{s}" stroke-width="1.5"{dash}/>')
        o.append(f'<text x="{x + 28:.1f}" y="42" font-size="12.5" fill="#323F4B">{lab}</text>')
        x += w
    o.append(f'<text x="{x:.1f}" y="42" font-size="12.5" fill="{OUT_COLOR}" font-weight="700">★ 交付结果</text>')

    # stage bands
    for r, name in enumerate(stages):
        fill = "#F7F8FA" if r % 2 == 0 else "#FFFFFF"
        o.append(f'<rect x="12" y="{row_y[r]}" width="{width - 24}" height="{row_h[r]}" fill="{fill}"/>')
        lines = name.split("\n")
        y0 = row_y[r] + row_h[r] / 2 - (len(lines) - 1) * 9 + 5
        for j, ln in enumerate(lines):
            o.append(f'<text x="{GUTTER - 14}" y="{y0 + j * 18:.1f}" text-anchor="end" font-size="14" '
                     f'font-weight="700" fill="#616E7C">{escape(ln)}</text>')
    o.append(f'<line x1="{GUTTER}" y1="{TOP}" x2="{GUTTER}" y2="{height - 20}" stroke="#D9DEE4" stroke-width="1.5"/>')

    # edges (under nodes)
    for e in edges:
        a, b = nodes[e[0]], nodes[e[1]]
        opts = e[2] if len(e) > 2 else {}
        o.append(edge_svg(a, b, opts))

    for n in nodes.values():
        o.append(node_svg(n))
    o.append("</svg>")
    with open(path, "w") as fh:
        fh.write("\n".join(o) + "\n")


def text_w(t, size):
    return sum(size if ord(c) > 0x2E80 else size * 0.56 for c in t)


def port(n, p, off=0.0):
    return {"t": (n.cx + off, n.y), "b": (n.cx + off, n.y + n.h),
            "l": (n.x, n.cy + off), "r": (n.x + n.w, n.cy + off)}[p]


def edge_svg(a, b, opts):
    pa, pb = opts.get("ports", "bt")
    (x1, y1), (x2, y2) = port(a, pa, opts.get("aoff", 0)), port(b, pb, opts.get("boff", 0))
    if pa == "b" and pb == "t":
        if abs(x1 - x2) < 1:
            pts = [(x1, y1), (x2, y2)]
        else:
            my = opts.get("my", y2 - 18)
            pts = [(x1, y1), (x1, my), (x2, my), (x2, y2)]
    elif pa in "lr" and pb == "t":
        pts = [(x1, y1), (x2, y1), (x2, y2)]
    elif pa == "b" and pb in "lr":
        pts = [(x1, y1), (x1, y2), (x2, y2)]
    elif pa in "lr" and pb in "lr":
        if abs(y1 - y2) < 1:
            pts = [(x1, y1), (x2, y2)]
        else:
            mx = opts.get("mx", (x1 + x2) / 2)
            pts = [(x1, y1), (mx, y1), (mx, y2), (x2, y2)]
    else:  # t -> b etc. not used
        pts = [(x1, y1), (x2, y2)]
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    dashed = opts.get("dashed", False)
    color = "#B7791F" if dashed else "#5A6470"
    dash = 'stroke-dasharray="6 4" ' if dashed else ""
    marker = "arro" if dashed else "arr"
    (px, py), (ex, ey) = pts[-2], pts[-1]
    L = max(((ex - px) ** 2 + (ey - py) ** 2) ** 0.5, 1e-6)
    ux, uy = (ex - px) / L, (ey - py) / L
    bx, by = ex - ux * 10, ey - uy * 10
    head = f"{ex:.1f},{ey:.1f} {bx - uy * 5:.1f},{by + ux * 5:.1f} {bx + uy * 5:.1f},{by - ux * 5:.1f}"
    pts[-1] = (bx, by)
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    s = (f'<path d="{d}" fill="none" stroke="{color}" stroke-width="1.8" {dash}/>'
         f'<polygon points="{head}" fill="{color}"/>')
    if "label" in opts:
        # label next to the longest segment
        segs = list(zip(pts, pts[1:]))
        (sx1, sy1), (sx2, sy2) = max(segs, key=lambda s: abs(s[0][0] - s[1][0]) + abs(s[0][1] - s[1][1]))
        mx, my = (sx1 + sx2) / 2, (sy1 + sy2) / 2
        if abs(sx1 - sx2) < 1:  # vertical segment: label on the right
            lx, ly, anchor = mx + 7, my + 4, "start"
        else:
            lx, ly, anchor = mx, my - 6, "middle"
        lx += opts.get("ldx", 0)
        ly += opts.get("ldy", 0)
        tw = text_w(opts["label"], 11.5)
        rx = lx - (tw / 2 if anchor == "middle" else 0) - 3
        s += (f'<rect x="{rx:.1f}" y="{ly - 12:.1f}" width="{tw + 6:.1f}" height="16" rx="3" fill="#FFFFFF" opacity="0.92"/>'
              f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="{anchor}" font-size="11.5" font-style="italic" '
              f'fill="{color}">{escape(opts["label"])}</text>')
    return s


def node_svg(n):
    f, s, tc, d = STYLE[n.kind]
    dash = ' stroke-dasharray="6 4"' if d else ""
    o = [f'<rect x="{n.x:.1f}" y="{n.y:.1f}" width="{n.w:.1f}" height="{n.h:.1f}" rx="8" '
         f'fill="{f}" stroke="{s}" stroke-width="1.8"{dash}/>']
    y = n.y + 10 + 15
    o.append(f'<text x="{n.cx:.1f}" y="{y:.1f}" text-anchor="middle" font-size="15" font-weight="700" '
             f'fill="{tc}">{escape(n.title)}</text>')
    y += 3
    for ln in n.sub:
        y += LINE_S
        o.append(f'<text x="{n.cx:.1f}" y="{y:.1f}" text-anchor="middle" font-size="12" fill="#3E4C59">'
                 f'{escape(ln)}</text>')
    if n.out:
        y += LINE_O + 2
        o.append(f'<text x="{n.cx:.1f}" y="{y:.1f}" text-anchor="middle" font-size="11.5" font-weight="700" '
                 f'font-family="{escape(MONO)}" fill="{OUT_COLOR}">★ {escape(n.out)}</text>')
    return "".join(o)


# =====================================================================
# bulk RNA-seq
# =====================================================================
def rna(path):
    N = Node
    stages = ["输入", "UMI 预处理", "修剪", "比对", "UMI 去重", "链方向\n/ 信号", "定量", "汇总"]
    nodes = [
        N("in", 0, 1.5, "FASTQ", ["samples.tsv：SE / PE、链方向、", "protocol、umi 可按样本设置"], kind="input", span=1.4),
        N("umi", 1, 1.5, "umi_tools extract / 改写 read 名",
          ["umi: read → 从序列提取 UMI", "umi: header（comment）→ UMI 移入 read 名"], kind="opt", span=1.5),
        N("fqc", 2, 0, "FastQC", ["raw / trimmed"], kind="qc"),
        N("trim", 2, 1.5, "Trim Galore", ["--quality 20 --length 30"]),
        N("q3", 2, 3.3, "3' 端修剪", ["切除 read1 5' 端 N nt", "cutadapt 去 polyA"], kind="opt"),
        N("refchk", 3, 0, "参考文件检查", ["index ↔ GTF 染色体命名", "启动时运行"], kind="qc"),
        N("hisat2", 3, 1.5, "HISAT2", ["-t；链方向已知时加", "--rna-strandness", "samtools sort / index"],
          out="03_align/{sample}.bam"),
        N("gtf", 3, 3.3, "GTF 过滤", ["可选，默认不过滤", "如去除 miRNA / snoRNA / snRNA"], kind="opt"),
        N("dedup", 4, 1.5, "umi_tools dedup", ["仅 UMI 样本，之后使用 .dedup.bam"], kind="opt", span=1.3),
        N("bw", 5, 0, "bigwig", ["deeptools bamCoverage", "CPM、bin 100；可分正负链"], out="07_bigwig/"),
        N("strand", 5, 1.5, "链方向检查", ["featureCounts -s 1 / -s 2", "auto 时采用推断值"],
          out="strandedness.tsv"),
        N("fc", 6, 1.5, "featureCounts", ["exon：-t exon -g gene_id", "intron：未分给 exon 的 reads", "再按 -t gene 计数（exon 优先）"],
          out="*_exon.tsv、*_exon_intron.tsv", span=1.35),
        N("st", 6, 3.3, "StringTie", ["TPM / FPKM", "仅 full_length 样本"], kind="opt", out="stringtie_*.tsv"),
        N("rep", 6, 0, "repeat 计数", ["featureCounts", "RepeatMasker GTF（自备）"], kind="opt",
          out="repeats_counts.tsv"),
        N("mqc", 7, 1.65, "summary.tsv + MultiQC", ["reads 数、比对率、UMI 重复率、assigned 率、intron 率、链方向；汇总所有 QC 结果"],
          kind="report", span=4.3, out="tables/summary.tsv、qc/multiqc/multiqc_report.html"),
    ]
    e = [
        ("in", "umi", {}),
        ("umi", "trim", {}),
        ("trim", "fqc", {"ports": "lr"}),
        ("q3", "trim", {"ports": "lr", "dashed": True, "label": "protocol: 3prime"}),
        ("trim", "hisat2", {}),
        ("hisat2", "dedup", {}),
        ("dedup", "strand", {}),
        ("dedup", "bw", {"ports": "bt", "aoff": -50}),
        ("gtf", "strand", {"ports": "br", "label": "genes.gtf（定量共用）"}),
        ("strand", "fc", {"label": "-s 0 / 1 / 2"}),
        ("strand", "rep", {"ports": "bt", "aoff": -80, "dashed": True}),
        ("strand", "st", {"ports": "bt", "aoff": 80, "dashed": True}),
    ]
    render("snakemake_bulk_rna", "bulk RNA-seq：fastq → 比对 → 链方向 → 基因定量 / bigwig",
           stages, nodes, e, 4.3, path)


# =====================================================================
# bulk ChIP-seq / ATAC-seq（结构相同）
# =====================================================================
def peaks(path, atac):
    N = Node
    stages = ["输入", "修剪", "比对", "过滤去重", "call peak\n/ 信号", "后处理", "motif", "汇总"]
    if atac:
        inp = ["samples.tsv：sample、fq1、fq2", "SE / PE 自动识别"]
        trim_sub = ["--nextera（Tn5 接头）"]
        bt2_sub = ["--very-sensitive -X 2000", "PE：--no-mixed --no-discordant"]
        flt_sub = ["MAPQ ≥ 30、properly paired", "去除 chrM", "samtools markdup -r"]
        macs_sub = ["PE：-f BAMPE", "SE：--shift -100 --extsize 200", "narrowPeak，q < 0.01"]
        stats_sub = ["片段长度 / 核小体周期", "chrM 比例、flagstat"]
        homer_sub = ["summit ± 100 bp"]
        sumtxt = ["reads 数、比对率、chrM 比例、重复率、peak 数、FRiP、TSS enrichment；汇总所有 QC 结果"]
        title, sub = "snakemake_bulk_atac", "bulk ATAC-seq：fastq → bowtie2 → 过滤去重 → MACS3 → peak / bigwig"
    else:
        inp = ["samples.tsv：sample、fq1、fq2、", "control（input / IgG）、peak_type"]
        trim_sub = ["自动识别接头"]
        bt2_sub = ["-N 1", "PE：--no-mixed --no-discordant"]
        flt_sub = ["MAPQ ≥ 30、properly paired", "chrM 可选去除", "samtools markdup -r"]
        macs_sub = ["有对照：-c control", "无对照：--nolambda --nomodel", "narrow / broad；PE 用 BAMPE"]
        stats_sub = ["片段长度、chrM 比例", "flagstat / idxstats"]
        homer_sub = ["summit ± 100 bp", "仅 narrow、非 control 样本"]
        sumtxt = ["reads 数、比对率、重复率、peak 数、FRiP、TSS enrichment；汇总所有 QC 结果"]
        title, sub = "snakemake_bulk_chip", "ChIP-seq / CUT&Tag：fastq → bowtie2 → 过滤去重 → MACS3 → peak / bigwig"
    nodes = [
        N("in", 0, 1.5, "FASTQ", inp, kind="input", span=1.4),
        N("fqc", 1, 0, "FastQC", ["raw / trimmed"], kind="qc"),
        N("trim", 1, 1.5, "Trim Galore", trim_sub),
        N("refchk", 2, 0, "参考文件检查", ["index / blacklist 染色体命名", "GTF / fasta 自动转换命名"], kind="qc"),
        N("bt2", 2, 1.5, "Bowtie2", bt2_sub),
        N("contam", 2, 3.2, "污染检查", ["比对到其它物种（如 hg38）", "只报告比对率"], kind="opt"),
        N("stats", 3, 0, "samtools stats", stats_sub, kind="qc"),
        N("flt", 3, 1.5, "过滤 + 去重", flt_sub, out="02_align/{sample}.dedup.bam", span=1.25),
        N("macs", 4, 1.5, "MACS3", macs_sub),
        N("bw", 4, 3.2, "bigwig", ["deeptools bamCoverage", "CPM、bin 100、排除 blacklist"], out="04_bigwig/"),
        N("frip", 5, 0, "FRiP", ["peak 内 reads 比例"], kind="qc"),
        N("bl", 5, 1.5, "去除 blacklist", ["bedtools，ENCODE v2"], out="03_peaks/"),
        N("tss", 5, 3.2, "TSS profile", ["TSS 来自 GTF", "TSS enrichment"], kind="qc"),
        N("homer", 6, 1.5, "HOMER motif", homer_sub, kind="opt", out="05_motif/"),
        N("mqc", 7, 1.6, "summary.tsv + MultiQC", sumtxt, kind="report", span=4.2,
          out="tables/summary.tsv、qc/multiqc/multiqc_report.html"),
    ]
    if not atac:
        nodes.append(N("fp", 4, 0, "fingerprint", ["deeptools", "IP 相对 input 的富集"], kind="qc"))
    e = [
        ("in", "trim", {}),
        ("trim", "fqc", {"ports": "lr"}),
        ("trim", "bt2", {}),
        ("trim", "contam", {"ports": "rt", "dashed": True}),
        ("bt2", "flt", {}),
        ("flt", "stats", {"ports": "lr"}),
        ("flt", "macs", {"label": "" if atac else "IP 样本（+ control BAM）"}),
        ("flt", "bw", {"ports": "rt", "label": "所有样本" if not atac else ""}),
        ("macs", "bl", {}),
        ("bl", "frip", {"ports": "lr"}),
        ("bw", "tss", {}),
        ("bl", "homer", {"dashed": True, "label": "narrowPeak"}),
    ]
    e = [(a, b, {k: v for k, v in o.items() if not (k == "label" and not v)}) for a, b, o in e]
    if not atac:
        e.append(("flt", "fp", {"ports": "bt", "aoff": -50}))
    render(title, sub, stages, nodes, e, 4.2, path)


def to_png(svg):
    if shutil.which("rsvg-convert"):
        subprocess.run(["rsvg-convert", "-z", "2", svg, "-o", svg[:-4] + ".png"], check=True)
    else:
        print("未找到 rsvg-convert，跳过 PNG", file=sys.stderr)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    jobs = [("rna", rna), ("chip", lambda p: peaks(p, atac=False)), ("atac", lambda p: peaks(p, atac=True))]
    for name, fn in jobs:
        out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", f"snakemake_bulk_{name}", "docs", "images")
        os.makedirs(out, exist_ok=True)
        svg = os.path.join(out, f"bulk_{name}_workflow.svg")
        fn(svg)
        to_png(svg)
        print(os.path.normpath(svg))
