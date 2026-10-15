#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Usage: nohup bash download_sra_iseq.sh <metadata_file> > download.log 2>&1 &
#
# metadata: 带表头，制表符或空白分隔；按列名取值，列顺序任意，多余列忽略。
# 必需列：run_accession, sample_title
#   run_accession  study_accession  library_strategy  sample_title
#   SRR22796946    PRJNA912643      ATAC-seq          Brain-endothelial_O_1
#
# 同一 sample_title 对应多个 run（如 Neurons_O_3）时，按出现顺序改名为
# <title>_batch1, <title>_batch2 ...，再由 merge_all_v2.py 合并为 <title>。

PYTHON_BIN="python"
MERGE_SCRIPT="/home/LAB/lab_ab_1/pipeline/merge_all_v3.py"
ISEQ_OPTS=(-g)

log()  { echo "[$(date '+%F %T')] [INFO] $*"; }
warn() { echo "[$(date '+%F %T')] [WARN] $*" >&2; }
die()  { echo "[$(date '+%F %T')] [ERROR] $*" >&2; exit 1; }

###############################################################################
# 参数检查
###############################################################################

[ $# -ge 1 ] || { echo "Usage: $0 <metadata_file>"; exit 1; }

metadata="$1"
[ -f "$metadata" ] || die "Metadata file '$metadata' not found."

metadata_dir=$(dirname "$metadata")
metadata_base=$(basename "$metadata")
cd "$metadata_dir" || die "Failed to cd to $metadata_dir"

log "Working dir: $PWD"
log "Metadata: $metadata_base"

command -v iseq >/dev/null || die "iseq not found in PATH."

###############################################################################
# 解析 metadata → run_map.tsv（run_accession <TAB> 最终文件名前缀）
###############################################################################

run_map="run_map.tsv"

# tr 去掉 Windows 换行；按表头列名定位，不依赖列顺序
tr -d '\r' < "$metadata_base" | awk '
    function trim(x) { gsub(/^[ \t]+|[ \t]+$/, "", x); return x }
    NF == 0 { next }
    !hdr {
        hdr = 1
        if (index($0, "\t")) FS = "\t"          # 有 tab 用 tab，否则按空白分隔
        for (i = 1; i <= NF; i++) col[trim($i)] = i
        if (!("run_accession" in col) || !("sample_title" in col)) {
            print "missing column run_accession or sample_title in header" > "/dev/stderr"
            exit 2
        }
        rc = col["run_accession"]; sc = col["sample_title"]
        next
    }
    {
        run = trim($rc); title = trim($sc)
        if (run == "" || title == "") { print "skip incomplete line: " $0 > "/dev/stderr"; next }
        if (run !~ /^[SDE]RR[0-9]+$/) print "unusual run_accession: " run > "/dev/stderr"
        if (run in seen) { print "duplicate run_accession skipped: " run > "/dev/stderr"; next }
        gsub(/[^A-Za-z0-9._-]/, "_", title)      # 文件名安全
        seen[run] = 1
        n++; runs[n] = run; titles[n] = title; cnt[title]++
    }
    END {
        if (!hdr) exit 2
        for (i = 1; i <= n; i++) {
            t = titles[i]
            if (cnt[t] > 1) t = t "_batch" (++idx[t])
            print runs[i] "\t" t
        }
    }
' > "$run_map" || die "Failed to parse metadata."

[ -s "$run_map" ] || die "No runs found in metadata."
cut -f1 "$run_map" > sra_list.txt
log "$(wc -l < sra_list.txt) run(s) -> sra_list.txt, mapping -> $run_map"

if grep -q '_batch[0-9]' "$run_map"; then
    [ -f "$MERGE_SCRIPT" ] || die "Need merging but merge script not found: $MERGE_SCRIPT"
    log "Samples with multiple runs (will be merged):"
    grep '_batch[0-9]' "$run_map" | sed 's/^/        /'
fi

# 某个 run 的所有 fastq：SRRxxx.fastq.gz / SRRxxx_1.fastq.gz / SRRxxx_2.fastq.gz
# 精确匹配，避免 SRR123* 误匹配 SRR1234
run_files() {
    files=()
    local f
    for f in "$1".fastq.gz "$1"_[0-9].fastq.gz; do
        [ -e "$f" ] && files+=( "$f" )
    done
    return 0
}

###############################################################################
# 下载
###############################################################################

# 旧的 fail.log 会导致重跑时误判失败，先归档
[ -f fail.log ] && mv fail.log "fail.log.$(date +%Y%m%d_%H%M%S)"

log "Running iseq ${ISEQ_OPTS[*]} ..."
iseq -i sra_list.txt "${ISEQ_OPTS[@]}" >> iseq.log 2>&1 || die "iseq failed. Check iseq.log."
log "iseq finished."

if [ -f fail.log ]; then
    cat fail.log >&2
    die "fail.log detected. Some downloads failed; fix and rerun (finished runs are skipped by iseq)."
fi

###############################################################################
# 检查每个 run 是否下载完整
###############################################################################

log "Checking fastq.gz for each run ..."
missing=0
while IFS=$'\t' read -r run title; do
    run_files "$run"
    if [ ${#files[@]} -eq 0 ]; then
        warn "No FASTQ found for $run ($title)"
        missing=$((missing + 1))
    elif [ -e "${run}_1.fastq.gz" ] && [ ! -e "${run}_2.fastq.gz" ]; then
        warn "$run has _1 but no _2"
        missing=$((missing + 1))
    fi
done < "$run_map"
[ "$missing" -eq 0 ] || die "$missing run(s) incomplete; aborting before rename/merge."

###############################################################################
# MD5（重命名前）
###############################################################################

log "Computing MD5 before renaming ..."
raw_files=()
while IFS=$'\t' read -r run _; do
    run_files "$run"
    raw_files+=( "${files[@]}" )
done < "$run_map"
md5sum "${raw_files[@]}" > raw.md5sum
log "raw.md5sum written (${#raw_files[@]} files)."

###############################################################################
# 重命名：run_accession 前缀 → sample_title
###############################################################################

log "Renaming files (run_accession -> sample_title) ..."
while IFS=$'\t' read -r run title; do
    run_files "$run"
    for file in "${files[@]}"; do
        new_name="${title}${file#"$run"}"
        if [ -e "$new_name" ]; then
            warn "Target '$new_name' already exists, skip '$file'."
            continue
        fi
        mv -- "$file" "$new_name"
        log "$file -> $new_name"
    done
done < "$run_map"

###############################################################################
# 合并 batch
###############################################################################

batch_files=( *_batch[0-9]*.fastq.gz )
if [ ${#batch_files[@]} -gt 0 ]; then
    mkdir -p unmerged
    log "Moving ${#batch_files[@]} batch file(s) to unmerged/ ..."
    mv -- "${batch_files[@]}" unmerged/

    log "Running merge_all_v2.py ..."
    "$PYTHON_BIN" "$MERGE_SCRIPT" \
        --type filename-based \
        -i unmerged/ \
        -o ./ \
        --remove-pattern '_batch[0-9]' \
        || die "merge_all_v2.py failed."
    log "Merging done."
else
    log "No *_batch* files, skip merging."
fi

###############################################################################
# 最终检查 + MD5
###############################################################################

log "Verifying final files ..."
lost=0
while read -r sample; do
    run_files "$sample"
    if [ ${#files[@]} -eq 0 ]; then
        warn "No final FASTQ for sample: $sample"
        lost=$((lost + 1))
    fi
done < <(cut -f2 "$run_map" | sed 's/_batch[0-9]*$//' | awk '!s[$0]++')
[ "$lost" -eq 0 ] || die "$lost sample(s) missing final FASTQ."

final_files=( ./*.fastq.gz )
md5sum "${final_files[@]}" > final.md5sum
log "final.md5sum written (${#final_files[@]} files)."

log "All done."
