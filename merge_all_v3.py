#!/usr/bin/env python3
# eg: python merge_all_v2.py --type filename-based -i ../rawdata/unmerged -o ./ --suffix .qc_for_fastq_gz.sh --remove-pattern="_batch[0-9]"

import argparse
import os
import sys
import re
import subprocess
import logging
from pathlib import Path
from collections import defaultdict


def setup_logging():
    """配置日志格式和级别。"""
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def parse_args():
    parser = argparse.ArgumentParser(description="Merge 指定后缀文件 with validation. ")
    parser.add_argument("--type", required=True, choices=["dir-based", "filename-based"],
                        help="Merge mode: dir-based or filename-based")
    parser.add_argument("-i", "--input", required=True, help="Input directory")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--suffix", default=".qc_for_fastq_gz.sh",
                        help="File suffix to merge (e.g. .qc_for_fastq_gz.sh, .fq.gz) [default: %(default)s]")
    parser.add_argument("--r1-id", default="_R1",
                        help="Identifier for R1 files (dir-based only) [default: %(default)s]")
    parser.add_argument("--r2-id", default="_R2",
                        help="Identifier for R2 files (dir-based only) [default: %(default)s]")
    parser.add_argument("--remove-pattern",
                        help="String or regex to remove from filenames (filename-based only eg: _rep[0-9])")
    parser.add_argument("-n", "--dry-run", action="store_true", help="Dry run: print commands without executing")
    parser.add_argument("--skip-lines", action="store_true",default=True,
                        help="Skip line counting (useful for huge files to save time)")
    return parser.parse_args()


def format_size(size_in_bytes):
    """将字节数转化为人类可读的格式 (类似于 ls -h)"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_in_bytes < 1024.0:
            return f"{size_in_bytes:.2f} {unit}"
        size_in_bytes /= 1024.0
    return f"{size_in_bytes:.2f} PB"


def count_lines(filepath):
    """智能统计行数：根据扩展名决定是否使用 zcat"""
    filepath_str = str(filepath)
    # 增加对普通文本文件的兼容，防止非 gz 文件使用 zcat 报错
    if filepath_str.endswith('.gz'):
        cmd = f"zcat '{filepath_str}' | wc -l"
    else:
        cmd = f"wc -l < '{filepath_str}'"

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    return int(result.stdout.strip())


def merge_files(input_files, output_file, dry_run=False, skip_lines=False):
    logging.info(f"Merging {len(input_files)} files to: {output_file.name}")

    # 提取所有合并前文件的绝对路径以便记录
    input_paths = [str(f.resolve()) for f in input_files]

    if dry_run:
        logging.info(f"[DRY-RUN] Would run: cat {' '.join(input_paths)} > {output_file}")
        return

    # 1. 记录合并前的文件大小
    before_sizes = [os.path.getsize(f) for f in input_files]
    total_before_size = sum(before_sizes)

    # 2. 执行合并
    # 使用标准库进行安全的 cat 操作，完美规避带空格的文件名导致的 shell 报错
    with open(output_file, 'wb') as outfile:
        # 仅调用 cat 命令，通过 stdout 写入文件，更加安全高效
        subprocess.run(["cat"] + input_paths, stdout=outfile, check=True)

    # 3. 记录合并后的文件大小
    after_size = os.path.getsize(output_file)

    # 4. 大小与文件展示
    logging.info("File Size Check:")
    for f, s in zip(input_files, before_sizes):
        logging.info(f"   [Input]  {format_size(s):>10}  {f.name}")
    logging.info(f"   [Output] {format_size(after_size):>10}  {output_file.name}")

    if after_size != total_before_size:
        logging.error(f"× Size mismatch! Expected {format_size(total_before_size)} but got {format_size(after_size)}")
    else:
        logging.info("√ File size strictly matches (Output == sum of Inputs)!")
    logging.info("-" * 50)

    # 5. 验证行数一致性 (可选，如果是超大文件，可能耗时较长)
    if not skip_lines:
        logging.info(f"Counting lines for validation... (This may take a while for large files)")
        before_lines = sum(count_lines(f) for f in input_files)
        after_lines = count_lines(output_file)
        if before_lines != after_lines:
            logging.error(f"× Line count mismatch: before={before_lines}, after={after_lines}")
            sys.exit(1)
        else:
            logging.info(f"√ Line count matched: {after_lines} lines")

def run_dir_based(input_dir, output_dir, r1_id, r2_id, suffix, dry_run, skip_lines):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for subdir in sorted(input_dir.iterdir()):
        if not subdir.is_dir():
            continue

        sample = subdir.name
        all_files = sorted(subdir.glob(f"*{suffix}"))
        r1_files = [f for f in all_files if r1_id in f.name]
        r2_files = [f for f in all_files if r2_id in f.name]
        r1_set = set(r1_files)
        r2_set = set(r2_files)
        paired_set = r1_set.union(r2_set)

        if r1_files:
            out_r1 = output_dir / f"{sample}_R1{suffix}"
            merge_files(r1_files, out_r1, dry_run=dry_run, skip_lines=skip_lines)

        if r2_files:
            out_r2 = output_dir / f"{sample}_R2{suffix}"
            merge_files(r2_files, out_r2, dry_run=dry_run, skip_lines=skip_lines)

        unpaired_files = [f for f in all_files if f not in paired_set]
        if not r1_files and not r2_files and unpaired_files:
            logging.info(f"SE mode detected for {sample}")
            out_se = output_dir / f"{sample}{suffix}"
            merge_files(unpaired_files, out_se, dry_run=dry_run, skip_lines=skip_lines)
        elif not all_files:
            logging.warning(f"No files found in {subdir.name} with suffix {suffix}")


def run_filename_based(input_dir, output_dir, suffix, remove_pattern, dry_run, skip_lines):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    groups = defaultdict(list)
    for file in sorted(input_dir.glob(f"*{suffix}")):
        name = file.name
        key = name
        if remove_pattern:
            key = re.sub(remove_pattern, "", name)
        groups[key].append(file)

    for key, files in groups.items():
        if len(files) == 1:
            logging.info(f"Only 1 file found for {key}, skipping merge.")
            continue

        output_file = output_dir / key
        merge_files(files, output_file, dry_run=dry_run, skip_lines=skip_lines)


def main():
    setup_logging()
    args = parse_args()

    logging.info("============= ^.^ Start merge process =============")
    logging.info(f"Mode: {args.type}, Suffix: {args.suffix}")

    if args.type == "dir-based":
        run_dir_based(args.input, args.output, args.r1_id, args.r2_id, args.suffix, args.dry_run, args.skip_lines)
    elif args.type == "filename-based":
        run_filename_based(args.input, args.output, args.suffix, args.remove_pattern, args.dry_run, args.skip_lines)
    else:
        logging.error("Unknown merge type.")
        sys.exit(1)

    logging.info("All operations completed successfully!")


if __name__ == "__main__":
    main()