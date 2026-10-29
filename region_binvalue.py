#!/usr/bin/env python3
"""
region_binvalue.py

Region Bin Value Calculator - 计算区域内分箱加权信号值

核心功能：
    以 Region 为索引，通过与 Value File 的 Overlap 比例，计算加权平均信号强度。
    输入 Value 视为"均匀分布的信号强度 (Uniform Signal Density)"。

计算逻辑（模式 B：未覆盖区域不参与平均值计算）：
    - BigWig:  sum(valid_values) / count(valid_values)，使用逐碱基精确计算，NaN 区域跳过
    - BED: sum(overlap * value) / sum(overlap)，无交集区域跳过
    - 完全无覆盖时输出 NA

数据预处理：
    1) 限制范围：过滤非法染色体，裁剪超出染色体边界的坐标
    2) 区间扩展：支持定长 (--extend-bp) 或比例扩展 (--extend-ratio)
       注意：比例扩展基于裁剪后的区间长度计算
    3) 分 Bin：将区间划分为固定数量的子窗口，引入空间信息

输出格式：
    Wide Format TSV，适配 R/ggplot2 热图绘制
    列结构：chrom, start, end, bin_id, region_uid, <retain_cols... >, <value_files...>

Requirements:
    - Python 3.8+
    - pandas, numpy, pybedtools
    - pyBigWig (optional, required for bigWig files)
    - bedtools installed and in PATH

Usage:
    # 定长扩展 (上下游各扩展 2000bp)
    python region_binvalue.py \\
        --region regions.bed \\
        --values a.bed,b.bw \\
        --chrom_size mm10.chrom.sizes \\
        --value_col 3 \\
        --binnum 100 \\
        --extend-bp 2000 \\
        --retain-cols "3,4" \\
        --outdir ./out \\
        --threads 4

    # 比例扩展 (上下游各扩展 0.5 倍区间长度)
    python region_binvalue.py \\
        --region regions.bed \\
        --values a.bed,b.bw \\
        --chrom_size mm10.chrom.sizes \\
        --value_col 3 \\
        --binnum 100 \\
        --extend-ratio 0.5 \\
        --retain-cols "3,4" \\
        --save-extended-region \\
        --outdir ./out \\
        --threads 4

Output:
    <region_basename>_binvalue.tsv
    <region_basename>_extended.bed (if --save-extended-region)
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from pybedtools import BedTool

# Try to import pyBigWig; if not installed, we'll raise when needed
try:
    import pyBigWig
except ImportError:
    pyBigWig = None

# -----------------------
# Defaults
# -----------------------
DEFAULT_CHROM_SIZE = "/home_data/home/slst/leixy2023/data/database/mm10/archive/mm10_filtered.chrom.sizes"
DEFAULT_MAX_BINS_PER_BATCH = 2_000_000


# -----------------------
# Logging Setup
# -----------------------
def setup_logging(verbose: bool = False) -> None:
    """
    Configure global logging settings.

    Args:
        verbose:  If True, set logging level to DEBUG; otherwise INFO.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# -----------------------
# Argument Parsing
# -----------------------
def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parse CLI arguments and return Namespace.

    Args:
        argv: Command line arguments.  If None, uses sys.argv.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(
        description="Region fixed binnum weighted-sum calculation with support for BED and BigWig.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples: 
    region的前三列一定是chr,start,end
    .bed的计算会另外调用bedtools，i/o会稍微慢一些，建议使用.bw文件以提升效率。

    # Fixed bp extension (2000bp on each side)
    python region_binvalue.py --region peaks.bed --values signal.bw \\
        --chrom_size mm10.sizes --value_col 3 --binnum 100 --extend-bp 2000

    # Proportional extension (0.5x region length on each side)
    python region_binvalue.py --region peaks.bed --values signal.bw \\
        --chrom_size mm10.sizes --value_col 3 --binnum 100 --extend-ratio 0.5
        """,
    )
    parser.add_argument(
        "--region",
        required=True,
        help="Region BED file(s), comma or space separated.",
    )
    parser.add_argument(
        "--values",
        required=True,
        help="Value files (BED or BigWig), comma or space separated.",
    )
    parser.add_argument(
        "--chrom_size",
        default=DEFAULT_CHROM_SIZE,
        help=f"Chromosome sizes file (default: {DEFAULT_CHROM_SIZE}).",
    )
    parser.add_argument(
        "--value_col",
        type=int,
        help="0-based column index for value in BED files.",
    )
    parser.add_argument(
        "--extend-bp",
        type=int,
        default=None,
        help="Extend region by fixed bp on each side. Mutually exclusive with --extend-ratio.",
    )
    parser.add_argument(
        "--extend-ratio",
        type=float,
        default=None,
        help="Extend region by ratio of clipped region length on each side. "
        "E.g., 0.5 means extend 0.5x length; 3 means extend 3x length. "
        "Mutually exclusive with --extend-bp.",
    )
    parser.add_argument(
        "--binnum",
        type=int,
        required=True,
        help="Number of bins per region.",
    )
    parser.add_argument(
        "--retain-cols",
        default="",
        help="Comma-separated 0-based extra region columns to retain (>=3).",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Worker processes to parallelize over value files.",
    )
    parser.add_argument(
        "--max-bins-per-batch",
        type=int,
        default=DEFAULT_MAX_BINS_PER_BATCH,
        help=f"Max bins per batch to limit memory (default: {DEFAULT_MAX_BINS_PER_BATCH}).",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Do not write header line to output file.",
    )
    parser.add_argument(
        "--save-extended-region",
        action="store_true",
        help="Save extended region file with retained columns.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Set logging level to DEBUG.",
    )

    args = parser.parse_args(argv)

    # Validate mutually exclusive extend parameters
    if args.extend_bp is not None and args.extend_ratio is not None:
        parser.error(
            "--extend-bp and --extend-ratio are mutually exclusive.  Please specify only one."
        )

    # Validate extend values
    if args.extend_bp is not None and args.extend_bp < 0:
        parser.error("--extend-bp must be >= 0.")
    if args.extend_ratio is not None and args.extend_ratio < 0:
        parser.error("--extend-ratio must be >= 0.")

    return args


# -----------------------
# Utilities
# -----------------------
def parse_multi_files(arg: str) -> List[str]:
    """Parse comma or whitespace separated file lists."""
    if "," in arg:
        files = [f.strip() for f in arg.split(",") if f.strip()]
    else:
        files = [f.strip() for f in arg.split() if f.strip()]
    return files


def read_bed_file_auto_sep(bed_file: str) -> pd.DataFrame:
    """
    Read a BED-like file trying common separators.

    Args:
        bed_file: Path to BED file.

    Returns:
        DataFrame with BED data (no header).

    Raises:
        ValueError: If file cannot be parsed.
    """
    for sep in ["\t", ",", " "]:
        try:
            if sep == " ":
                df = pd.read_csv(bed_file, delim_whitespace=True, header=None)
            else:
                df = pd.read_csv(bed_file, sep=sep, header=None)
            if df.shape[1] >= 3:
                return df
        except Exception:
            continue
    raise ValueError(f"Cannot parse BED file separators: {bed_file}")


def read_chrom_size_file(chrom_size_file: str) -> pd.DataFrame:
    """
    Read chromosome sizes into DataFrame.

    Args:
        chrom_size_file: Path to chromosome sizes file.

    Returns:
        DataFrame with columns ['chrom', 'size'].
    """
    for sep in ["\t", ",", " "]:
        try:
            if sep == " ":
                df = pd.read_csv(chrom_size_file, delim_whitespace=True, header=None)
            else:
                df = pd.read_csv(chrom_size_file, sep=sep, header=None)
            if df.shape[1] >= 2:
                df = df.iloc[:, :2].copy()
                df.columns = ["chrom", "size"]
                df["size"] = df["size"].astype(int)
                return df
        except Exception:
            continue
    raise ValueError(f"Cannot parse chrom_size file: {chrom_size_file}")


def output_filename(region_file: str) -> str:
    """Generate output filename from region file."""
    base = Path(region_file).stem
    return f"{base}_binvalue.tsv"


def extended_region_filename(region_file: str) -> str:
    """Generate extended region filename from region file."""
    base = Path(region_file).stem
    return f"{base}_extended.bed"


def _safe_remove(path: Optional[str]) -> None:
    """Safely remove a file if it exists."""
    try:
        if path and os.path.exists(path):
            os.remove(path)
    except Exception:
        pass


def is_bigwig_file(filepath: str) -> bool:
    """Check if file is a BigWig file based on extension."""
    return filepath.lower().endswith((".bw", ".bigwig"))


# -----------------------
# Core Classes
# -----------------------
class RegionProcessor:
    """
    Read/filter/extend regions and extract extra columns.

    区域预处理器：
    1. 过滤非法染色体
    2. 裁剪超出边界的坐标
    3. 扩展区间（支持定长 bp 或比例模式）
    4. 提取用户指定的元数据列
    """

    def __init__(self, chrom_sizes:  pd.DataFrame):
        """
        Initialize RegionProcessor.

        Args:
            chrom_sizes: DataFrame with columns ['chrom', 'size'].
        """
        self.chrom_sizes = chrom_sizes
        self.size_map:  Dict[str, int] = dict(
            zip(self.chrom_sizes["chrom"], self.chrom_sizes["size"])
        )

    def filter_and_clip(self, region_df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter to valid chromosomes and clip coordinates to chromosome boundaries.

        过滤非法染色体，并将坐标裁剪到染色体范围内。

        Args:
            region_df: Input region DataFrame.

        Returns:
            Filtered and clipped DataFrame.
        """
        # Filter to valid chromosomes
        valid_chroms = set(self.chrom_sizes["chrom"])
        df = region_df[region_df.iloc[:, 0].isin(valid_chroms)].copy()

        if df.empty:
            return df.reset_index(drop=True)

        # Clip start to >= 0
        df.iloc[:, 1] = df.iloc[:, 1].astype(int).clip(lower=0)

        # Clip end to <= chrom_size
        sizes = df.iloc[:, 0].map(self.size_map).astype(int)
        df.iloc[:, 2] = np.minimum(df.iloc[:, 2].astype(int), sizes)

        # Filter out invalid regions (start >= end)
        df = df[df.iloc[:, 2] > df.iloc[:, 1]]

        return df.reset_index(drop=True)

    def extend_regions(
        self,
        region_df: pd.DataFrame,
        extend_bp: Optional[int] = None,
        extend_ratio: Optional[float] = None,
    ) -> pd.DataFrame:
        """
        Extend each region by fixed bp or by ratio.

        区间扩展逻辑：
        - 定长模式 (extend_bp): 上下游各扩展 extend_bp bp
        - 比例模式 (extend_ratio): 上下游各扩展 (裁剪后区间长度 × extend_ratio) bp

        注意：扩展基于裁剪后的区间长度计算，扩展后会再次裁剪到染色体边界。

        Args:
            region_df: Input region DataFrame (already filtered and clipped).
            extend_bp: Fixed bp to extend on each side (mutually exclusive with extend_ratio).
            extend_ratio: Ratio of region length to extend on each side.

        Returns:
            Extended DataFrame.
        """
        # No extension if both are None or zero
        if extend_bp is None and extend_ratio is None:
            return region_df.copy()
        if extend_bp == 0 or extend_ratio == 0:
            return region_df.copy()

        df = region_df.copy()
        starts = df.iloc[:, 1].astype(int)
        ends = df.iloc[:, 2].astype(int)

        if extend_ratio is not None:
            # Ratio mode: extend by (clipped_length * ratio) on each side
            lengths = ends - starts
            deltas = (lengths * extend_ratio).astype(int)

            # Calculate new coordinates
            new_starts = starts - deltas
            new_ends = ends + deltas
        else:
            # Fixed bp mode
            deltas = np.full(len(df), extend_bp, dtype=int)

            # Calculate new coordinates
            center = (starts + ends) // 2
            new_starts = center - deltas
            new_ends = center + deltas

        # Clip to chromosome boundaries
        new_starts = new_starts.clip(lower=0)
        chrom_sizes = df.iloc[:, 0].map(self.size_map).astype(int)
        new_ends = np.minimum(new_ends, chrom_sizes)

        df.iloc[:, 1] = new_starts
        df.iloc[:, 2] = new_ends

        # Filter out invalid regions after extension
        df = df[df.iloc[:, 2] > df.iloc[:, 1]]

        return df.reset_index(drop=True)

    @staticmethod
    def extract_metadata(
            region_df: pd.DataFrame,
            retain_cols: List[int]
    ) -> Optional[pd.DataFrame]:
        """
        Extract extra columns to keep.

        Args:
            region_df:  Input region DataFrame.
            retain_cols:  List of 0-based column indices to retain.

        Returns:
            DataFrame with 'region_uid' and retained columns, or None if no columns to retain.
        """
        if not retain_cols:
            return None

        n = region_df.shape[0]
        out = pd.DataFrame({"region_uid": list(range(n))})

        for c in retain_cols:
            colname = f"region_col_{c}"
            if c < region_df.shape[1]:
                out[colname] = region_df.iloc[:, c].values
            else:
                out[colname] = pd.NA

        return out


class BinCalculator:
    """
    Integer-safe binning operations.

    分箱计算器：将区间划分为固定数量的子窗口。
    使用"基数+余数"算法确保整数精度。
    """

    @staticmethod
    def split_region_to_bins(start: int, end: int, binnum: int) -> List[Tuple[int, int]]:
        """
        Split [start, end) into `binnum` integer bins.

        使用基数+余数算法：前 rem 个 bin 长度为 base+1，其余为 base。

        Args:
            start: Region start coordinate.
            end: Region end coordinate.
            binnum: Number of bins.

        Returns:
            List of (bin_start, bin_end) tuples.
        """
        length = int(end) - int(start)
        if length <= 0 or binnum < 1:
            return []

        base = length // binnum
        rem = length % binnum

        # Allocate bin lengths
        lens = np.empty(binnum, dtype=int)
        if rem > 0:
            lens[:rem] = base + 1
            lens[rem:] = base
        else:
            lens[:] = base

        # Calculate bin coordinates
        starts = np.empty(binnum, dtype=int)
        starts[0] = start
        if binnum > 1:
            starts[1:] = start + np.cumsum(lens[:-1])
        ends = starts + lens

        return list(zip(starts.tolist(), ends.tolist()))

    @staticmethod
    def regions_df_to_bins_df(region_df: pd.DataFrame, binnum: int) -> pd.DataFrame:
        """
        Convert region DataFrame into bin DataFrame.

        Args:
            region_df: DataFrame with columns [chrom, start, end, ... ].
            binnum: Number of bins per region.

        Returns:
            DataFrame with columns [chrom, start, end, bin_id, region_uid].
            - region_uid: Global index referring to original region row.
            - bin_id: Bin index within the region (0 to binnum-1).
        """
        rows = []
        for idx, row in region_df.iterrows():
            chrom = row.iloc[0]
            s = int(row.iloc[1])
            e = int(row.iloc[2])
            for bin_id, (bs, be) in enumerate(
                BinCalculator.split_region_to_bins(s, e, binnum)
            ):
                rows.append([chrom, int(bs), int(be), int(bin_id), int(idx)])

        if not rows:
            return pd.DataFrame(
                columns=["chrom", "start", "end", "bin_id", "region_uid"]
            )

        return pd.DataFrame(
            rows, columns=["chrom", "start", "end", "bin_id", "region_uid"]
        )


class BigWigRunner:
    """
    Wrapper for pyBigWig access.

    BigWig 处理器（模式 B）：使用 bw.values() 逐碱基精确计算。

    计算逻辑：
    1. 使用 bw.values() 获取区间内每个碱基的值
    2. 过滤 NaN（未覆盖区域）
    3. 计算 sum(valid_values) / count(valid_values)
    4. 如果完全无覆盖，返回 NA
    """

    def __init__(self, bw_path: str):
        """
        Initialize BigWigRunner.

        Args:
            bw_path: Path to BigWig file.

        Raises:
            ImportError: If pyBigWig is not installed.
            IOError: If file cannot be opened.
        """
        if pyBigWig is None:
            raise ImportError(
                "pyBigWig is required to read BigWig files.  "
                "Install with `pip install pyBigWig`."
            )

        self.bw_path = bw_path
        try:
            self.bw = pyBigWig.open(bw_path)
        except (RuntimeError, IOError) as e:
            raise IOError(f"Failed to open BigWig file {bw_path}:  {e}")

    def query_bins_density(self, bin_df: pd.DataFrame) -> List[Optional[float]]:
        """
        Calculate average signal density for each bin (模式 B).

        使用 bw.values() 逐碱基计算，未覆盖区域不参与平均值计算。
        计算公式:  sum(valid_values) / count(valid_values)

        Args:
            bin_df: DataFrame with columns [chrom, start, end].

        Returns:
            List of density values (float) or None if no coverage.
        """
        densities:  List[Optional[float]] = []

        for _, row in bin_df.iterrows():
            chrom = str(row["chrom"])
            start = int(row["start"])
            end = int(row["end"])

            try:
                # 逐碱基查询
                values = self.bw.values(chrom, start, end)

                if values is None:
                    # 染色体不存在
                    densities.append(None)
                    continue

                # 过滤 NaN（NaN != NaN）
                valid_values = [v for v in values if v == v]

                if len(valid_values) == 0:
                    # 完全无覆盖 - 返回 NA
                    densities.append(None)
                else:
                    # 模式 B:  sum / covered_length
                    mean_val = sum(valid_values) / len(valid_values)
                    densities.append(round(mean_val, 4))

            except Exception:
                densities.append(None)

        return densities

    def close(self) -> None:
        """Close the BigWig file handle."""
        try:
            self.bw.close()
        except Exception:
            pass


# -----------------------
# Worker Functions (for BED processing)
# -----------------------
def compute_weighted_for_bed_value(
    bin_bed_path: str,
    value_bed_path: str,
    value_col: int,
    tmp_dir: str,
) -> str:
    """
    Compute weighted density for a BED value file (模式 B).

    BED 加权计算 Worker 函数：
    1. 使用 bedtools intersect 计算几何交集
    2. 计算信号贡献:  overlap × value
    3. 按 Bin 聚合:  sum(overlap × value) / sum(overlap)
    4. 无交集的 Bin 标记为 NA

    Args:
        bin_bed_path: Path to temporary bin BED file.
        value_bed_path: Path to value BED file.
        value_col:  0-based column index for value in value BED file.
        tmp_dir:  Temporary directory for intermediate files.

    Returns:
        Path to grouped result TSV file.
    """
    tmp_intersect:  Optional[str] = None
    tmp_grouped: Optional[str] = None

    try:
        # Perform intersection with overlap length (-wo)
        a = BedTool(bin_bed_path)
        b = BedTool(value_bed_path)
        res = a.intersect(b, wo=True)

        # Save intersection result to temp file
        with tempfile.NamedTemporaryFile(
            dir=tmp_dir, delete=False, mode="w", suffix=".bed"
        ) as fh:
            tmp_intersect = fh.name
        res.saveas(tmp_intersect)

        # Read bin file to get all bins (for NA filling)
        bin_df = pd.read_csv(
            bin_bed_path,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "bin_id", "region_uid"],
        )

        # Prepare output file
        with tempfile.NamedTemporaryFile(
            dir=tmp_dir, delete=False, mode="w", suffix=".tsv"
        ) as gh:
            tmp_grouped = gh.name

        # Handle empty intersection
        if os.path.getsize(tmp_intersect) == 0:
            # No overlaps - all bins get NA
            result_df = bin_df.copy()
            result_df["weighted_value"] = pd.NA
            result_df.to_csv(tmp_grouped, sep="\t", header=True, index=False)
            return tmp_grouped

        # Read intersection results
        df = pd.read_csv(tmp_intersect, sep="\t", header=None)

        # Column indices:
        # 0-4: bin columns (chrom, start, end, bin_id, region_uid)
        # 5+: value file columns
        # -1: overlap length

        if df.shape[1] <= 5:
            # Malformed result
            result_df = bin_df.copy()
            result_df["weighted_value"] = pd.NA
            result_df.to_csv(tmp_grouped, sep="\t", header=True, index=False)
            return tmp_grouped

        # Extract overlap length
        overlap = df.iloc[:, -1].astype(float)

        # Get value column (offset by 5 bin columns)
        val_col_idx = 5 + value_col
        if val_col_idx >= df.shape[1] - 1:  # -1 for overlap column
            result_df = bin_df.copy()
            result_df["weighted_value"] = pd.NA
            result_df.to_csv(tmp_grouped, sep="\t", header=True, index=False)
            return tmp_grouped

        values = pd.to_numeric(df.iloc[:, val_col_idx], errors="coerce").astype(float)

        # 模式 B：未覆盖区域不参与平均值计算
        # 计算公式：sum(overlap × value) / sum(overlap)
        # Create aggregation DataFrame
        agg_df = df.iloc[:, :5].copy()
        agg_df.columns = ["chrom", "start", "end", "bin_id", "region_uid"]
        agg_df["signal_area"] = (overlap * values).values  # 分子部分：overlap × value
        agg_df["covered_length"] = overlap.values  # 分母部分：overlap

        # Group by bin and sum both components
        grouped = agg_df.groupby(
            ["chrom", "start", "end", "bin_id", "region_uid"], as_index=False
        ).agg({
            "signal_area": "sum",
            "covered_length": "sum"
        })

        # Calculate weighted mean:  signal_area / covered_length
        grouped["weighted_value"] = np.where(
            grouped["covered_length"] > 0,
            (grouped["signal_area"] / grouped["covered_length"]).round(4),
            np.nan
        )

        # 只保留需要的列
        grouped = grouped[["chrom", "start", "end", "bin_id", "region_uid", "weighted_value"]]

        # Merge with all bins to fill NA for non-overlapping bins
        result_df = bin_df.merge(
            grouped,
            on=["chrom", "start", "end", "bin_id", "region_uid"],
            how="left",
        )

        result_df.to_csv(tmp_grouped, sep="\t", header=True, index=False)
        return tmp_grouped

    except Exception:
        logging.exception("Error in compute_weighted_for_bed_value worker.")
        # Return empty result on error
        if tmp_grouped is None:
            with tempfile.NamedTemporaryFile(
                dir=tmp_dir, delete=False, mode="w", suffix=".tsv"
            ) as gh:
                tmp_grouped = gh.name
                gh.write("chrom\tstart\tend\tbin_id\tregion_uid\tweighted_value\n")
        return tmp_grouped

    finally:
        _safe_remove(tmp_intersect)


# -----------------------
# Main Processing Flow
# -----------------------
def process_region_file(
    region_file: str,
    value_files: List[str],
    chrom_size_df: pd.DataFrame,
    value_col: int,
    extend_bp: Optional[int],
    extend_ratio: Optional[float],
    binnum:  int,
    retain_cols: List[int],
    outdir: str,
    threads: int,
    max_bins_per_batch: int,
    no_header: bool,
    save_extended_region: bool,
    tmp_dir: str,
) -> None:
    """
    Process a single region file end-to-end.

    主处理流程：
    1. 读取并预处理 Region 文件
    2. 按染色体分批处理
    3. 计算 BigWig 和 BED 的加权值
    4. 合并元数据并输出

    Args:
        region_file: Path to region BED file.
        value_files: List of value file paths (BED or BigWig).
        chrom_size_df:  Chromosome sizes DataFrame.
        value_col: 0-based column index for value in BED files.
        extend_bp: Fixed bp to extend (or None).
        extend_ratio: Ratio to extend (or None).
        binnum: Number of bins per region.
        retain_cols: List of column indices to retain.
        outdir: Output directory.
        threads: Number of worker threads.
        max_bins_per_batch: Maximum bins per batch.
        no_header: If True, don't write header to output.
        save_extended_region: If True, save extended region file.
        tmp_dir: Temporary directory.
    """
    #logging.info("Processing region file: %s", region_file)

    # Read region file
    try:
        region_df = read_bed_file_auto_sep(region_file)
        logging.info("Loaded region file '%s' shape=%s", region_file, region_df.shape)
    except Exception as exc:
        logging.error("Failed to read region file %s: %s", region_file, exc)
        return

    # Initialize processor
    rp = RegionProcessor(chrom_size_df)

    # Step 1: Filter and clip
    region_df = rp.filter_and_clip(region_df)
    logging.info("After filtering and clipping:  %d regions", region_df.shape[0])

    if region_df.empty:
        logging.warning("No valid regions after filtering.  Skipping.")
        return

    # Step 2: Extend regions (based on clipped length)
    if extend_bp is not None or extend_ratio is not None:
        region_df = rp.extend_regions(region_df, extend_bp, extend_ratio)
        if extend_bp is not None:
            logging.info("After extension (bp=%d): %d regions", extend_bp, region_df.shape[0])
        else:
            logging.info("After extension (ratio=%.2f): %d regions", extend_ratio, region_df.shape[0])

    region_df = region_df.reset_index(drop=True)

    # Step 3: Save extended region file if requested
    if save_extended_region:
        extended_out = region_df.iloc[:, :3].copy()
        extended_out.columns = ["chrom", "start", "end"]
        extended_out["region_uid"] = range(len(extended_out))

        # 添加 retain_cols
        for c in retain_cols:
            colname = f"region_col_{c}"
            if c < region_df.shape[1]:
                extended_out[colname] = region_df.iloc[:, c].values

        extended_path = os.path.join(outdir, extended_region_filename(region_file))
        extended_out.to_csv(extended_path, sep="\t", header=False, index=False)
        logging.info("Saved extended region file: %s", extended_path)

    # Step 4: Extract metadata
    region_metadata = rp.extract_metadata(region_df, retain_cols)

    # Prepare chromosome order
    chrom_order = chrom_size_df["chrom"].to_list()

    # Separate value files by type
    bw_files = [v for v in value_files if is_bigwig_file(v)]
    bed_files = [v for v in value_files if not is_bigwig_file(v)]

    logging.info("Value files: %d BigWig, %d BED", len(bw_files), len(bed_files))

    # Collect output chunks
    output_chunks: List[pd.DataFrame] = []

    # Process by chromosome
    for chrom in chrom_order:
        chrom_mask = region_df.iloc[:, 0] == chrom
        chrom_regions = region_df[chrom_mask]

        if chrom_regions.empty:
            continue

        logging.info("Processing %s:  %d regions", chrom, chrom_regions.shape[0])

        # Batch processing
        regions_per_batch = max(1, max_bins_per_batch // max(1, binnum))
        n_regions = chrom_regions.shape[0]
        start_idx = 0

        while start_idx < n_regions:
            end_idx = min(n_regions, start_idx + regions_per_batch)
            batch = chrom_regions.iloc[start_idx:end_idx].copy()

            # Preserve original global region indices
            orig_indices = batch.index.tolist()
            batch = batch.reset_index(drop=True)

            # Generate bins
            bin_df = BinCalculator.regions_df_to_bins_df(batch, binnum)

            if bin_df.empty:
                start_idx = end_idx
                continue

            # Map region_uid back to global indices
            bin_df["region_uid"] = bin_df["region_uid"].apply(
                lambda i: int(orig_indices[int(i)])
            )

            logging.debug(
                "Chrom %s batch [%d:%d) produced %d bins",
                chrom,
                start_idx,
                end_idx,
                bin_df.shape[0],
            )

            # Initialize result container
            batch_out = bin_df.copy().reset_index(drop=True)
            for val_file in value_files:
                batch_out[Path(val_file).stem] = pd.NA

            # Process BigWig files (in main process)
            for bw_file in bw_files:
                vname = Path(bw_file).stem
                try:
                    logging.debug(
                        "Querying BigWig %s for %d bins", bw_file, bin_df.shape[0]
                    )
                    runner = BigWigRunner(bw_file)
                    densities = runner.query_bins_density(
                        bin_df[["chrom", "start", "end"]]
                    )
                    runner.close()
                    batch_out[vname] = densities
                except Exception:
                    logging.exception("Failed to process BigWig %s; filling NA.", bw_file)
                    batch_out[vname] = pd.NA

            # Process BED files (in worker processes)
            if bed_files:
                # Write bin_df to temporary file
                bin_bed_path:  Optional[str] = None
                try:
                    with tempfile.NamedTemporaryFile(
                        dir=tmp_dir, delete=False, mode="w", suffix=".bed"
                    ) as rf:
                        bin_bed_path = rf.name
                    bin_df.to_csv(
                        bin_bed_path,
                        sep="\t",
                        header=False,
                        index=False,
                        columns=["chrom", "start", "end", "bin_id", "region_uid"],
                    )
                except Exception:
                    logging.exception("Failed to write temporary bin BED for workers.")
                    _safe_remove(bin_bed_path)
                    raise

                # Launch workers
                grouped_paths:  Dict[str, Optional[str]] = {}
                with ProcessPoolExecutor(
                    max_workers=max(1, min(len(bed_files), threads))
                ) as exe:
                    futures = {}
                    for bed_file in bed_files:
                        fut = exe.submit(
                            compute_weighted_for_bed_value,
                            bin_bed_path,
                            bed_file,
                            value_col,
                            tmp_dir,
                        )
                        futures[fut] = bed_file

                    for fut in futures:
                        bed_file = futures[fut]
                        try:
                            grouped_path = fut.result()
                            grouped_paths[bed_file] = grouped_path
                        except Exception:
                            logging.exception("Worker failed for %s", bed_file)
                            grouped_paths[bed_file] = None

                # Merge grouped results
                for bed_file in bed_files:
                    vname = Path(bed_file).stem
                    gp = grouped_paths.get(bed_file)

                    if not gp:
                        logging.warning(
                            "No grouped output for %s; leaving NA", bed_file
                        )
                        continue

                    try:
                        gdf = pd.read_csv(gp, sep="\t", header=0)
                        if not gdf.empty:
                            merged = batch_out.merge(
                                gdf[
                                    [
                                        "chrom",
                                        "start",
                                        "end",
                                        "bin_id",
                                        "region_uid",
                                        "weighted_value",
                                    ]
                                ],
                                on=["chrom", "start", "end", "bin_id", "region_uid"],
                                how="left",
                            )
                            batch_out[vname] = merged["weighted_value"].values
                    except Exception:
                        logging.exception("Failed to read/merge grouped file %s", gp)
                    finally:
                        _safe_remove(gp)

                # Cleanup
                _safe_remove(bin_bed_path)

            output_chunks.append(batch_out)
            start_idx = end_idx

    # Combine all batches
    if output_chunks:
        out_df = pd.concat(output_chunks, ignore_index=True)
    else:
        out_df = pd.DataFrame(
            columns=["chrom", "start", "end", "bin_id", "region_uid"]
        )
        for v in value_files:
            out_df[Path(v).stem] = pd.NA

    # Merge region metadata if any
    if region_metadata is not None:
        out_df = out_df.merge(region_metadata, on="region_uid", how="left")

    # Write output file
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, output_filename(region_file))

    try:
        out_df.to_csv(
            out_path,
            sep="\t",
            header=not no_header,
            index=False,
            na_rep="NA",
            float_format="%.4f",
        )
        logging.info("Wrote output:  %s", out_path)
    except Exception:
        logging.exception("Failed to write output file %s", out_path)


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Main entry point."""
    args = parse_arguments(argv)
    setup_logging(args.verbose)

    # Parse input files
    region_files = parse_multi_files(args.region)
    value_files = parse_multi_files(args.values)

    # Parse retain columns
    retain_cols:  List[int] = []
    if args.retain_cols:
        retain_cols = [int(x) for x in args.retain_cols.split(",") if x.strip()]
        if any(i < 3 for i in retain_cols):
            logging.error(
                "retain-cols should refer to extra columns (0-based index >= 3).Aborting."
            )
            sys.exit(1)

    # Read chromosome sizes
    try:
        chrom_size_df = read_chrom_size_file(args.chrom_size)
    except Exception as exc:
        logging.error("Failed to read chrom size file %s: %s", args.chrom_size, exc)
        sys.exit(1)

    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)

    # Process each region file
    with tempfile.TemporaryDirectory(prefix="region_binvalue_tmp_") as tmp_dir:
        logging.debug("Using tmp dir: %s", tmp_dir)

        for region_file in region_files:
            try:
                process_region_file(
                    region_file=region_file,
                    value_files=value_files,
                    chrom_size_df=chrom_size_df,
                    value_col=args.value_col,
                    extend_bp=args.extend_bp,
                    extend_ratio=args.extend_ratio,
                    binnum=args.binnum,
                    retain_cols=retain_cols,
                    outdir=args.outdir,
                    threads=args.threads,
                    max_bins_per_batch=args.max_bins_per_batch,
                    no_header=args.no_header,
                    save_extended_region=args.save_extended_region,
                    tmp_dir=tmp_dir,
                )
            except Exception:
                logging.exception("Unexpected error processing %s", region_file)


if __name__ == "__main__":
    main()