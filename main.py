import pyChIPProm
import argparse
from pathlib import Path


def main(bw_path: str, gff_path: str, output_dir: str, threshold: float, output_type: str):
    chiprom = pyChIPProm.PyChIPProm(bw_path, gff_path)
    for chrom in chiprom.bw.chroms().keys():
        res = chiprom.analyze(chrom, 0, chiprom.bw.chroms(chrom), threshold)
        output_path = Path(output_dir) / f"{chrom}.{output_type}"
        chiprom.output(res, str(output_path), output_type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to analyze BigWig files with pyChIPProm.")
    parser.add_argument("--bw_path", type=str, help="Path to the BigWig file")
    parser.add_argument("--gff_path", type=str, help="Path to the GFF file")
    parser.add_argument("-o", "--output_dir", type=str, help="Path to the output dir where results will be saved")
    parser.add_argument("-t", "--threshold", type=float, help="Threshold value for filtering peaks")
    parser.add_argument("--output_type", type=str, default="csv", help="Output file format. default=csv")

    args = parser.parse_args()
    main(args.bw_path, args.gff_path, args.output_dir, args.threshold, args.output_type)
