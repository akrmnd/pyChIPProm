import pyChIPProm
import argparse


def main(bw_path: str, output_path: str, threshold: float):
    chiprom = pyChIPProm.PyChIPProm("***", "***", bw_path)
    for chrom in chiprom.bw.chroms().keys():
        if chrom != "chr11":
            continue
        res = chiprom.analyze(chrom, 0, chiprom.bw.chroms(chrom), threshold)
        print(res)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to analyze BigWig files with pyChIPProm.")
    parser.add_argument("-i", "--bw_path", type=str, help="Path to the BigWig file")
    parser.add_argument("-o", "--output_path", type=str, help="Path to the output file where results will be saved")
    parser.add_argument("-t", "--threshold", type=float, help="Threshold value for filtering peaks")

    args = parser.parse_args()
    main(args.bw_path, args.output_path, args.threshold)
