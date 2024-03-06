import pyChIPProm
import argparse


def main(bw_path: str, output_path: str, threshold: float):
    chiprom = pyChIPProm.PyChIPProm(bw_path)

    for chrom in chiprom.bw.chroms().keys():
        if chrom != "chr11":
            continue
        peaks = chiprom.filter_peaks(chrom, 0, chiprom.bw.chroms(chrom), threshold)
        for peak in peaks:
            result = chiprom.get_transcript_regulator(chrom, peak[0], peak[1])
            # [...{"gene_id": int, "gene_name": str, "tss": Optional(int)}]
            if result is not None:
                print(f"chrom: {chrom} start: {peak[0]} val: {peak[2]} result: {result}")
        break


if __name__ == "__main__":
    main()
