from pyChIPProm.get_transcript_regulator import get_transcript_regulator
from pyChIPProm.open import open
from pyChIPProm.filter_peaks import filter_peaks
from pprint import pprint


def main():
    # とりあえずのデータ
    bw = open("Sample.bw")

    for chrom in bw.chroms().keys():
        if chrom != "chr11":
            continue
        peaks = filter_peaks(bw, chrom, 0, bw.chroms(chrom), 15.0)

        pprint(f"{chrom}: {peaks}")
        for peak in peaks:
            result = get_transcript_regulator(chrom, peak[0], peak[1])
            if result is not None:
                print(f"chrom: {chrom} start: {peak[0]} val: {peak[2]} result: {result}")
        break


if __name__ == "__main__":
    main()
