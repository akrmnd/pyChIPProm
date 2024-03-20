import json
from pathlib import Path
from typing import Any, Optional
import pyBigWig
import pandas as pd
from tqdm import tqdm


class PyChIPProm:
    def __init__(self, bw_path: str, gff_path: str, window: int = 150):
        self.bw = self.load_bw(bw_path)
        self.df = self.load_gff(gff_path)
        self.window = window

    def load_bw(self, path: Optional[str]):
        if path is None:
            raise ValueError("BigWig file path must not be None")
        bw_path = Path(path)

        # ファイルの存在とファイルタイプをチェック
        if not bw_path.exists():
            raise FileNotFoundError(f"The specified BigWig file was not found: {path}")
        if not bw_path.is_file():
            raise ValueError(f"The specified path is not a file: {path}")

        # pyBigWigライブラリを使用してBigWigファイルを開く
        try:
            bw_file = pyBigWig.open(str(bw_path))
        except Exception as e:
            raise IOError(f"Failed to open BigWig file {path}: {e}")

        return bw_file

    def load_gff(self, path: Optional[str]) -> pd.DataFrame:
        """GFF3ファイルを読み込み、pandas DataFrameに変換する"""
        if path is None:
            raise ValueError("GFF file path must not be None")
        col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        df = pd.read_csv(path, sep="\t", comment="#", names=col_names)
        # 属性フィールドをパースして新しい列を作成する（例：ID, Nameなど）
        # NCBIのseqidをchr形式に変換するマッピング
        seqid_mapping = {
            "NC_000001.11": "chr1",
            "NC_000002.12": "chr2",
            "NC_000003.12": "chr3",
            "NC_000004.12": "chr4",
            "NC_000005.10": "chr5",
            "NC_000006.12": "chr6",
            "NC_000007.14": "chr7",
            "NC_000008.11": "chr8",
            "NC_000009.12": "chr9",
            "NC_000010.11": "chr10",
            "NC_000011.10": "chr11",
            "NC_000012.12": "chr12",
            "NC_000013.11": "chr13",
            "NC_000014.9": "chr14",
            "NC_000015.10": "chr15",
            "NC_000016.10": "chr16",
            "NC_000017.11": "chr17",
            "NC_000018.10": "chr18",
            "NC_000019.10": "chr19",
            "NC_000020.11": "chr20",
            "NC_000021.9": "chr21",
            "NC_000022.11": "chr22",
            "NC_000023.11": "chrX",
            "NC_000024.10": "chrY",
            "NC_012920.1": "chrMT",
        }
        # seqid列をマッピングに基づいて変換
        df["seqid"] = df["seqid"].map(seqid_mapping)
        df["gene_id"] = df["attributes"].str.extract(r"ID=([^;]+)")
        df["gene_name"] = df["attributes"].str.extract(r"gene=([^;]+)")
        df["gene_type"] = df["attributes"].str.extract(r"gene_biotype=([^;]+)")
        # 不要な列を削除
        df.drop("attributes", axis=1, inplace=True)
        return df[(df["gene_type"] == "protein_coding") & (df["type"] == "gene")]

    def filter_peaks(self, chrom: str, start: int, end: int, threshold: float):
        """
        BigWigファイルから指定したchromosomeに関して、startとendの間のcountのピークを取得します。
        取得するピークはthresholdを超える連続するbin内で最大の位置です。
        args:
            data: BigWig file
            chrom: string chromosome number
            start: int start index
            end: int end index
            threshold: 取得したいピークの閾値
        return:
            [...(start, end, value)]: 抽出された各ピークを格納したリスト
        """

        df = self._get_peak_data(chrom, start, end)
        # 特定の値を超える行のみをフィルタリングし、コピーを作成してWarningを回避
        filtered_df = self._filter_peaks_above_threshold(df, threshold)
        # グループ化のための新しい列を安全に追加
        grouped_peaks = self._group_peaks(filtered_df)
        # 各グループで最大valueを持つ行を取得
        grouped_peaks = self._group_peaks(filtered_df)
        # 結果をタプルのリストとして返す
        max_peaks = self._select_max_peak_per_group(grouped_peaks)

        return max_peaks

    def _get_peak_data(self, chrom, start, end):
        return pd.DataFrame(self.bw.intervals(chrom, start, end), columns=["start", "end", "value"])

    def _filter_peaks_above_threshold(self, df, threshold):
        return df[df["value"] > threshold].copy()

    def _group_peaks(self, filtered_df):
        filtered_df.loc[:, "group"] = (filtered_df["start"] > filtered_df["end"].shift()).cumsum()
        return filtered_df

    def _select_max_peak_per_group(self, grouped_df):
        max_peaks = grouped_df.groupby("group").apply(lambda x: x.loc[x["value"].idxmax()])
        return [(int(row.start), int(row.end), float(row.value)) for index, row in max_peaks.iterrows()]

    def get_transcript_regulator(self, chrom: str, start: int, end: int):
        pass

    def output(self, results: list[dict[str, Any]], output_path: str, output_type: str = "csv"):
        df = pd.DataFrame(results)
        if output_type == "csv":
            # pandasを使用して結果をcsvファイルに出力
            df.to_csv(output_path, index=False)
        elif output_type == "json":
            with open(output_path, "w") as file:
                json.dump({"result": results}, file)
        else:
            raise ValueError("Output type error. must be 'csv' or 'json'.")

    def analyze(self, chrom: str, start: int, end: int, threshold: float):
        peaks = self.filter_peaks(chrom, start, end, threshold)
        res = []
        genes_chrom = self.df[(self.df["seqid"] == chrom)]
        for peak_start, peak_end, peak_value in tqdm(peaks):
            genes_target = genes_chrom[
                ~((genes_chrom["end"] < (peak_start - self.window)) | (genes_chrom["start"] > (peak_end + self.window)))
            ]

            if len(genes_target) == 0:
                res.append(
                    {
                        "peak_start": peak_start,
                        "peak_end": peak_end,
                        "peak_value": peak_value,
                        "gene_id": None,
                        "gene_name": None,
                        "gene_start": None,
                        "tss": False,
                    }
                )
                continue

            for _, gene_t in genes_target.iterrows():
                if (gene_t["start"] >= peak_start - self.window) and (gene_t["start"] <= peak_end + self.window):
                    res.append(
                        {
                            "peak_start": peak_start,
                            "peak_end": peak_end,
                            "peak_value": peak_value,
                            "gene_id": gene_t.get("gene_id"),
                            "gene_name": gene_t.get("gene_name"),
                            "gene_start": gene_t.get("start"),
                            "tss": True,
                        }
                    )
                else:
                    res.append(
                        {
                            "peak_start": peak_start,
                            "peak_end": peak_end,
                            "peak_value": peak_value,
                            "gene_id": gene_t.get("gene_id"),
                            "gene_name": gene_t.get("gene_name"),
                            "gene_start": gene_t.get("start"),
                            "tss": False,
                        }
                    )

        return res
