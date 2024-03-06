from pathlib import Path
from typing import Optional
import pyBigWig
from Bio import Entrez
import pandas as pd
from tqdm import tqdm


class PyChIPProm:
    def __init__(self, email: str, api_key: str, bw_path: Optional[str], window: int = 150):
        self.bw_path = bw_path
        self.bw = self.load_bw(bw_path)
        self.window = window
        self.EnterzGS = EntrezGeneSearch(email, api_key)
        self.results = {}

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

    def output_results_to_csv(self, results, output_path):
        # pandasを使用して結果をcsvファイルに出力
        df = pd.DataFrame(results)
        df.to_csv(output_path, index=False)

    def analyze(self, chrom: str, start: int, end: int, threshold: float):
        peaks = self.filter_peaks(chrom, start, end, threshold)
        res = []
        for peak_start, peak_end, peak_value in tqdm(peaks):
            gene_ids = self.EnterzGS.gene_search(
                "".join(filter(str.isdigit, chrom)), peak_start - self.window, peak_end + self.window
            )
            for gene_id in gene_ids:
                gene_detail = self.EnterzGS.fetch_gene_detail(gene_id, peak_start - self.window, peak_end + self.window)
                result = {
                    "peak_start": peak_start,
                    "peak_end": peak_end,
                    "peak_value": peak_value,
                    "gene_id": gene_id,
                    "gene_name": gene_detail.get("name"),
                    "tss": gene_detail.get("tss"),
                }
                res.append(result)
        self.results["chrom"] = res

        return res


class EntrezGeneSearch:
    def __init__(self, email: str, api_key: str):
        Entrez.email = email
        Entrez.api_key = api_key

    def gene_search(self, chrm: str, start: int, end: int, organism: str = "Homo sapiens"):
        search_term = f"({start}[Base Position]:{end}[Base Position]) AND {chrm}[Chromosome] AND {organism}[Organism]"
        with Entrez.esearch(db="gene", term=search_term) as handle:
            record = Entrez.read(handle)
        return record.get("IdList", [])

    def fetch_gene_detail(self, gene_id: str, peak_start: int, peak_end: int):
        try:
            with Entrez.efetch(db="gene", id=gene_id, retmode="xml") as handle:
                gene_record = Entrez.read(handle)
        except Exception as e:
            print(f"Error fetching gene detail for gene ID {gene_id}: {e}")
            return {}

        try:
            gene_name = gene_record[0].get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_locus")
            tss = (
                gene_record[0]
                .get("Entrezgene_locus", [{}])[0]
                .get("Gene-commentary_seqs", [{}])[0]
                .get("Seq-loc_int", {})
                .get("Seq-interval", {})
                .get("Seq-interval_from")
            )
            if gene_name and tss and peak_start <= int(tss) <= peak_end:
                return {"name": gene_name, "tss": tss}
            else:
                print(f"TSS for gene ID {gene_id} is out of the specified peak range or not found.")
                return {}
        except (IndexError, KeyError) as e:
            print(f"Error parsing gene detail for gene ID {gene_id}: {e}")
            return {}
