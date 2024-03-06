import pyBigWig
import pandas as pd


class PyChIPProm:
    def __init__(self, bw_path: str):
        self.bw_path = bw_path
        self.bw = self._open_bw(bw_path)

    def _open_bw(self, path: str):
        # pyBigWigライブラリを使用してBigWigファイルを開く
        return pyBigWig.open(path)

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

    def query_ncbi_for_gene_id(self, gene_name):
        # Entrezを使用してNCBIからgene_idを取得
        # 省略: Entrezへのクエリと結果の処理
        pass

    def get_transcript_regulator(self, chrom: str, start: int, end: int):
        pass

    def output_results_to_csv(self, results, output_path):
        # pandasを使用して結果をcsvファイルに出力
        df = pd.DataFrame(results)
        df.to_csv(output_path, index=False)

    def analyze(self, chrom, start, end, threshold, output_path):
        # 解析の全プロセスを実行
        peaks = self.filter_peaks(chrom, start, end, threshold)
        results = []
        for peak in peaks:
            gene_id = self.query_ncbi_for_gene_id(peak["gene_name"])
            peak["gene_id"] = gene_id
            results.append(peak)
        self.output_results_to_csv(results, output_path)
