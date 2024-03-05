import pandas as pd


def filter_peaks(data, chrom: str, start: int, end: int, threshold: float) -> list[tuple[int, int, float]]:
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

    df = pd.DataFrame(data.intervals(chrom, start, end), columns=["start", "end", "value"])
    # 特定の値を超える行のみをフィルタリングし、コピーを作成してWarningを回避
    filtered_df = df[df["value"] > threshold].copy()
    # グループ化のための新しい列を安全に追加
    filtered_df.loc[:, "group"] = (filtered_df["start"] > filtered_df["end"].shift()).cumsum()
    # 各グループで最大valueを持つ行を取得
    grouped_max_value = filtered_df.groupby("group").apply(lambda x: x.loc[x["value"].idxmax()])
    # 不要になった"group"列を削除
    grouped_max_value = grouped_max_value.drop(columns=["group"])
    # 結果をタプルのリストとして返す
    return [(int(start), int(end), float(value)) for start, end, value in grouped_max_value.values]
