from Bio import Entrez
import time

Entrez.email = "***"


def get_transcript_regulator(chrm: str, start: int, end: int, window: int = 150):
    # 検索範囲を設定
    search_range_start = start - window
    search_range_end = end + window
    chromosome = "".join(filter(str.isdigit, chrm))
    search_term = f"""({search_range_start}[Base Position] : 
        {search_range_end}[Base Position]) AND {chromosome}[Chromosome] 
        AND Homo sapiens[Organism]
    """
    time.sleep(2)
    # NCBI Geneデータベースで検索
    with Entrez.esearch(db="gene", term=search_term) as handle:
        record = Entrez.read(handle)

    if record is None:
        print("No genes found in the specified range.")
        return

    if record["IdList"]:
        gene_ids = record["IdList"]
        print("Found gene IDs:", gene_ids)

        result = {}
        for gene_id in gene_ids:
            time.sleep(2)
            with Entrez.efetch(db="gene", id=gene_id, retmode="xml") as handle:
                gene_record = Entrez.read(handle)

                # 遺伝子名を表示
                gene_name = gene_record[0]["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
                # TSS
                tss = gene_record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"][
                    "Seq-interval_from"
                ]

                print(gene_record[0]["Entrezgene_locus"][0]["Gene-commentary_accession"])
                print(f"Gene ID: {gene_id}, Gene Name: {gene_name}")
                result[gene_id] = {"name": gene_name, "tss": tss}
        return result


if __name__ == "__main__":
    peaks = [(47821818, 47821822, 2.5)]  # 例: [(start, end, value)]

    for peak in peaks:
        get_transcript_regulator(peak[0], peak[1])
