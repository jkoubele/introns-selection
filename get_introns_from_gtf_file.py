from pathlib import Path

import pandas as pd
from pybedtools import BedTool, Interval
from tqdm import tqdm

if __name__ == "__main__":
    data_folder = Path("./reference_genome")
    input_file_path = data_folder / "gencode.vM10.primary_assembly.annotation.gtf"

    # Fast check of the file with reference genome:
    gtf_df = pd.read_csv(input_file_path,
                         header=4,
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                                'attribute'],
                         delimiter="\t")

    print(f"Num. genes: {sum(gtf_df['feature'] == 'gene')}")
    print(f"Possible features: {set(gtf_df['feature'])}")

    # Processing of the reference genome to obtain introns
    gtf_records = BedTool(input_file_path)
    genes = BedTool(
        [record for record in tqdm(gtf_records, desc='Extracting genes') if record.fields[2] == 'gene']).sort()
    exons = BedTool(
        [record for record in tqdm(gtf_records, desc='Extracting exons') if record.fields[2] == 'exon']).sort()
    utr = BedTool(
        [record for record in tqdm(gtf_records, desc="Extracting UTR") if record.fields[2] == 'UTR']).sort()
    introns = genes.subtract(exons, s=True).subtract(utr, s=True).sort()
    introns = BedTool([Interval(chrom=x.chrom, start=x.start, end=x.end,
                                name=x.fields[-1].split()[1][1:-2], score='.',
                                strand=x.strand)
                       for x in tqdm(introns)])

    introns = introns.merge(s=True, c=[4, 5, 6], o='distinct').sort()
    introns.saveas(data_folder / 'introns.bed')
