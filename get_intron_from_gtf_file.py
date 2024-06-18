from pybedtools import BedTool
from tqdm import tqdm
from collections import Counter
from pathlib import Path
import pandas as pd

data_folder = Path("/cellfile/datapublic/jkoubele/introns-selection/data")
input_file_path = data_folder / "gencode.vM10.primary_assembly.annotation.gtf"

gtf_df = pd.read_csv(input_file_path,
                     header=4,
                     names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
                     delimiter="\t")

print(f"Num. genes: {sum(gtf_df['feature'] == 'gene')}")
print(f"Possible features: {set(gtf_df['feature'])}")

# %%
gtf_records = BedTool(input_file_path)

genes = BedTool([record for record in tqdm(gtf_records, desc='Extracting genes') if record.fields[2] == 'gene']).sort()
exons = BedTool([record for record in tqdm(gtf_records, desc='Extracting exons') if record.fields[2] == 'exon']).sort()
utr = BedTool(
    [record for record in tqdm(gtf_records, desc="Extracting UTR") if record.fields[2] == 'UTR']).sort()
introns_before_merge = genes.subtract(exons, s=True).subtract(utr, s=True).sort()

introns_after_merge = introns_before_merge.merge(s=True) #, c=[1], o='distinct')
introns_after_merge = introns_before_merge.merge(s=True, c=[3,6,7],o='distinct')

introns_after_merge.saveas(data_folder /'introns_after_merge.bed')

introns = genes.subtract(exons, s=True).subtract(utr, s=True).sort().merge(s=True, c='4,5,6', o='distinct')

#%%
introns.saveas(data_folder /'introns.bed')

# %%
bash_script_genes = BedTool(data_folder / 'gencode.vM10.primary_assembly.annotation.gene.bed')
bash_script_exons = BedTool(data_folder / 'gencode.vM10.primary_assembly.annotation.exon.bed')
bash_script_utr = BedTool(data_folder / 'gencode.vM10.primary_assembly.annotation.utr.bed')
bash_script_introns = BedTool(data_folder / 'gencode.vM10.primary_assembly.annotation.intron.bed')

# %%
assert len(genes) == len(bash_script_genes)
assert len(exons) == len(bash_script_exons)
assert len(utr) == len(bash_script_utr)
assert len(introns) == len(bash_script_introns)
