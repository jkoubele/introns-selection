from collections import defaultdict
from pathlib import Path

import pandas as pd
from pybedtools import BedTool, Interval
from tqdm import tqdm

# sj_column_names = ['chromosome', 'start', 'end', 'strand', 'intron_motif', 'annotated', 'reads_unique', 'reads_multimapped', 'max_overhang']

paths = ['/cellfile/datapublic/jkoubele/FLI_total_RNA/SJ/no003-1_OA3.SJ.out.tab',
         '/cellfile/datapublic/jkoubele/FLI_total_RNA/SJ/no004-0_OD1.SJ.out.tab']


def load_sj_dataframe(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path,
                     delimiter='\t',
                     names=['chromosome', 'start', 'end', 'strand', 'intron_motif',
                            'annotated', 'reads_unique', 'reads_multimapped', 'max_overhang'])
    # STAR SJ files are 1-based and interval end is inclusive, while .bed files are 0-based and interval end is exclusive.
    # We convert SJ files to .bed-style indexing (note that interval end is unchanged as 0-/1- indexing cancels out with exclusiveness / inclusiveness)
    df['start'] -= 1

    # We omit SJ with unspecified strand    
    df = df[df['strand'] != 0]
    # 
    df['strand'] = df['strand'].apply(lambda x: '+' if x == 1 else '-')

    # We keep track how many samples support this SJ - this will be useful when aggregating dataframes from multiple samples
    df['present_in_samples'] = 1
    return df


def quality_control_sj(df: pd.DataFrame, min_max_overhang=20, min_unique_reads=5) -> pd.DataFrame:
    df = df[df['max_overhang'] >= min_max_overhang]
    df = df[df['reads_unique'] >= min_unique_reads]
    return df


sj_all = [load_sj_dataframe(path) for path in paths]
sj_all = [quality_control_sj(sj) for sj in sj_all]

df_aggregated = pd.concat(sj_all).groupby(['chromosome', 'start', 'end', 'strand', 'intron_motif', 'annotated']).agg(
    {'reads_unique': 'sum',
     'reads_multimapped': 'sum',
     'max_overhang': 'max',
     'present_in_samples': 'sum'}).reset_index()

df_aggregated['present_in_samples_fraction'] = df_aggregated['present_in_samples'] / len(sj_all)

# %%
df_1 = df_aggregated
df_2 = df_aggregated

bedtool_1 = BedTool([Interval(chrom=row['chromosome'],
                              start=row['start'],
                              end=row['end'],
                              name=index,
                              strand=row['strand']) for index, row in df_1.iterrows()]).sort()

bedtool_2 = BedTool([Interval(chrom=row['chromosome'],
                              start=row['start'],
                              end=row['end'],
                              name=index,
                              strand=row['strand']) for index, row in df_2.iterrows()]).sort()

intersection = bedtool_1.intersect(bedtool_2, sorted=True, s=True, wa=True, loj=True)
intersections_by_query_id = defaultdict(list)
for x in intersection:
    intersections_by_query_id[int(x.fields[3])].append(int(x.fields[9]))

# %%
has_nested_intron: list[bool] = []
nested_intron_indices: list[list[int]] = []

for index, row in tqdm(df_1.iterrows(), desc='Finding nested introns', total=len(df_1)):
    intersections_indices = intersections_by_query_id.get(index)
    if not intersections_indices:
        has_nested_intron.append(False)
        nested_intron_indices.append([])
        continue
    intersections_df = df_2.loc[intersections_indices]
    nested_introns_df = intersections_df[intersections_df['start'] > row['start']]
    nested_introns_df = nested_introns_df[nested_introns_df['end'] < row['end']]

    nested_introns_indices = list(nested_introns_df.index)
    has_nested_intron.append(bool(nested_introns_indices))
    nested_intron_indices.append(nested_introns_indices)

df_1['has_nested_intron'] = has_nested_intron
df_1['nested_intron_indices'] = nested_intron_indices

# %%
df_introns_with_nested_introns = df_1[df_1['has_nested_intron']]

example_intron = df_introns_with_nested_introns.iloc[1]
nested_introns = df_2.loc[example_intron.nested_intron_indices]

# gtf_mouse = pd.read_csv('/cellfile/datapublic/jkoubele/STAR_2.7.11a/reference_genomes/GRCm38/gencode.vM10.primary_assembly.annotation.gtf' ,
#                   delimiter='\t', header=5,
#                   names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

# gtf_drosophila = pd.read_csv("/cellfile/datapublic/rmarel_1/Internship/Pol_II_RvdM/data/gtf_reference_files/Drosophila_melanogaster.BDGP6.90.gtf" ,
#                   delimiter='\t', header=5,
#                   names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
