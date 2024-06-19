# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Finding Nested Introns

# +
from collections import defaultdict
from pathlib import Path

import pandas as pd
from pybedtools import BedTool, Interval
from tqdm import tqdm

# -

# We define functions for loading different file types with SJ / intron data:

# +
SJ_FILE_COLUMN_NAMES = ['chromosome', 'start', 'end', 'strand', 'intron_motif',
                        'annotated', 'reads_unique', 'reads_multimapped', 'max_overhang']
BED_FILE_COLUMN_NAMES = ['chromosome', 'start', 'end', 'name', 'score', 'strand']


def load_sj_dataframe(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path,
                     delimiter='\t',
                     names=SJ_FILE_COLUMN_NAMES)
    # STAR SJ files are 1-based and interval end is inclusive, while .bed files are 0-based and interval end is exclusive.
    # We convert SJ files to .bed-style indexing (note that interval end is unchanged as 0-/1- indexing cancels out with exclusiveness / inclusiveness)
    df['start'] -= 1

    # We omit SJ with unspecified strand    
    df = df[df['strand'] != 0]
    df['strand'] = df['strand'].apply(lambda x: '+' if x == 1 else '-')

    # We keep track how many samples support this SJ # this will be useful when aggregating dataframes from multiple samples
    df['present_in_samples'] = 1
    return df


def load_bed_file_as_dataframe(file_path: Path) -> pd.DataFrame:
    return pd.read_csv(file_path, delimiter='\t', names=BED_FILE_COLUMN_NAMES)


# -

# We load SJ files, filter them by some quality control, and aggregate them to a single dataframe:

# +
sj_folder_path = Path('./splice_junction_files')
sj_dataframes = [load_sj_dataframe(path) for path in tqdm(sj_folder_path.iterdir(), desc='Loading SJ files')]


def quality_control_sj(df: pd.DataFrame, min_max_overhang=20, min_unique_reads=5) -> pd.DataFrame:
    df = df[df['max_overhang'] >= min_max_overhang]
    df = df[df['reads_unique'] >= min_unique_reads]
    return df


sj_dataframes = [quality_control_sj(sj) for sj in sj_dataframes]

sj_aggregated = pd.concat(sj_dataframes).groupby(['chromosome', 'start', 'end',
                                                  'strand', 'intron_motif', 'annotated']).agg(
    {'reads_unique': 'sum',
     'reads_multimapped': 'sum',
     'max_overhang': 'max',
     'present_in_samples': 'sum'}).reset_index()

sj_aggregated['present_in_samples_fraction'] = sj_aggregated['present_in_samples'] / len(sj_dataframes)
sj_aggregated
# -

# We may further e.g. filter only SJ present in all samples:

sj_aggregated = sj_aggregated[sj_aggregated['present_in_samples_fraction'] == 1]
sj_aggregated


# We define function for finding nested introns and apply it to dataframe with SJ:

def find_nested_introns(df_1: pd.DataFrame, df_2: pd.DataFrame) -> pd.DataFrame:
    """
    Find nested introns (SJ). The potential 'outside' and 'nested' introns are specified by different arguments (
      df_1 resp. df_2), which may be useful if we wish to apply different pre-filtering logic on them.
    :param df_1: Dataframe with potential 'outside' introns.
    :param df_2: Dataframe with potential 'nested' introns.
    :return: Dataframe same as df_1, with added columns 'has_nested_intron' and 'nested_intron_indices'. The column
    'has_nested_intron' contain boolean flag specifying whether this intron contains any nested introns. The column
    'nested_intron_indices' contain indices (from df_2) of nested introns.
    """
    df_output = df_1.copy()

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

    intersection_bedtool = bedtool_1.intersect(bedtool_2, sorted=True, s=True, wa=True, wb=True)
    intersections_by_query_id = defaultdict(list)
    for intersection_interval in intersection_bedtool:
        intersections_by_query_id[int(intersection_interval.fields[3])].append(int(intersection_interval.fields[9]))

    has_nested_intron: list[bool] = []
    nested_intron_indices: list[list[int]] = []

    for index, row in tqdm(df_1.iterrows(), desc='Finding nested introns', total=len(df_output)):
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

    df_output['has_nested_intron'] = has_nested_intron
    df_output['nested_intron_indices'] = nested_intron_indices
    return df_output


sj_aggregated = find_nested_introns(df_1=sj_aggregated, df_2=sj_aggregated)

# Some SJ have nested SJ inside them:

sj_with_nested_introns = sj_aggregated[sj_aggregated['has_nested_intron']]
sj_with_nested_introns

# The column 'nested_intron_indices' keeps indices (of second dataframe passed to find_nested_introns() function) which we may use to retrieve the nested SJ:

example_sj = sj_with_nested_introns.iloc[3]
example_nested_introns = sj_aggregated.loc[example_sj['nested_intron_indices']]
print(f"Example SJ: \n{example_sj}")
print(50 * '-')
print("Example nested introns of the SJ above:")
example_nested_introns

# We may find nested SJ inside different reference introns, if we wish to apply different filtering logic on the nested and 'outside' introns.
# We may e.g. find SJ that are nested inside introns from the reference genome:

reference_introns = load_bed_file_as_dataframe('./reference_genome/introns.bed')
reference_introns = find_nested_introns(reference_introns, sj_aggregated)

reference_introns_with_nested_sj = reference_introns[reference_introns['has_nested_intron']]
reference_introns_with_nested_sj

example_intron = reference_introns_with_nested_sj.iloc[3]
example_nested_sj = sj_aggregated.loc[example_intron['nested_intron_indices']]
print(f"Example intron with nested SJ: \n{example_intron}")
print(50 * '-')
print("Nested SJ of the intron above:")
example_nested_sj
