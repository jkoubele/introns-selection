import pandas as pd
from pathlib import Path


# sj_column_names = ['chromosome', 'start', 'end', 'strand', 'intron_motif', 'annotated', 'reads_unique', 'reads_multimapped', 'max_overhang']

paths =  ['/cellfile/datapublic/jkoubele/FLI_total_RNA/SJ/no003-1_OA3.SJ.out.tab',
          '/cellfile/datapublic/jkoubele/FLI_total_RNA/SJ/no004-0_OD1.SJ.out.tab']


def load_sj_dataframe(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, 
                     delimiter='\t', 
                     names=['chromosome', 'start', 'end', 'strand', 'intron_motif', 
                            'annotated', 'reads_unique', 'reads_multimapped', 'max_overhang'])
    # STAR SJ files are 1-based and interval end is inclusive, while .bed files are 0-based and interval end is exclusive.
    # We convert SJ files to .bed-style indexing (note that interval end is unchanged as 0-/1- indexing cancels out with exclusiveness / inclusiveness)
    df['start'] -= 1
    
    # We keep track how many samples support this SJ
    df['present_in_samples'] = 1
    return df
    

def quality_control_sj(df: pd.DataFrame, min_max_overhang = 20, min_unique_reads = 5) -> pd.DataFrame:
    df = df[df['max_overhang']>=min_max_overhang]
    df = df[df['reads_unique']>=min_unique_reads]
    
    # We omit SJ with undefined strand
    df = df[df['strand'] != 0]
    return df


sj_all = [load_sj_dataframe(path) for path in paths]

df = pd.concat([quality_control_sj(sj) for sj in sj_all])

df_aggregated = df.groupby(['chromosome', 'start', 'end', 'strand', 'intron_motif', 'annotated']).agg({'reads_unique': 'sum', 
                                                                                                            'reads_multimapped': 'sum',
                                                                                                            'max_overhang': 'max',
                                                                                                            'present_in_samples': 'sum'})

df_aggregated['present_in_samples_fraction'] = df_aggregated['present_in_samples'] / len(sj_all)

gtf_mouse = pd.read_csv('/cellfile/datapublic/jkoubele/STAR_2.7.11a/reference_genomes/GRCm38/gencode.vM10.primary_assembly.annotation.gtf' ,
                  delimiter='\t', header=5,
                  names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

gtf_drosophila = pd.read_csv("/cellfile/datapublic/rmarel_1/Internship/Pol_II_RvdM/data/gtf_reference_files/Drosophila_melanogaster.BDGP6.90.gtf" ,
                  delimiter='\t', header=5,
                  names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
