import pandas as pd


df_intron = pd.read_csv('/cellfile/datapublic/jkoubele/introns-selection/data/introns_after_merge.bed',
                 delimiter='\t',
                 names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])

df_bash_script = df = pd.read_csv('/cellfile/datapublic/jkoubele/introns-selection/data/gencode.vM10.primary_assembly.annotation.intron.bed',
                                  delimiter='\t',
                                  names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])

for col_name in ['chrom', 'chromStart', 'chromEnd', 'strand']:
    assert all(df_intron[col_name] == df_bash_script[col_name])