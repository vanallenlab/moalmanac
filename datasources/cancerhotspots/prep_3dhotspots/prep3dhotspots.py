import pandas as pd
df2_cols = ['Gene', 'Residue', 'p-value', 'Class']
df5_cols = ['Gene', 'Reference amino acid', 'Variant amino acid', 'Amino_Acid_Position']
df2 = pd.read_csv('3d_hotspots_T2.txt', sep='\t', usecols = df2_cols)
df5 = pd.read_csv('3d_hotspots_T5.txt', sep = '\t', usecols = df5_cols)

df5['Residue'] = df5.loc[:,'Reference amino acid'] + df5.loc[:,'Amino_Acid_Position'].astype(str)
df5['alteration'] = 'p.' + df5.loc[:,'Reference amino acid'] + df5.loc[:,'Amino_Acid_Position'].astype(str) + df5.loc[:,'Variant amino acid']
df5 = df5.drop(['Reference amino acid', 'Variant amino acid', 'Amino_Acid_Position'], axis = 1)

df = pd.merge(df2, df5, on=['Gene', 'Residue'], how='left')

class_map = {
    'Cluster-exclusive': 1,
    'Hotspot-linked': 2,
    'Hotspot': 3
}

df['cancerhotspots3D_bin'] = df['Class'].map(class_map)

df.to_csv('hotspots3d.txt', sep = '\t', index = False)
