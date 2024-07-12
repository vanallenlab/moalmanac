import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--exac', help='Tab delimited ExAC', required=True)
args = parser.parse_args()

exac = pd.read_csv(args.exac, sep='\t', low_memory=False)

cols = ['ALT', 'AF', 'AC', 'AC_AFR', 'AC_AMR', 'AC_EAS', 'AC_FIN', 'AC_NFE', 'AC_OTH', 'AC_SAS']
fillcols = list(set(exac.columns.tolist()) - set(cols))

idx = exac['ALT'].astype(str).str.contains(',')
idx_multiallele = exac[idx].index
idx_singleallele = exac[~idx].index

expanded_list = []
for i in idx_multiallele:
    expand = exac.loc[i, cols].str.split(',', expand=True).T
    fill = pd.DataFrame(exac.loc[i, fillcols]).T.reset_index(drop=True)
    for j in expand.index[1:]:
        fill = fill.append(exac.loc[i, fillcols], ignore_index=True)
    expanded = pd.concat([expand, fill], axis=1)
    expanded_list.append(expanded)
expanded_exac = pd.concat(expanded_list, ignore_index=True)

df = pd.concat([exac.loc[idx_singleallele, :], expanded_exac], ignore_index=True)

mincols = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'AF', 'AC', 'AN']
df.to_csv('exac.expanded.r1.txt', sep='\t', index=False)
df.loc[:, mincols].to_csv('exac.expanded.min.r1.txt', sep='\t', index=False)
