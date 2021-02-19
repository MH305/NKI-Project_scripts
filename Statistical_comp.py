from scipy import stats
import pandas as pd
import numpy as np

organ_1 = pd.read_csv('/Users/mashrurhaidernew/Desktop/MacPro/NKI/Speed_prof/speed_profiles/Liver_speed_prof.csv', sep=',', skip_blank_lines=True)
organ_2 = pd.read_csv('/Users/mashrurhaidernew/Desktop/MacPro/NKI/Speed_prof/speed_profiles/Skeletal_speed_prof.csv', sep=',', skip_blank_lines=True)
organ_1 = organ_1.rename(columns={"Unnamed: 0": "transcript"})
organ_2 = organ_2.rename(columns={"Unnamed: 0": "transcript"} )

#find unique transcripts in each
trans_org1 = organ_1['transcript'].unique().tolist()
#
trans_org2 = organ_2['transcript'].unique().tolist()

#find overlapping transcripts
def common(list1,list2):
	return list(set(list1) & set(list2))

#getting all the positions from the data
pos_perc = []
transcripts = common(trans_org1, trans_org2)

for col in organ_1.columns:
	pos_perc.append(col)

pos_perc.remove('transcript')

#do the analysis and return the D and P value append the D and P value to a dataFrame and sort it by the D value
kidn_livr = pd.DataFrame(index=transcripts, columns=['Distance', 'P_value'])
for trans in transcripts:
#First organ data
	org1 = organ_1[organ_1['transcript'] == trans]
	org1 = org1.iloc[:,1:]
	#org1 = org1.transpose()
	org1 = org1.to_numpy()
	org1 = org1.flatten()
#second organ data
	org2 = organ_2[organ_2['transcript'] == trans]
	org2 = org2.iloc[:,1:]
	#org2 = org2.transpose()
	org2 = org2.to_numpy()
	org2 = org2.flatten()
	dist, p_val = stats.ks_2samp(org1, org2)
	kidn_livr.loc[trans, 'Distance'] = dist
	kidn_livr.loc[trans, 'P_value'] = p_val

kidn_livr = kidn_livr.sort_values(by='Distance')
print(kidn_livr)

##Do svm ??


kidn_livr.to_csv(r'/Users/mashrurhaidernew/Desktop/MacPro/NKI/Speed_prof/Liver & Skeletal comparison.csv', encoding='utf-8-sig')
df_new = pd.read_csv('/Users/mashrurhaidernew/Desktop/MacPro/NKI/Speed_prof/Liver & Skeletal comparison.csv', sep=',', skip_blank_lines=True)
print(df_new)

