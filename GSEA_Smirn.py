import glob
import os
import numpy as np
import types
import pandas as pd
import matplotlib.pyplot as plt
import plotly.io as pio
from matplotlib_venn import venn3, venn3_circles
from functools import reduce
import seaborn as sns
import seaborn


## Smirnoff test results
kidney_liver_sdf = pd.read_csv('Kidney_Liver_comparison.csv', sep=',', skip_blank_lines=True)
kidney_skeletal_sdf = pd.read_csv('Kidney_Skeletal_comparison.csv', sep=',', skip_blank_lines=True)
liver_skeletal_sdf = pd.read_csv('Liver_Skeletal_comparison.csv', sep=',', skip_blank_lines=True)

# Filtering out transcripts with a P_value of more than 0.05
kidney_liver_sdf_filtered = kidney_liver_sdf[kidney_liver_sdf.P_value < 0.05]
kidney_skeletal_sdf_filtered = kidney_skeletal_sdf[kidney_skeletal_sdf.P_value < 0.05]
liver_skeletal_sdf_filtered = liver_skeletal_sdf[liver_skeletal_sdf.P_value < 0.05]
print(len(kidney_liver_sdf_filtered))
# Sorting transcripts based on distance
kidney_liver_sdf_filtered = kidney_liver_sdf_filtered.sort_values(by='Distance', ascending=False)
kidney_liver_sdf_filtered.rename(columns={'Unnamed: 0': 'Transcript'}, inplace=True)
kidney_skeletal_sdf_filtered = kidney_skeletal_sdf_filtered.sort_values(by='Distance', ascending=False)
kidney_skeletal_sdf_filtered.rename(columns={'Unnamed: 0': 'Transcript'}, inplace=True)
liver_skeletal_sdf_filtered = liver_skeletal_sdf_filtered.sort_values(by='Distance', ascending=False)
liver_skeletal_sdf_filtered.rename(columns={'Unnamed: 0': 'Transcript'}, inplace=True)


# calculating the mean and standard deviation of the results
kidney_liver_mean = kidney_liver_sdf_filtered.iloc[:, 1].mean()
kidney_liver_stdev = kidney_liver_sdf_filtered.iloc[:, 1].std()
kidney_liver_median = kidney_liver_sdf_filtered.iloc[:, 1].median()

kidney_skeletal_mean = kidney_skeletal_sdf_filtered.iloc[:, 1].mean()
kidney_skeletal_stdev = kidney_skeletal_sdf_filtered.iloc[:, 1].std()
kidney_skeletal_median = kidney_skeletal_sdf_filtered.iloc[:, 1].median()

liver_skeletal_mean = liver_skeletal_sdf_filtered.iloc[:, 1].mean()
liver_skeletal_stdev = liver_skeletal_sdf_filtered.iloc[:, 1].std()
liver_skeletal_median = liver_skeletal_sdf_filtered.iloc[:, 1].median()

# filtering transcripts above the mean distance in each instance
kidney_liver_mean_filter = kidney_liver_sdf_filtered[kidney_liver_sdf_filtered.Distance > 0.315162]
kidney_skeletal_mean_filter = kidney_skeletal_sdf_filtered[kidney_skeletal_sdf_filtered.Distance > 0.291418]
liver_skeletal_mean_filter = liver_skeletal_sdf_filtered[liver_skeletal_sdf_filtered.Distance > 0.285379]

## Top 10 transcripts with the highest distance in each comparison
plotting_kid_liv = kidney_liver_sdf_filtered.iloc[:5, :].set_index("Transcript")
plotting_kid_skel = kidney_skeletal_sdf_filtered.iloc[:5, :].set_index("Transcript")
plotting_liv_skel = liver_skeletal_sdf_filtered.iloc[:5, :].set_index("Transcript")

plotting_kid_liv = plotting_kid_liv.T
plotting_kid_skel = plotting_kid_skel.T
plotting_liv_skel = plotting_liv_skel.T

plotting_kid_liv.rename(
	columns={"ENSMUST00000068322": 'SEC14L3', "ENSMUST00000111827": "SLCO1A6-201", "ENSMUST00000203597": "SLCO1B2-204",
	         "ENSMUST00000077788": "TNFAIP8L1", "ENSMUST00000058777": "ANGPTL8"}, inplace=True)
plotting_kid_skel.rename(
	columns={"ENSMUST00000037964": "TXLNB", "ENSMUST00000092048": "NACA", "ENSMUST00000168279": "SLC15A2",
	         "ENSMUST00000131723": "SLC47A1-002", "ENSMUST00000108486": "ATP2A3"}, inplace=True)
plotting_liv_skel.rename(
	columns={"ENSMUST00000092048": "NACA", "ENSMUST00000087050": "COL4a4", "ENSMUST00000211150": "GYS1",
	         "ENSMUST00000102970": "PRKCQ", "ENSMUST00000088246": "TUBA3A"}, inplace=True)

plotting_kid_liv = plotting_kid_liv.T.reset_index(drop=False)
plotting_kid_skel = plotting_kid_skel.T.reset_index(drop=False)
plotting_liv_skel = plotting_liv_skel.T.reset_index(drop=False)

print("kid:skel")
print(plotting_kid_skel)
print("liv:skel")
print(plotting_liv_skel)
print("kid:liv")
print(plotting_kid_liv)

# Plot the different comparison with X= transcript and Y=Distance
fig2, axs2 = plt.subplots(nrows=3, sharey=True)

fig_1 = sns.scatterplot(ax=axs2[0], x=plotting_kid_liv['Transcript'], y=plotting_kid_liv['Distance'], color='navy')
fig_2 = sns.scatterplot(ax=axs2[1], x=plotting_kid_skel['Transcript'], y=plotting_kid_skel['Distance'], color='red')
fig_3 = sns.scatterplot(ax=axs2[2], x=plotting_liv_skel['Transcript'], y=plotting_liv_skel['Distance'], color='green')

plt.show()

#


# Comparing the transcripts within each comparison. First we need to convert it to sets
kidney_liver_set = set(kidney_liver_sdf_filtered['Transcript'])
kidney_skeletal_set = set(kidney_skeletal_sdf_filtered['Transcript'])
liver_skeletal_set = set(liver_skeletal_sdf_filtered['Transcript'])

# Union and distinct transcripts
union_set = kidney_liver_set | kidney_skeletal_set | liver_skeletal_set
intersect_all = kidney_skeletal_set & kidney_liver_set & liver_skeletal_set
intersect_kidney = (kidney_liver_set & kidney_skeletal_set) - liver_skeletal_set
intersect_skeletal = (kidney_skeletal_set & liver_skeletal_set) - kidney_liver_set
intersect_liver = (liver_skeletal_set & kidney_liver_set) - kidney_skeletal_set

print("all set: ")
print(len(union_set))
print("intersect all:")
print(len(intersect_all))
print('Intersect for Liver: ')
print(intersect_liver)
print('Intersect for kidney: ')
print(intersect_kidney)
print('Intersect for skeletal: ')
print(intersect_skeletal)

# Venn diagram
venn3([liver_skeletal_set, kidney_liver_set, kidney_skeletal_set],
      set_labels=('liver & skeletal', 'kidney & liver', 'kidney & skeletal'),
      set_colors=('SteelBlue', 'Seagreen', 'SlateGrey'), alpha=0.8)

plt.show()

# Top 10 transcripts with the highest distance in each comparison
t_trans_kid_liv = kidney_liver_sdf_filtered.iloc[:5, :]
t_trans_kid_skel = kidney_skeletal_sdf_filtered.iloc[:5, :]
t_trans_liv_skel = liver_skeletal_sdf_filtered.iloc[:5, :]

# Saving the file of top transcripts


print('Top 10 transcripts with the highest distance in liver_skeletal comparison: ')
print(t_trans_liv_skel)
print('Top 10 transcripts with the highest distance in kidney_skeletal comparison: ')
print(t_trans_kid_skel)
print('Top 10 transcripts with the highest distance in kidney_liver comparison: ')
print(t_trans_kid_liv)
'''

## GSEA results
kidney_liver_gdf = pd.read_csv('GSEA_Kidney_Liver_results.csv', sep=',', skip_blank_lines=True)
kidney_skeletal_gdf = pd.read_csv('GSEA_Kidney_Skeletal_results.csv', sep=',', skip_blank_lines=True)
liver_skeletal_gdf = pd.read_csv('GSEA_Liver_Skeletal_results.csv', sep=',', skip_blank_lines=True)

kidney_liver_gdf.rename(columns={'p.adjust': 'p_adjust'}, inplace=True)
kidney_skeletal_gdf.rename(columns={'p.adjust': 'p_adjust'}, inplace=True)
liver_skeletal_gdf.rename(columns={'p.adjust': 'p_adjust'}, inplace=True)

# Filtering out transcripts with a P_value of more than 0.05
kidney_liver_gdf_filtered = kidney_liver_gdf[kidney_liver_gdf.p_adjust < 0.05]
kidney_skeletal_gdf_filtered = kidney_skeletal_gdf[kidney_skeletal_gdf.p_adjust < 0.05]
liver_skeletal_gdf_filtered = liver_skeletal_gdf[liver_skeletal_gdf.p_adjust < 0.05]

# sorting the data based on NES value
kidney_liver_gdf_filtered = kidney_liver_gdf_filtered.sort_values(by='NES', ascending=False)
kidney_skeletal_gdf_filtered = kidney_skeletal_gdf_filtered.sort_values(by='NES', ascending=False)
liver_skeletal_gdf_filtered = liver_skeletal_gdf_filtered.sort_values(by='NES', ascending=False)

# Comparing the transcripts within each comparison. First we need to convert it to sets
kidney_liver_gset = set(kidney_liver_gdf_filtered['Description'])
kidney_skeletal_gset = set(kidney_skeletal_gdf_filtered['Description'])
liver_skeletal_gset = set(liver_skeletal_gdf_filtered['Description'])

# Union and distinct gsea
union_gset = kidney_liver_gset | kidney_skeletal_gset | liver_skeletal_gset
intersect_gall = kidney_skeletal_gset & kidney_liver_gset & liver_skeletal_gset
intersect_gkidney = (kidney_liver_gset & kidney_skeletal_gset)
intersect_gskeletal = (kidney_skeletal_gset & liver_skeletal_gset)
intersect_gliver = (liver_skeletal_gset & kidney_liver_gset)


print('Common gene sets for Liver_skeletal and Kidney_liver GSEA: ')
print(intersect_gliver)
print('Common gene sets for Kidney_skeletal and Kidney_liver GSEA: ')
print(intersect_gkidney)
print('Common gene sets for Kidney_skeletal and Liver_skeletal GSEA: ')
print(intersect_gskeletal)
print('Common gene sets for all: ')
print(intersect_gall)
print(kidney_liver_gdf_filtered)

#printing the important GO's
print('GSEA result for kidney and liver, with significant p value: ')
print(kidney_liver_gdf_filtered.iloc[:, 2:])
print('GSEA result for kidney and skeletal, with significant p value: ')
print(kidney_skeletal_gdf_filtered.iloc[:, 2:])
print('GSEA result for liver and skeletal, with significant p value: ')
print(liver_skeletal_gdf_filtered.iloc[:, 2:])

#writing the results of the GSEA in the excel and csv
kidney_liver_gdf_filtered = kidney_liver_gdf_filtered.round(2)
liver_skeletal_gdf_filtered = liver_skeletal_gdf_filtered.round(2)
kidney_skeletal_gdf_filtered = kidney_skeletal_gdf_filtered.round(2)

print('kidney_liver')
print(kidney_liver_gdf_filtered)
print('kidney_skel')
print(kidney_skeletal_gdf_filtered)
print('liver_skel')
print(liver_skeletal_gdf_filtered)

kidney_skeletal_gdf_filtered.to_csv('Kidney_skeletal_gsea.csv', index=False)
kidney_liver_gdf_filtered.to_csv('Kidney_liver_gsea.csv', index=False)
liver_skeletal_gdf_filtered.to_csv('Liver_skeletal_gsea.csv', index=False)



# plotting x=description y=NES
# Plot the different comparison with X= transcript and Y=Distance
fig_4 = plt.scatter(data_frame=kidney_liver_gdf_filtered,
                    x=kidney_liver_gdf_filtered['Description'], y=kidney_liver_gdf_filtered['NES'],
                    title='Kidney_liver_significant_transcripts', opacity=0.8)
fig_5 = plt.scatter(data_frame=kidney_skeletal_gdf_filtered,
                    x=kidney_skeletal_gdf_filtered['Description'], y=kidney_skeletal_gdf_filtered['NES'],
                    title='Kidney_skeletal_significant_transcripts', opacity=0.8)
fig_6 = plt.scatter(data_frame=liver_skeletal_gdf_filtered,
                    x=liver_skeletal_gdf_filtered['Description'], y=liver_skeletal_gdf_filtered['NES'],
                    title='Liver_skeletal_significant_transcripts', opacity=0.8)

#fig_4.show()
#fig_5.show()
#fig_6.show()

'''
## Find the speed profiles of relevent transcripts and plotting
kidney_speed_prof = pd.read_csv('Kidney_speed_prof.csv', sep=',', skip_blank_lines=True)
liver_speed_prof = pd.read_csv('Liver_speed_prof.csv', sep=',', skip_blank_lines=True)
skeletal_speed_prof = pd.read_csv('Skeletal_speed_prof.csv', sep=',', skip_blank_lines=True)

kidney_speed_prof.rename(columns={'Unnamed: 0': 'Transcript'}, inplace=True)
liver_speed_prof.rename(columns={'Unnamed: 0': 'Transcript'}, inplace=True)
skeletal_speed_prof.rename(columns={'Unnamed: 0': 'Transcript'}, inplace=True)

print(kidney_speed_prof)


# Get the top ten transcripts based on distance

# Get the transcript speed profile from each organ

# t_trans_kid_liv = list(t_trans_kid_liv.Transcript)


def create_top_df(Organ1sp, Organ2sp, top_trans):
	# Changing the top_trans set will change the comparison being processed.
	Organ_1 = pd.DataFrame()
	Organ_2 = pd.DataFrame()

	organ_set_1 = set(Organ1sp['Transcript'])
	organ_set_2 = set(Organ2sp['Transcript'])
	top_set_comparison = set(top_trans['Transcript'])
	common_set_2 = organ_set_1 & top_set_comparison & organ_set_2
	Organ1 = Organ1sp.set_index('Transcript', drop=True)
	Organ2 = Organ2sp.set_index('Transcript', drop=True)

	for i in common_set_2:
		Organ_1[i] = Organ1.loc[i, '0':'100']
		Organ_2[i] = Organ2.loc[i, '0':'100']

	Organ_1["Position"] = Organ_1.index.astype(int)
	Organ_1.set_index('Position', inplace=True)

	Organ_2["Position"] = Organ_2.index.astype(int)
	Organ_2.set_index('Position', inplace=True)

	new_df1 = Organ_1
	new_df2 = Organ_2

	print("Being compared: " + str(top_trans))
	print(Organ_1, Organ_2)

	return Organ_1, Organ_2


# Plotting Kidney Liver top transcripts. Speed profile

top_speed_kidney, top_speed_liver = create_top_df(kidney_speed_prof, liver_speed_prof, t_trans_kid_liv)

top_speed_kidney['Organ'] = 'Kidney'
top_speed_liver['Organ'] = 'Liver'

top_speed_kidney.rename(
	columns={"ENSMUST00000068322": 'SEC14L3', "ENSMUST00000111827": "SLCO1A6-201", "ENSMUST00000203597": "SLCO1B2-204",
	         "ENSMUST00000077788": "TNFAIP8L1", "ENSMUST00000058777": "ANGPTL8"}, inplace=True)

top_speed_liver.rename(
	columns={"ENSMUST00000068322": 'SEC14L3', "ENSMUST00000111827": "SLCO1A6-201", "ENSMUST00000203597": "SLCO1B2-204",
	         "ENSMUST00000077788": "TNFAIP8L1", "ENSMUST00000058777": "ANGPTL8"}, inplace=True)

print("printing Kidney: ")
print(top_speed_kidney)
print("printing Liver: ")
print(top_speed_liver)

fig, axs = plt.subplots(3, 2)

# comparison between Kidney and Liver
f1 = sns.violinplot(ax=axs[0, 0], data=top_speed_kidney, inner='box', color='navy', saturation=0.6)
f2 = sns.violinplot(ax=axs[0, 0], data=top_speed_liver, inner='box', color='red', saturation=0.6)
axs[0, 0].set_ylabel("Read Counts")
plt.legend([f1, f2], ["Kidney", "Liver"])

f3 = sns.lineplot(ax=axs[0, 1], x='Position', y=top_speed_kidney['Sec14l3'], data=top_speed_kidney, color='navy',
                  marker='o', label="Kidney")
f4 = sns.lineplot(ax=axs[0, 1], x='Position', y=top_speed_liver['Sec14l3'], data=top_speed_liver, color='red',
                  marker='o', label="Liver")
axs[0, 1].legend(loc="upper right")

# Plotting top comparison of kidney and skeletal
top_speed_kidney, top_speed_skeletal = create_top_df(kidney_speed_prof, skeletal_speed_prof, t_trans_kid_skel)

top_speed_kidney['Organ'] = 'Kidney'
top_speed_skeletal['Organ'] = 'Skeletal'

print("printing Kidney: ")
print(top_speed_kidney)
print("printing Skeletal")
print(top_speed_skeletal)

top_speed_skeletal.rename(
	columns={"ENSMUST00000037964": "TXLNB", "ENSMUST00000092048": "NACA", "ENSMUST00000168279": "SLC15A2",
	         "ENSMUST00000131723": "SLC47A1-002", "ENSMUST00000108486": "ATP2A3"}, inplace=True)
top_speed_kidney.rename(
	columns={"ENSMUST00000037964": "TXLNB", "ENSMUST00000092048": "NACA", "ENSMUST00000168279": "SLC15A2",
	         "ENSMUST00000131723": "SLC47A1-002", "ENSMUST00000108486": "ATP2A3"}, inplace=True)

# Plotting
# comparison between Kidney and skeletal
f5 = sns.violinplot(ax=axs[1, 0], data=top_speed_kidney, inner='box', color='navy', saturation=0.6)
f6 = sns.violinplot(ax=axs[1, 0], data=top_speed_skeletal, inner='box', color='grey', saturation=0.6)
axs[1, 0].set_ylabel("Read Counts")

f7 = sns.lineplot(ax=axs[1, 1], x='Position', y=top_speed_kidney['Txlnb'], data=top_speed_kidney, color='navy',
                  marker='o', label="Kidney")
f8 = sns.lineplot(ax=axs[1, 1], x='Position', y=top_speed_skeletal['Txlnb'], data=top_speed_skeletal, color='grey',
                  marker='s', label="Skeletal")
axs[1, 1].legend(loc="upper right")

# Plotting top comparison of kidney and skeletal
top_speed_liver, top_speed_skeletal = create_top_df(liver_speed_prof, skeletal_speed_prof, t_trans_liv_skel)

top_speed_liver["Organ"] = "Liver"
top_speed_skeletal["Skeletal"] = "Skeletal"

top_speed_liver.rename(
	columns={"ENSMUST00000092048": "NACA", "ENSMUST00000087050": "COL4a4", "ENSMUST00000211150": "GYS1",
	         "ENSMUST00000102970": "PRKCQ", "ENSMUST00000088246": "TUBA3A"}, inplace=True)
top_speed_skeletal.rename(
	columns={"ENSMUST00000092048": "NACA", "ENSMUST00000087050": "COL4a4", "ENSMUST00000211150": "GYS1",
	         "ENSMUST00000102970": "PRKCQ", "ENSMUST00000088246": "TUBA3A"}, inplace=True)

print("printing Skeletal")
print(top_speed_skeletal)
print("printing Liver: ")
print(top_speed_liver)

# Plotting
# comparison between Liver and Skeletal muscle
f9 = sns.violinplot(ax=axs[2, 0], data=top_speed_liver, inner='box', color='red', saturation=0.6)
f10 = sns.violinplot(ax=axs[2, 0], data=top_speed_skeletal, inner='box', color='grey', saturation=0.6)
axs[2, 0].set_ylabel("Read Counts")

f11 = sns.lineplot(ax=axs[2, 1], x='Position', y=top_speed_liver['Naca'], data=top_speed_liver, color='red', marker='o',
                   label="Liver")
f12 = sns.lineplot(ax=axs[2, 1], x='Position', y=top_speed_skeletal['Naca'], data=top_speed_skeletal, color='grey',
                   marker='s', label="Skeletal")
axs[2, 1].legend(loc="upper right")

# Putting in legend


plt.show()

