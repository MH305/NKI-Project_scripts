#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt

# read the data from different files into one data frame
'''
path = r'/Users/mashrurhaidernew/Desktop/NKI/data'    # use your path
all_files = glob.glob(os.path.join(path, "*.counts"))     # advisable to use os.path.join as this makes concatenation OS independent
df_from_each_file = (pd.read_csv(f, sep='\t') for f in all_files)
data = pd.concat((pd.read_csv(f, sep='\t') for f in all_files))
data = pd.DataFrame(data)
'''
columns = ['transcript', 'gene', 'sample', 'pos_perc', 'cds_norm_count', 'alpha']

# opens each individual file
data_0 = pd.read_csv('/Users/mashrurhaidernew/Downloads/Skeletal_1/26-32-0TX.profile.percentage.filtered.tsv', sep='\t', skip_blank_lines=True, names=columns, header=None)
data_15 = pd.read_csv('/Users/mashrurhaidernew/Downloads/Skeletal_1/26-32-15TX.profile.percentage.filtered.tsv', sep='\t', skip_blank_lines=True, names=columns, header=None)
data_30 = pd.read_csv('/Users/mashrurhaidernew/Downloads/Skeletal_1/26-32-30TX.profile.percentage.filtered.tsv', sep='\t', skip_blank_lines=True, names=columns, header=None)
data_45 = pd.read_csv('/Users/mashrurhaidernew/Downloads/Skeletal_1/26-32-45TX.profile.percentage.filtered.tsv', sep='\t', skip_blank_lines=True, names=columns, header=None)
data_60 = pd.read_csv('/Users/mashrurhaidernew/Downloads/Skeletal_1/26-32-60TX.profile.percentage.filtered.tsv', sep='\t', skip_blank_lines=True, names=columns, header=None)
data_300 = pd.read_csv('/Users/mashrurhaidernew/Downloads/Skeletal_1/26-32-300TX.profile.percentage.filtered.tsv', sep='\t', skip_blank_lines=True, names=columns, header=None)

# Formatting data and calculating mean for each transcript
ax_0 = pd.pivot_table(data_0.reset_index(), index=['pos_perc', 'transcript'], columns='sample',
                      values='cds_norm_count', fill_value=0)
ax_15 = pd.pivot_table(data_15.reset_index(), index=['pos_perc', 'transcript'], columns='sample',
                       values='cds_norm_count', fill_value=0)
ax_30 = pd.pivot_table(data_30.reset_index(), index=['pos_perc', 'transcript'], columns='sample',
                       values='cds_norm_count', fill_value=0)
ax_45 = pd.pivot_table(data_45.reset_index(), index=['pos_perc', 'transcript'], columns='sample',
                       values='cds_norm_count', fill_value=0)
ax_60 = pd.pivot_table(data_60.reset_index(), index=['pos_perc', 'transcript'], columns='sample',
                       values='cds_norm_count', fill_value=0)
ax_300 = pd.pivot_table(data_300.reset_index(), index=['pos_perc', 'transcript'], columns='sample',
                        values='cds_norm_count', fill_value=0)


# mean with zero elements and adding a new column
ax_0['avg_cds_norm_count'] = ax_0.mean(axis=1, skipna=False)
ax_0['Time Point'] = 0
ax_0 = ax_0.reset_index()
ax_15['avg_cds_norm_count'] = ax_15.mean(axis=1, skipna=False)
ax_15['Time Point'] = 15
ax_15 = ax_15.reset_index()
ax_30['avg_cds_norm_count'] = ax_30.mean(axis=1, skipna=False)
ax_30['Time Point'] = 30
ax_30 = ax_30.reset_index()
ax_45['avg_cds_norm_count'] = ax_45.mean(axis=1, skipna=False)
ax_45['Time Point'] = 45
ax_45 = ax_45.reset_index()
ax_60['avg_cds_norm_count'] = ax_60.mean(axis=1, skipna=False)
ax_60['Time Point'] = 60
ax_60 = ax_60.reset_index()
ax_300['avg_cds_norm_count'] = ax_300.mean(axis=1, skipna=False)
ax_300['Time Point'] = 300
ax_300 = ax_300.reset_index()

print(ax_0)
'''
#calculating cumulative avg_cds_norm_count
ax_0['cum_sum_cds_count'] = ax_0['avg_cds_norm_count'].cumsum()
ax_15['cum_sum_cds_count'] = ax_15['avg_cds_norm_count'].cumsum()
ax_30['cum_sum_cds_count'] = ax_30['avg_cds_norm_count'].cumsum()
ax_45['cum_sum_cds_count'] = ax_45['avg_cds_norm_count'].cumsum()
ax_60['cum_sum_cds_count'] = ax_60['avg_cds_norm_count'].cumsum()
ax_300['cum_sum_cds_count'] = ax_300['avg_cds_norm_count'].cumsum()
ax_0 = ax_0.reset_index()
ax_15 = ax_15.reset_index()
ax_30 = ax_30.reset_index()
ax_45 = ax_45.reset_index()
ax_60 = ax_60.reset_index()
ax_300 = ax_300.reset_index()
'''

#creating one DF with the new column and create a new csv file
df = ax_0.append([ax_15, ax_30, ax_45, ax_60, ax_300], sort=True).fillna(value=0)
df = pd.DataFrame(df)
df.to_csv(r'/Users/mashrurhaidernew/Downloads/Skeletal_1/Skeletal_Proc_new.csv', encoding='utf-8-sig')
df_new = pd.read_csv('/Users/mashrurhaidernew/Downloads/Skeletal_1/Skeletal_Proc_new.csv', sep=',', skip_blank_lines=True)
print(df_new)

'''
#subplot according to timepoint

fig, axes = plt.subplots(nrows=2, ncols=3)
plot1_a = ax_0.plot(x='position', y='avg_cds_norm_count', ax=axes[0, 0], label='0', legend=False, color="red", sharex=True, sharey=True, xlim=(0,2100), ylim=(0,0.004))
plot2_b = ax_15.plot(x='position', y='avg_cds_norm_count', ax=axes[0, 1], label='15', legend=False, color="green", sharex=True, sharey=True, xlim=(0,2100), ylim=(0,0.004))
plot3_c = ax_30.plot(x='position', y='avg_cds_norm_count', ax=axes[0, 2], label='30', legend=False, color="blue", sharex=True, sharey=True, xlim=(0,2100), ylim=(0,0.004))
plot4_d = ax_45.plot(x='position', y='avg_cds_norm_count', ax=axes[1, 0], label='45', legend=False, color="orange",sharex=True, sharey=True, xlim=(0,2100), ylim=(0,0.004))
plot5_e = ax_60.plot(x='position', y='avg_cds_norm_count', ax=axes[1, 1], label='60', legend=False, color="magenta", sharex=True, sharey=True, xlim=(0,2100), ylim=(0,0.004))
plot6_f = ax_300.plot(x='position', y='avg_cds_norm_count', ax=axes[1, 2], label='300', legend=False, color="purple", sharex=True, sharey=True, xlim=(0,2100), ylim=(0,0.004))

line_labels = ["Time 0", "Time 15", "Time 30", "Time 45", "Time 60", "Time 300"]

fig.legend([plot1_a, plot2_b, plot3_c, plot4_d, plot5_e, plot6_f], labels=line_labels, loc="center right", borderaxespad=0.1)
plt.subplots_adjust(right=0.85, wspace=0.4, hspace=0.4)
plt.savefig('/Volumes/Haider/Data_haider_NKI/RiboSAK_results/Liver_prof_count/Liver_graph.pdf')


#cumulative plot
x_0 = ax_0['position']
y_0 = ax_0['cum_sum_cds_count']
x_15 = ax_15['position']
y_15 = ax_15['cum_sum_cds_count']
x_30 = ax_30['position']
y_30 = ax_30['cum_sum_cds_count']
x_45 = ax_45['position']
y_45 = ax_45['cum_sum_cds_count']
x_60 = ax_60['position']
y_60 = ax_60['cum_sum_cds_count']
x_300 = ax_300['position']
y_300 = ax_300['cum_sum_cds_count']

plt.plot(x_0, y_0)
plt.plot(x_15, y_15)
plt.plot(x_30, y_30)
plt.plot(x_45, y_45)
plt.plot(x_60, y_60)
plt.plot(x_300, y_300)
plt.legend(['time = 0', 'time = 15', 'time = 30', 'time = 45', 'time = 60', 'time = 300'], loc='upper right')
plt.show()
plt.savefig('/Volumes/Haider/Data_haider_NKI/RiboSAK_results/Liver_prof_count/Liver_graph_cumsum.pdf')
'''
















