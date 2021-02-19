#!/usr/bin/env python

import glob
import os
import numpy as np
import types
import pandas as pd
#import matplotlib.pyplot as plt


#Read csv
df_data_1 = pd.read_csv('/Users/mashrurhaidernew/Desktop/MacPro/NKI/Speed_prof/Kidney/Kidney_Proc_new.csv', sep=',', skip_blank_lines=True)
#check columns in the data
for cols in df_data_1.columns:
    print(cols)
print(df_data_1)
#get the relevant features
df_data_1 = df_data_1[['transcript', 'Time Point', 'pos_perc', 'avg_cds_norm_count']]


#list all the transcripts in the data
transcripts = df_data_1['transcript'].unique().tolist()
#df = df_data_1[df_data_1['transcript'] == 'ENSMUST00000217608']
#df = pd.pivot_table(df, index=['transcript','pos_perc'], columns='Time Point', values='avg_cds_norm_count', fill_value=0)#do .reset_index() here
#print(df)
pos_perc = df_data_1['pos_perc'].unique().tolist()


#reducing the features to those that are need for the elongation calculation
def reduce_data(df_data_1):
    unique_time = df_data_1['Time Point'].unique().tolist()
    df = pd.pivot_table(df_data_1, index=['transcript','pos_perc'], columns='Time Point', values='avg_cds_norm_count', fill_value=0)#do .reset_index() here
    df = df.reset_index()
    df2 = df.drop(['transcript', 'pos_perc'], axis = 1)
    df_new = pd.DataFrame(columns = unique_time)
    for col in unique_time:
        df_new[col] = df[col]
    #Total sum per column: 
    df_new.loc['col_sum',:]= df2.sum(axis=0)
    #Total sum per row: 
    df_new.loc[:,'row_sum'] = df2.sum(axis=1)
    df_new['pos_perc'] = df['pos_perc']  
    #df_new = df_new.set_index('position')
    df_new = pd.DataFrame(df_new)
    #unique_time = pd.DataFrame(data = unique_time, columns = ['time'])
    position = pd.DataFrame(df['pos_perc'])
    position = position.reset_index()
    position = position.iloc[:,-1]
    return df_new, position


#normalization function for accounting for the differences in transcript length
def bipercentile_norm(df):
    df = pd.DataFrame(df)
    col_sum = df.iloc[-1,:-2]
    #col_sum = col_sum.drop(columns = ['index'])
    col_sum = col_sum.iloc[0:6].to_numpy()
    row_sum = df.iloc[:-1,-2].to_numpy()
    s = (len(row_sum), len(col_sum))
    new_matrix = np.zeros(s)
    df = df.iloc[:-1, :-2]
    df = df.to_numpy()
    #print(df)
    #print(unique_time)
    #print(col_sum)
    #print(new_matrix)

    for row in range(len(row_sum)):
        for col in range(len(col_sum)):
            if df[row, col] != 0:
                new_matrix[row, col] = ((df[row, col]/col_sum[col])*(df[row, col]/row_sum[row]))/((df[row, col]/col_sum[col])+(df[row, col]/row_sum[row]))
    return new_matrix        


#pass the bipercentile function and speed profile calc inside the for loop
#changing back into dataframe
#df_2 = bipercentile_norm(df_new)
def speed_profile(df_2):
    unique_time = df_data_1['Time Point'].unique().tolist()    
    norm_df = pd.DataFrame(df_2, columns = unique_time)
    #norm_df = norm_df.replace(0, np.NaN)
    #calc position and time means
    norm_df.loc['col_mean',:] = norm_df.mean(axis = 0, skipna=True)
    norm_df.loc[:-1,'row_mean'] = norm_df.mean(axis = 1, skipna=True)
    norm_df['pos_perc'] = position
    #append the data to one file with the transcript name in each column and pos_perc index
    return norm_df


#running functions
speed_prof_df = pd.DataFrame(columns = pos_perc)
unique_time = df_data_1['Time Point'].unique().tolist()
length_u = len(unique_time)
remaining_trans = transcripts
for trans in transcripts:
    df = df_data_1[df_data_1['transcript'] == trans]
    new_time = df['Time Point'].unique().tolist()
    length_n = len(new_time)
    remaining_trans.remove(trans)
    print('remaining num: ' + str(len(remaining_trans)))
    print(new_time)
    if length_u == length_n:
        print('working on: ' + str(trans))
        print(new_time)
        df = df.sort_values(by='pos_perc')
        df_new, position = reduce_data(df)
        new_matrix = bipercentile_norm(df_new)
        norm_df = speed_profile(new_matrix)
        norm_df = norm_df.sort_values(by='pos_perc')
        row_mean = norm_df.iloc[:-1, -2]
        print(row_mean)
        speed_prof_df.loc[trans, :] = row_mean
    else:
        print('Ignored: ' + str(trans))
print('Final DF: ')
print(speed_prof_df)


#save file to csv
# Import the lib
print('saving')
speed_prof_df = pd.DataFrame(speed_prof_df)
speed_prof_df.to_csv(r'/Users/mashrurhaidernew/Desktop/MacPro/NKI/Speed_prof/Kidney/Kidney_speed_prof.csv', encoding='utf-8-sig')
print('saved')


