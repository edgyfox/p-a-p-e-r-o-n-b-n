#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 13:54:08 2018
OBJECTIVE: COUNT THE OCCURENCE OF OUTPUT VECTORS IN EACH FAULT DATAFRAME
@author: arghanandan
"""

import pandas as pd

count_df=pd.read_csv("outs/outv.csv",index_col=0,header=None)
count_df.columns=["output vector"]
count=[0]*len(count_df.index)

k=0
for i in range(4):
    ifile="outs/output_"+str(i+1)+"f.csv"
    df=pd.read_csv(ifile,index_col=0)
    count_df["fault" + str(i+1)] = count
    for j in range(len(df.columns)-1):
        out=" ".join(map(str,list(df.iloc[:,j+1])))
        row=list(count_df["output vector"]).index(out)
        count_df.iloc[row,i+1]=count_df.iloc[row,i+1] + 1
        
count_df.loc[20]=["total"] + [0]*4
for i in range(4):
    count_df.iloc[20,i+1]=count_df.iloc[:,i+1].sum()
    
#write dataframe to file
count_df.to_csv("outs/output_count.csv")