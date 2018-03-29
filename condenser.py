#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 01:24:56 2018
OBJECTIVE: CALCULATE CONDENSED ENCODED WEIGHTS FOR ALL FAULT SCENARIOS PER DRUG
@author: arghanandan
"""

def condense(converted,convert,case):
    df_con=converted
    df=convert
    max=0
    col_df=len(df.columns[1:])
    len_df=len(df.index)
    for i in range(len_df):
        weight=0
        for j in range(col_df):
            weight=weight+df.iloc[i,j+1]
        df_con.iloc[i,1]=weight
        if i==0:
            max=weight
#    print(df_con)
    df_con.iloc[:,1]=(df_con.iloc[:,1]/max)*100
    return df_con