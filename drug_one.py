#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 01:23:56 2018
OBJECTIVE: DRUG VECTOR APPLICATION ON SINGLE FAULT
@author: arghanandan
"""

import pathway_drugged as drugpath
import pandas as pd
import combination as cmb
import encoder as en
import matplotlib.pyplot as plt
import time

#common settings
plt.style.use("ggplot")
pd.set_option("display.max_columns",None)
pd.set_option("display.max_rows",None)
        
#reading protein file
inp=pd.read_csv("ins/inp.csv",delimiter=",",index_col=0)
path=pd.read_csv("ins/path.csv",delimiter=",",index_col=0)
out=pd.read_csv("ins/out.csv",delimiter=",",index_col=0)
inp["values"]=[0]*len(inp.index)
path["values"]=[0]*len(path.index)
out["values"]=[0]*len(out.index)

#input,pathway and output vectors
f=open("outs/output_unq.txt","r")
unq=f.readline()
unq=list(map(int,unq.split(" ")))
inpv=unq
pathv=list(path["values"])
outv=list(out["values"])

#reading drug file
df_drug=pd.read_csv("ins/drug.csv",header=None)
df_drug.columns=["drugs"]
df_drug["values"]=[0] * len(df_drug.index)

drugv=list(df_drug["values"])

#creating output file
cols=["drug vector"] + [i for i in range(28)]
output_drugone=pd.DataFrame(columns=cols)

j=0
start_time=time.clock()
while True:
    encoded=[]
    print("Drug scenario:",drugv,"starts")
    for i in range(28):
        drugpath.pathway([i],drugv,inpv,pathv,outv)
        encoded.append(float(en.encode(outv)))
        inpv=unq
    print("%0.3f"%(time.clock()-start_time),"Drug scenario:",drugv,"ends")
    output_drugone.loc[j,"drug vector"]=' '.join(map(str,drugv))
    output_drugone.iloc[j,1:]=encoded
    drugv=cmb.combination(drugv)
    if drugv==False:
        break
    j=j+1

print("Execution time: ","%0.3f"%(time.clock()-start_time)," seconds")
#print(output_drugone)

#write to output_drugone.csv
output_drugone.to_csv("outs/output_drugone.csv")