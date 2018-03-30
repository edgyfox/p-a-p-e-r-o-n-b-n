#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 01:30:56 2018
OBJECTIVE: DRUG VECTOR APPLICATION ON TWO FAULTS
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
df_drug.columns=["Drugs"]
df_drug["Values"]=[0] * len(df_drug.index)

drugv=list(df_drug["Values"])

#creating output file
cols=["Drug Vector"]
for i in range(1,28):
    for j in range(i+1,28):
        cols=cols+[str(i)+","+str(j)]
output_drugtwo=pd.DataFrame(columns=cols)

k=0
start_time=time.clock()
while True:
    output_drugtwo.loc[k,"Drug Vector"]=' '.join(map(str,drugv))
    l=1
    print("Drug scenario:",drugv,"starts")
    for i in range(1,28):
        for j in range(i+1,28):
            drugpath.pathway([i,j],drugv,inpv,pathv,outv)
            inpv=unq
            output_drugtwo.iloc[k,l]=en.encode(outv)
            l=l+1
    print("%0.3f"%(time.clock()-start_time),"Drug scenario:",drugv,"ends")
    drugv=cmb.combination(drugv)
    if drugv==False:
        break
    k=k+1
    
print("Execution time: ","%0.3f"%(time.clock()-start_time)," seconds")

#print(output_drugtwo)

#write to outs/output_drugtwo.csv
output_drugtwo.to_csv("outs/output_drugtwo.csv")
