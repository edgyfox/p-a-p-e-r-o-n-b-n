#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 00:37:02 2018
OBJECTIVE: DRUG VECTOR APPLICATION ON THREE FAULTS
@author: arghanandan
"""

import pathway_drugged as drugpath
import pandas as pd
import combination as cmb
import encoder as en
import matplotlib.pyplot as plt
import time

start_time=time.clock()

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
nodes=28

#reading drug file
df_drug=pd.read_csv("ins/drug.csv",header=None)
df_drug.columns=["Drugs"]
df_drug["Values"]=[0] * len(df_drug.index)

drugv=list(df_drug["Values"])

#creating output file
cols=["Drug Vector"]
for i in range(1,nodes):
    for j in range(i+1,nodes):
        for k in range(j+1,nodes):
            cols=cols+[str(i)+","+str(j)+","+str(k)]
output_drugthree=pd.DataFrame(columns=cols)

m=0
while True:
    output_drugthree.loc[m,"Drug Vector"]=' '.join(map(str,drugv))
    l=1
    print("Drug scenario:",drugv,"starts")
    for i in range(1,nodes):
        for j in range(i+1,nodes):
            for k in range(j+1,nodes):
                drugpath.pathway([i,j,k],drugv,inpv,pathv,outv)
                inpv=unq
                output_drugthree.iloc[m,l]=en.encode(outv)
#                print("Drug scenario:",m+1,drugv,"  Fault scenario:",l,"   Fault locations:",[i,j,k])
                l=l+1
    print("%0.3f"%(time.clock()-start_time),"Drug scenario:",drugv,"ends")
    drugv=cmb.combination(drugv)
    if drugv==False:
        break
    m=m+1

#output_drugthree.to_csv("outs/output_drugthree.csv")

print("Execution time: ","%0.3f"%(time.clock()-start_time)," seconds")
