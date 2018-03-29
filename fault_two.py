#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 01:15:56 2018
OBJECTIVE: INDUCE TWO FAULTS IN BOOLEAN NETWORK
@author: arghanandan
"""

import pandas as pd
import pathway_normal as pth
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

#output_double_fault dataframe
output_2f=pd.DataFrame(columns=["output proteins"])
output_2f["output proteins"]=out["proteins"]

#executing BN at ith gate
for i in range(1,28):
    for j in range(i+1,28):
        pth.pathway([i,j],inpv,pathv,outv)
        output_2f[str(i)+","+str(j)]=outv
        inpv=unq
            
#print(output_2f)

#write to .csv file   
output_2f.to_csv("outs/output_2f.csv")