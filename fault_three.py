#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 15:16:23 2018

@author: arghanandan
"""

import pandas as pd
import pathway_normal as pth
import time
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
output_3f=pd.DataFrame(columns=["output proteins"])
output_3f["output proteins"]=out["proteins"]

start_time=time.clock()
#executing BN at ith gate
for i in range(1,28):
    for j in range(i+1,28):
        for k in range(j+1,28):
            pth.pathway([i,j,k],inpv,pathv,outv)
            output_3f[str(i)+","+str(j)+","+str(k)]=outv
            inpv=unq
            
print("Execution time: ","%0.3f"%(time.clock()-start_time)," seconds")
#print(output_3f)

#write to .csv file   
output_3f.to_csv("outs/output_3f.csv")