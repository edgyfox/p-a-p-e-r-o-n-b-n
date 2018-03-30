#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 15:24:39 2018
OBJECTIVE: INDUCE FOUR FAULTS IN BOOLEAN NETWORK
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
output_4f=pd.DataFrame(columns=["output proteins"])
output_4f["output proteins"]=out["proteins"]

start_time=time.clock()
#executing BN at ith gate
for i in range(1,28):
    for j in range(i+1,28):
        for k in range(j+1,28):
            for l in range(k+1,28):
                pth.pathway([i,j,k,l],inpv,pathv,outv)
                output_4f[str(i)+","+str(j)+","+str(k)+","+str(l)]=outv
                inpv=unq
    print("Combinations of " + str(i) + " ends.")
            
print("Execution time: ","%0.3f"%(time.clock()-start_time)," seconds")
#print(output_4f)

#write to .csv file   
output_4f.to_csv("outs/output_4f.csv")