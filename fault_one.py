#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 01:10:56 2018
INDUCING SINGLE FAULTS IN BOOLEAN NETWORK
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

#output_single_fault dataframe
output_1f=pd.DataFrame(columns=["output proteins"])
output_1f["output proteins"]=out["proteins"]

start_time=time.clock()

#executing BN at ith gate
for i in range(1,28):
    pth.pathway([i],inpv,pathv,outv)
    output_1f[str(i)]=outv
    inpv=[0,0,0,0,1]
    
print("Execution time: ","%0.3f"%(time.clock()-start_time)," seconds")

#print(output_1f)

#write to .csv file   
output_1f.to_csv("outs/output_1f.csv")