#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 19:55:56 2018
CHECKING INPUT VECTORS ON SIGNALLING PATHWAY
@author: arghanandan
"""

#for DataFrames
import pandas as pd
import pathway_normal as pth
import combination as cmb

#reading protein file
inp=pd.read_csv("ins/inp.csv",delimiter=",",index_col=0)
path=pd.read_csv("ins/path.csv",delimiter=",",index_col=0)
out=pd.read_csv("ins/out.csv",delimiter=",",index_col=0)
inp["values"]=[0]*len(inp.index)
path["values"]=[0]*len(path.index)
out["values"]=[0]*len(out.index)

#input,pathway and output vectors
inpv=list(inp["values"])
pathv=list(path["values"])
outv=list(out["values"])

#output_faultless dataframe
cols=["input"] + list(out["proteins"])
output_fl=pd.DataFrame(columns=cols)

#pathway_faultless dataframe
path_fl=pd.DataFrame(columns=["proteins","values"])
path_fl.Proteins=path["proteins"]

#output vector for each input combination
i=0
while True:
    pth.pathway([0],inpv,pathv,outv)
    output_fl.loc[i,"input"]=' '.join(map(str,inpv))
    if i==0:
        path_fl["values"]=pathv
    if outv==[0]*len(out.index):
        unq=list(map(int,(output_fl.loc[i,"input"]).split(" ")))
    output_fl.iloc[i,1:]=outv
    inpv=cmb.combination(inpv)
    if inpv==False:
        break
    i=i+1
    
#write to .csv file   
output_fl.to_csv("outs/output_fl.csv")

#write unique_input_vector to file
f=open("outs/output_unq.txt","w")
f.write(" ".join(map(str,unq)))
f.close()