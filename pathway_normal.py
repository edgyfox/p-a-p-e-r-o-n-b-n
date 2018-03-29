#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  6 16:34:13 2018

@author: arghanandan
"""

#Boolean Network introducing fault at ith gate
def pathway(faultv,inpv,pathv,outv):
    
    if 1 in faultv:
        pathv[0]=1
    else:
        pathv[0]=inpv[0]
    
    if 2 in faultv:
        pathv[1]=1
    else:
        pathv[1]=inpv[0] or inpv[1]
    
    if 3 in faultv:
        pathv[2]=1
    else:
        pathv[2]=inpv[2]
    
    if 4 in faultv:
        pathv[3]=1
    else:
        pathv[3]=inpv[3]
    
    if 5 in faultv:
        pathv[4]=1
    else:
        pathv[4]=pathv[0] or pathv[1] or pathv[2] or pathv[3]
    
    if 8 in faultv:
        pathv[5]=1
    else:
        pathv[5]=pathv[2]
    
    if 6 in faultv:
        pathv[6]=1
    else:
        pathv[6]=pathv[4]
    
    if 7 in faultv:
        inpv[4]=0
    
    if 13 in faultv:
        pathv[7]=1
    else:
        pathv[7]=pathv[6]
    
    if 16 in faultv:
        pathv[8]=1
    else:
        pathv[8]=pathv[6]
    
    if 9 in faultv:
        pathv[9]=1
    else:
        pathv[9]=pathv[6] or pathv[5] or pathv[3]
    
    if 14 in faultv:
        pathv[10]=1
    else:
        pathv[10]=pathv[7]
    
    if 17 in faultv:
        pathv[11]=1
    else:
        pathv[11]=pathv[8]
    
    if 10 in faultv:
        pathv[12]=1
    else:
        pathv[12]=pathv[9] or (int)(not inpv[4])
    
    if 15 in faultv:
        pathv[13]=1
    else:
        pathv[13]=pathv[10]
    
    if 18 in faultv:
        pathv[14]=1
    else:
        pathv[14]=pathv[11]
    
    if 11 in faultv:
        pathv[15]=1
    else:
        pathv[15]=pathv[12]
    
    if 12 in faultv:
        pathv[16]=1
    else:
        pathv[16]=pathv[15]
    
    if 20 in faultv:
        pathv[17]=0
    else:
        pathv[17]=(int)(not pathv[16])
    
    if 19 in faultv:
        pathv[18]=0
    else:
        pathv[18]=(int)(not pathv[16])
    
    if 21 in faultv:
        pathv[19]=1
    else:
        pathv[19]=(int)(not pathv[17])
    
    if 22 in faultv:
        pathv[20]=1
    else:
        pathv[20]=pathv[19]
    
    if 23 in faultv:
        pathv[21]=1
    else:
        pathv[21]=pathv[14] or pathv[15] or pathv[20]
    
    if 24 in faultv:
        pathv[22]=0
    else:
        pathv[22]=(int)(not(pathv[21] or pathv[16]))
    
    if 25 in faultv:
        pathv[23]=1
    else:
        pathv[23]=pathv[14] or pathv[20]
    
    if 26 in faultv:
        pathv[24]=0
    else:
        pathv[24]=(int)(not(pathv[14] or pathv[20]))
        
    if 27 in faultv:
        pathv[25]=1
    else:
        pathv[25]=(int)(not(pathv[24]))
    
    outv[0]=pathv[13] and pathv[21]
    outv[1]=pathv[14]
    outv[2]=pathv[14] and pathv[21]
    outv[3]=pathv[14] and pathv[21]
    outv[4]=pathv[23]
    outv[5]=pathv[25]
    outv[6]=(int)(not pathv[22])
    outv[7]=(int)(not pathv[22])
    outv[8]=(int)(not pathv[18])