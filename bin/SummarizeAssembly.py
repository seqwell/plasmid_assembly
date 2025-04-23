#!/usr/bin/env python3
import glob
import sys
import pandas as pd
import os
from datetime import datetime
from matplotlib import pyplot as plt


now = datetime.now()
timestamp = now.strftime("%Y-%m-%d-%H-%M")

out_name = f"assembly-report-{timestamp}"


files=glob.glob('*.info.csv')

info=pd.DataFrame(columns=['ID', 'Length', 'Circle', 'Ambiguity'])
for f in sorted(files):

    filename=f.split('.')[0]

    try:
        df1=pd.read_csv(f)
 
        df2=pd.read_csv(filename+'.depth.csv', sep='\t',names=['', 'position','coverage'])
        plt.plot(df2.position, df2.coverage)
        plt.title(filename)
        plt.savefig(filename+'.png', dpi=200)
        plt.close()
    
        df3=pd.read_csv(filename + '.count.csv', sep='\t', names=['count'])

        df1[['# Sequences']]=df3['count']
        df1[['MEAN_coverage']]=df2.coverage.mean()
        df1[['MEDIAN_coverage']]=df2.coverage.median()
        df1[['SD_coverage']]=df2.coverage.std()
        df1[['PCT_50X']]=sum(df2.coverage>50)/len(df2)
        df1[['PCT_75X']]=sum(df2.coverage>75)/len(df2)
        df1[['PCT_100X']]=sum(df2.coverage>100)/len(df2)
        df1[['PCT_125X']]=sum(df2.coverage>125)/len(df2)
    
    except:
        columns = ['ID','Length','Circle',
                    'Ambiguity','# Sequences','MEAN_coverage',
                    'MEDIAN_coverage','SD_coverage',
                    'PCT_50X','PCT_75X','PCT_100X','PCT_125X']

        data = [filename] + ['']*11
        df1 = pd.DataFrame([data], columns = columns)

    info= info.append(df1)
info.to_csv(out_name + '.csv', index=False)


