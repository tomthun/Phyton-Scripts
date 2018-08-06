# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 10:22:37 2018

@author: heinzinger
"""
import argparse
from perseuspy import pd

parser = argparse.ArgumentParser(description='Parser Template')
# This is how one can add user based input parameters
parser.add_argument('--n_numbers',type = int, help = 'Parameter number one')
# Input and Output parameters are a must have
parser.add_argument('input',help='path to input file')
parser.add_argument('output',help='path to output file')

arg = parser.parse_args()

df = pd.read_perseus(arg.input)
# df_short = df.head(n_numbers) --> do operation with your dataframe here 
#                          e.g head() for top n rows
#                          OUTPUT MUST BE DATAFRAME
df_short.to_perseus(arg.output)
