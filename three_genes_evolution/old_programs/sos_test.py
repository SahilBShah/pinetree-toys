import pinetree as pt
import pandas
import random
import numpy as np
import string
import filecmp
import math
import os


#Removes unnecessary rows and columns in produced file
def edit_new_file(new_file):

  #new_file = new_file.drop("time", axis=1)
  new_file = new_file.drop(columns="protein", axis=1)
  new_file = new_file.drop(columns="ribo_density", axis=1)
  new_file = new_file[new_file.species != '__proteinX_rbs']
  new_file = new_file[new_file.species != '__proteinY_rbs']
  new_file = new_file[new_file.species != '__proteinZ_rbs']
  new_file = new_file[new_file.species != '__ribosome']
  new_file = new_file[new_file.species != '__rnase_site']
  new_file = new_file[new_file.species != '__rnase_site_ext']
  new_file = new_file[new_file.species != 'rnapol']
  new_file = new_file[new_file.species != 'p1']
  new_file = new_file[new_file.species != 'p2']
  return new_file

#Removes unnecessary rows and columns in target file
def edit_target_file(target_file, name_of_file):

  if name_of_file == "random_data.tsv":
      #target_file = target_file.drop("time", axis=1)
      return target_file
  else:
      #target_file = target_file.drop("time", axis=1)
      target_file = target_file.drop(columns="protein", axis=1)
      target_file = target_file.drop(columns="ribo_density", axis=1)
      target_file = target_file[target_file.species != '__proteinX_rbs']
      target_file = target_file[target_file.species != '__proteinY_rbs']
      target_file = target_file[target_file.species != '__proteinZ_rbs']
      target_file = target_file[target_file.species != '__ribosome']
      target_file = target_file[target_file.species != '__rnase_site']
      target_file = target_file[target_file.species != '__rnase_site_ext']
      target_file = target_file[target_file.species != 'rnapol']
      target_file = target_file[target_file.species != 'p1']
      target_file = target_file[target_file.species != 'p2']
      return target_file

def sum_of_squares(target_file):
#row[1][1] = species name
#row[1][2] = transcript
#row[1][0] = time
  sos = 0.0
  target_file = pandas.DataFrame(target_file)
  rows=1
  target_row = target_file.iloc[rows]
  print("THE ROW IS: " + str(target_file.iloc[rows][0]))
  rows+=1
  print("THE ROW IS: " + str(target_file.iloc[rows][1]))


df_name = input("Enter tsv file name: ")
df = pandas.read_table(df_name, delim_whitespace=True, header=0)
  #df = pandas.read_table("three_genes_test_file2.tsv", delim_whitespace=True, header=0)
df = edit_target_file(df, df_name)
sum_of_squares(df)
