###Ignore this file, none of this code is used#### 
##Playground file for testig pre processing steps##


import pandas as pd

df_CDH1 = pd.read_csv('data/CDH1.csv')
df_GATA3 = pd.read_csv('data/GATA3.csv')
df_PIK3CA = pd.read_csv('data/PIK3CA.csv')
df_TP53 = pd.read_csv('data/TP53.csv')

# print(df_CDH1)
# print("******************************")
# print(df_GATA3)
# print("******************************")
# print(df_PIK3CA)
# print("*********/*********************")
# print(df_TP53)


labl = "Mutation Type"
toSelect = df_CDH1.columns.to_list()
toSelect.remove(labl)
toSelect.remove("Reference")

# df2 = df_CDH1.loc[:, df_CDH1.columns != labl and df_CDH1.columns != "Reference"]
df2 = df_CDH1[toSelect]

# df2 = df_CDH1[labl]

# df2 = df_CDH1.loc[df_CDH1["Mutation Type"] != "Missense_Mutation"]

# st = len(df_CDH1["Reference"][0])

print(df2)



