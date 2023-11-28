import pandas as pd

df_CDH1 = pd.read_csv('data/CDH1.csv')
df_GATA3 = pd.read_csv('data/GATA3.csv')
df_PIK3CA = pd.read_csv('data/PIK3CA.csv')
df_TP53 = pd.read_csv('data/TP53.csv')

print(df_CDH1)
print("******************************")
print(df_GATA3)
print("******************************")
print(df_PIK3CA)
print("*********/*********************")
print(df_TP53)

# df2 = pd.DataFrame(df_CDH1["Mutation Type"] != "Missense_Mutation")
# df2 = df_CDH1.loc[df_CDH1["Mutation Type"] != "Missense_Mutation"]

# st = len(df_CDH1["Reference"][0])

# print(st)