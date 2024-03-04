import pandas as pd
import re
import sys

def process_motif_data(data, motif_column_index, filter_motifs):
    motif_column = data.columns[motif_column_index]
    motif = motif_column.split(" ")[0]
    print(motif)
    df_split = data[motif_column].str.split('\),', expand=True)
    df_split['original_row'] = data.index

    df_expanded = df_split.melt(id_vars=['original_row'], value_vars=list(range(df_split.shape[1] - 1)))
    df_expanded = df_expanded.dropna()[['original_row', 'value']]
    df_result = pd.merge(data, df_expanded, left_index=True, right_on='original_row')
    df_result = df_result.drop(motif_column, axis=1).rename(columns={'value': motif_column})
    df_result = df_result.reset_index(drop=True)



    df_split = df_result[motif_column].str.split('[(),]', expand=True)
    df_split.columns = [f'{motif}_Numeric', f'{motif}_Sequence', f'{motif}_Strand', f'{motif}_Value', f'{motif}_empty']
    df_split = df_split.drop([f'{motif}_Value', f'{motif}_empty'], axis=1)

    df_result = df_result.drop(motif_column, axis=1)
    df_result = pd.concat([df_result, df_split], axis=1)
    df_result = df_result.loc[~df_result[f'{motif}_Sequence'].isin(filter_motifs)]

    # Reset the index if needed
    df_result = df_result.reset_index(drop=True)    
    
    return df_result



args = sys.argv[1:]

if len(args) != 1:
    raise ValueError("Usage: python your_script.py <input_file>")

input_file = args[0]

name = input_file[:-28]
print(name)
data = pd.read_csv(input_file, sep='\t')
data = process_motif_data(data, 21, ["AAAAAA", "TTTTTT"])
data = process_motif_data(data, 21, ["TTTTTT", "AAAAAA"])
data = process_motif_data(data, 22, [])
#data = process_motif_data(data, 23, [])

for col in data.columns:
    print(col)



data["AAUAAA_Numeric"] = pd.to_numeric(data["AAUAAA_Numeric"])
data["UUGUUU_Numeric"] = pd.to_numeric(data["UUGUUU_Numeric"])
data["TGTA_Numeric"] = pd.to_numeric(data["TGTA_Numeric"])
#data["YA_Numeric"] = pd.to_numeric(data["YA_Numeric"])

print("4")

positive = data[(data["Strand"] == "+") & (131 > (data["AAUAAA_Numeric"] - data["TGTA_Numeric"])) & ((data["AAUAAA_Numeric"] - data["TGTA_Numeric"]) > 59) & (28 > (data["UUGUUU_Numeric"] - data["AAUAAA_Numeric"])) & ((data["UUGUUU_Numeric"] - data["AAUAAA_Numeric"]) > 5) ]
negative = data[(data["Strand"] == "-") & (131 > (data["TGTA_Numeric"] - data["AAUAAA_Numeric"])) & ((data["TGTA_Numeric"] - data["AAUAAA_Numeric"]) > 59) & (28 > (data["AAUAAA_Numeric"] - data["UUGUUU_Numeric"])) & ((data["AAUAAA_Numeric"] - data["UUGUUU_Numeric"]) > 5) ]

perfekte = pd.concat([positive, negative])
positionen_perfekt = perfekte[["Chr", "Start", "End", "Strand"]].drop_duplicates()
print("5")
positionen_perfekt.to_csv(f"{name}_perfect_motifs.bed", sep='\t', header=False, index=False)
