'''
Created on 17 dic. 2020

@author: jfsanchezherrero
'''
import pandas as pd
import numpy as np
from collections import defaultdict 

annot_table = "./pruebas/df.csv"
sort_csv = "./pruebas/filtered_results.csv"

annot_table = pd.read_csv(annot_table)
sort_csv = pd.read_csv(sort_csv)

qseqid = list(sort_csv["qseqid"])
sseqid =list(sort_csv["sseqid"])
qseqid.extend(sseqid)
# prot_id = np.sort(qseqid)
prot_id = set(qseqid)
print(prot_id)
 
 
 
# annot_table_filtered = annot_table[~prot_id.duplicated]
# print(annot_table_filtered)
filtered_annot = annot_table[annot_table.index.isin(prot_id)]
dup_annot = "./pruebas/dup_annot.csv"
print(filtered_annot)
filtered_annot.to_csv(dup_annot, header=True, index=False)

## 1a ronda
relations_df = defaultdict(list) 
for index, row in sort_csv.iterrows():
	relations_df[row['qseqid']].append(row['sseqid'])

print (relations_df)

## 2a ronda
new_relations_df = defaultdict(list)
dups=0
for key, value in relations_df.items():
	print ()
	print ("key: " + key)
	print ("value: " + str(value))

	stop=False
	for dup_id, new_value in new_relations_df.items():
		if key in new_value:
			stop=True
			print ("Belongs to group: " + dup_id)

	if not stop:
		for key2, value2 in relations_df.items():
			if (key == key2):
				continue
			else:
				if (key2 in value): 
					for i in value2: 
						if i not in value: 
							value.extend(i)

		dups += 1
		value.extend(key)
		new_relations_df["dup_"+str(dups)] = value
		print(new_relations_df)
		print("**")

print ()
print (new_relations_df)

filtered_annot.loc["duplicate_id"] = filtered_annot.locus_tag.map(new_relations_df)




## Filtrar tabla de anotacion con genes duplicados
## añadir columna con dup_id a tabla de anotacion.
## FIN





