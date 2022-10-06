import pandas as pd

#read in files
fungi_WT = pd.read_csv('Fungi_WT_gtf_species_full_list_Evo_Rates.csv')
fungi_WT.fillna(".", inplace = True)
fungi_ortho = pd.read_csv('Fungi_Ortholog_gtf_species_full_list_Evo_Rates.csv')
fungi_ortho.fillna('.', inplace = True)
fungi_novel = pd.read_csv('Fungi_Novel_gtf_species_full_list_Evo_Rates.csv')
fungi_novel.fillna('.', inplace = True)

#Add prefix to each column header
preserved_columns = ['species','commonName','group','database','taxid','link']
	

fungi_WT.rename(columns={col: 'WT_' + col for col in fungi_WT.columns if col not in preserved_columns}, inplace=True)
fungi_ortho.rename(columns={col: 'ortho_' + col for col in fungi_ortho.columns if col not in preserved_columns}, inplace=True)
fungi_novel.rename(columns={col: 'novel_' + col for col in fungi_novel.columns if col not in preserved_columns}, inplace=True)

# Combine dfs, while preserving only preserved_columns [0:6]
result = pd.merge(pd.merge(fungi_WT,fungi_ortho, how = "outer", on = preserved_columns),fungi_novel, how="outer", on = preserved_columns)
# ^^^ Place holder

final_grouping_results = result.groupby('group')[['WT_num_transcript','WT_num_gene','WT_num_exon','WT_num_uniqExon','ortho_num_transcript','ortho_num_gene','ortho_num_exon','ortho_num_uniqExon','novel_num_transcript','novel_num_gene','novel_num_exon','novel_num_uniqExon']].sum()

final_grouping_results["WT_groupEpG"] = final_grouping_results["WT_num_uniqExon"]/final_grouping_results["WT_num_gene"]
final_grouping_results["WT_groupEpT"] = final_grouping_results["WT_num_exon"]/final_grouping_results["WT_num_transcript"]
final_grouping_results["WT_groupTpG"] = final_grouping_results["WT_num_transcript"]/final_grouping_results["WT_num_gene"]
final_grouping_results["ortho_groupEpG"] = final_grouping_results["ortho_num_uniqExon"]/final_grouping_results["ortho_num_gene"]
final_grouping_results["ortho_groupEpT"] = final_grouping_results["ortho_num_exon"]/final_grouping_results["ortho_num_transcript"]
final_grouping_results["ortho_groupTpG"] = final_grouping_results["ortho_num_transcript"]/final_grouping_results["ortho_num_gene"]
final_grouping_results["novel_groupEpG"] = final_grouping_results["novel_num_uniqExon"]/final_grouping_results["novel_num_gene"]
final_grouping_results["novel_groupEpT"] = final_grouping_results["novel_num_exon"]/final_grouping_results["novel_num_transcript"]
final_grouping_results["novel_groupTpG"] = final_grouping_results["novel_num_transcript"]/final_grouping_results["novel_num_gene"]

# Write to a csv
final_grouping_results.to_csv('Fungi_Full_Trait.csv')