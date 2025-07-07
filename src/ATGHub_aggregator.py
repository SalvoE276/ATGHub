import csv
from gprofiler import GProfiler



def get_unique_genes(csvfile, delimiter):
    reader = csv.reader(csvfile, delimiter=delimiter)
    gene_list = []
    for gene in reader:
        gene_list.append(gene[0])

    gene_list = gene_list[1:] #remove header
    gene_list = list(set(gene_list)) # remove repeated genes
    return gene_list


# load AutophagyNet core genes
with open(r"databases_raw_data/autophagynet_810933fbcd1232140ef7_csv.csv", "r") as csvfile:
    core_autophagynet_genes = get_unique_genes(csvfile, ",")

# load AutophagyNet genes
with open(r"databases_raw_data/autophagynet_05e20118eb7667423678_csv.csv", "r") as csvfile:
    autophagynet_genes = get_unique_genes(csvfile, ",")

# load HADB genes
with open(r"databases_raw_data/HADB.csv", "r") as csvfile:
    HADB_genes = get_unique_genes(csvfile, ",")

# load HAMdb gene only human
with open(r"databases_raw_data/protein-role_filter.csv", "r") as csvfile:
    HAMdb_genes = get_unique_genes(csvfile, ";")

# load Gene Ontology genes GO:0061919
with open(r"databases_raw_data/GO_0061919.tsv", "r") as csvfile:
    GO_genes = get_unique_genes(csvfile, "\t")


# get unique gene list
final_gene_list = autophagynet_genes.copy()
final_gene_list.extend(HADB_genes) # merge HADB_genes
final_gene_list.extend(HAMdb_genes) # merge HAMdb_genes
final_gene_list.extend(GO_genes) # merge GO:0061919 genes
final_gene_list = list(set(final_gene_list)) # remove repeated genes
final_gene_list.sort()



#g:profier symbol conversion
gp = GProfiler(return_dataframe=False)
entrez_id_df = gp.convert(organism='hsapiens', query=final_gene_list, target_namespace='ENTREZGENE_ACC')
esnsembl_id_df = gp.convert(organism='hsapiens', query=final_gene_list, target_namespace='ENSG')




# search recurrences in all databases
gene_dataset = [["Gene Symbol", "Entrez ID", "Ensembl ID", "Core AutophagyNet", "AutophagyNet", "HADB", "HAMdb", "GO-0061919"]]
for e in final_gene_list:
    gene_vector = [e]
    
    for gene_dict in entrez_id_df:
        if gene_dict['incoming'] == e:
            gene_vector.append(gene_dict['converted'])
            break
    if len(gene_vector) < 2:
        gene_vector.append('')
    
    for gene_dict in esnsembl_id_df:
        if gene_dict['incoming'] == e:
            gene_vector.append(gene_dict['converted'])
            break
    if len(gene_vector) < 3:
        gene_vector.append('')
    
    
    if e in core_autophagynet_genes:
        gene_vector.append(1)
    else:
        gene_vector.append(0)
    
    
    if e in autophagynet_genes:
        gene_vector.append(1)
    else:
        gene_vector.append(0)

    if e in HADB_genes:
        gene_vector.append(1)
    else:
        gene_vector.append(0)

    if e in HAMdb_genes:
        gene_vector.append(1)
    else:
        gene_vector.append(0)

    if e in GO_genes:
        gene_vector.append(1)
    else:
        gene_vector.append(0)
    
    gene_dataset.append(gene_vector)



# export to csv
with open(r"ATGHub_dataset/ATGHub.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(gene_dataset)