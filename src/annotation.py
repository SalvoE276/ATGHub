import csv, json


# load ATGHub dataset
with open(r"ATGHub_dataset/ATGHub.json", "r") as jsonfile:
    ATGHub_json = json.load(jsonfile)


# load toolbox2021 annotation table
with open(r"databases_raw_data/toolbox2021.csv", "r") as csvfile:
    reader = csv.reader(csvfile, delimiter=";")
    toolbox_list = []
    for gene in reader:
        toolbox_list.append(gene)

    toolbox_list = toolbox_list[1:] #remove header


ATGHub_annotated = [["Gene Symbol", "Entrez ID", "Ensembl ID", "Core AutophagyNet", "AutophagyNet", "HADB", "HAMdb", "GO-0061919", "Category", "Biological Process"]]
for e in ATGHub_json:
    # filter only genes in all databases
    if bool(int(e["AutophagyNet"])) and bool(int(e["HADB"])) and bool(int(e["HAMdb"])) and bool(int(e["GO-0061919"])):
        e = list(e.values())
        
        # get Category
        categories = []
        for gene in toolbox_list:
            if e[0] == gene[0]:
                categories.append(gene[3])
        e.append(' | '.join(categories))

        # get Biological Function
        functions = []
        for gene in toolbox_list:
            if e[0] == gene[0]:
                functions.append(gene[2])
        e.append(' | '.join(functions))


        ATGHub_annotated.append(e)

# export to csv
with open(r"ATGHub_dataset/ATGHub_annotated.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(ATGHub_annotated)