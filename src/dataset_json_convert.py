import csv, json


#read csv ATGHub
with open(r"ATGHub_dataset/ATGHub.csv", "r") as csvfile:
    reader = csv.reader(csvfile, delimiter=";")
    ATGHub = []
    for gene in reader:
        ATGHub.append(gene)

    ATGHub = ATGHub[1:] #remove header

#from csv to json list of dictionaries
json_output = []
for e in ATGHub:
    json_output.append({"Gene Symbol" : e[0],
                        "Entrez ID" : e[1],
                        "Ensembl ID" : e[2],
                        "Core AutophagyNet" : e[3],
                        "AutophagyNet" : e[4],
                        "HADB" : e[5],
                        "HAMdb" : e[6],
                        "GO-0061919" : e[7]})


#save json file
with open(r"ATGHub_dataset/ATGHub.json", "w") as jsonfile:
    json.dump(json_output, jsonfile)