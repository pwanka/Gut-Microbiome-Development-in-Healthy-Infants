import pandas as pd 

data = pd.read_csv("Data/Sample_Metadata_combine.csv", sep="\t")
labeleddata = data[data["Label"].notna()]

#------ EMMA ------

regex_emma = "01\.[12]\.[SV]\.[0-9]{3}\.c\.[0-9]{3}"

emma_project = data[data["project"] == "EMMA"]
emma_regex = labeleddata[labeleddata["Label"].str.contains(regex_emma)]

index_list = list(set(emma_project.index) | set(emma_regex.index))
all_emma = data.iloc[index_list]

for index in index_list:
    data.loc[index,"project"] = "EMMA"

all_emma.loc[all_emma[all_emma["grouping"] == "15012S021c002"].index,"grouping"] = "012S021c002"
all_emma.loc[all_emma[all_emma["grouping"] == "15012V022c028"].index,"grouping"] = "012V022c028"
all_emma.loc[all_emma[all_emma["grouping"] == "15012V033c001"].index,"grouping"] = "012V033c001"
all_emma.loc[all_emma[all_emma["X.SampleID"] == "20012V025c028"].index,"grouping"] = "012V015c028"
all_emma.loc[all_emma[all_emma["grouping"] == "12V018c028"].index,"grouping"] = "012V018c028"

all_emma["Study_id"] = all_emma["grouping"]

#----- PRIMAL ------
#regex_primal = "[0-9]{2}\.[12]\.[78]\.[0-9]{3}\.[cm]\.[0-9]{3}"

primal_project = data[data["Study_Sample"]==True]

for index in list(primal_project.index):
    data.loc[index,"project"] = "PRIMAL"

primal_project = primal_project[primal_project["exclusion"]!=True]

primal_emma = pd.concat([all_emma,primal_project],ignore_index=True)

reads = data[["X.SampleID","input","filtered","denoised","merged","tabled","nonchim","NOarchaeaEukaryota"]].set_index("X.SampleID")

uniques = primal_emma[primal_emma.duplicated(subset=["Study_id"],keep=False)].sort_values(by=["Study_id"])[["X.SampleID","Study_id"]]
uniques = uniques.set_index("X.SampleID")
reads_of_interest = uniques.join(reads)

tokeep = reads_of_interest.sort_values(by="input").drop_duplicates(subset = ["Study_id"],keep="last").index
ix = list(reads_of_interest.drop(list(tokeep)).index)

data = data.set_index("X.SampleID")
for xid in ix:
    data.loc[xid,"exclusion"] = True

data.to_csv("Metadata_combined_adjusted.csv", sep = "\t")
