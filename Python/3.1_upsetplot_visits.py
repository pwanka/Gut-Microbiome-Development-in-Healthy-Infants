import pandas as pd
import upsetplot
import matplotlib.pyplot as plt

df = pd.read_csv("EMMA_Sample_Metadata.csv",sep="\t")
visits = df.set_index("EMMA_ID")[["visit"]]

visit1=list(visits[visits["visit"]=="Visit 1"].index)
visit2=list(visits[visits["visit"]=="Visit 2"].index)
visit3=list(visits[visits["visit"]=="Visit 3"].index)

upsetplotdata = upsetplot.from_contents({"Visit 1":visit1,"Visit 2":visit2,"Visit 3":visit3})
upsetplot.plot(upsetplotdata, sort_by="cardinality", sort_categories_by = None)
plt.tight_layout()
plt.savefig("3.1_upsetplot_visits.png")
plt.show()

print(upsetplotdata.index.value_counts())
print(visits.value_counts())
