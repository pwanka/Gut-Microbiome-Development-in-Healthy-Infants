import pandas as pd
import matplotlib.pyplot as plt 

data = "Data/PubMed_Timeline_Results_by_Year.csv"    # Downloaded Dec 1st 2022 10:45
df = pd.read_csv(data, header = 1, index_col= "Year").iloc[1:] # read file and remove 2023

# Calculate percentage of publications published from that year on
df["PercentSince"] = (df.cumsum()/df.sum())*100
df["PercentSince"] = df["PercentSince"].round(2)

# What percentage was published since 2017
print("{}% of all publications were published since 2017".format(df.loc[2017,"PercentSince"]))

# Plot publications per year 
plot_df = df.loc[range(1980,2023),"Count"]
plot_df.plot.bar() 
plt.xlabel("Year")
plt.ylabel("Number of Publications")
plt.xticks(range(0,len(plot_df),2),range(1980,2023,2))
plt.tight_layout()
figure = plt.gcf()
figure.set_size_inches(7,4)
plt.savefig("1.1_publications_human_microbiome.png")
plt.show()
