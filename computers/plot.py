import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

hardware = pd.read_csv('hardware.csv', sep=';')
print(hardware)
print (hardware.dtypes)

# Create one dataframe for memory (RAM + HD)
RAM = hardware.copy()
del RAM['Pixels']
del RAM['HD']
RAM.columns=['Year', 'value']
RAM['type'] = 'RAM'

HD = hardware.copy()
del HD['Pixels']
del HD['RAM']
# Remove unknown row
HD=HD[HD.Year != 2008]
HD.columns=['Year', 'value']
HD['type'] = 'HD'

mem = pd.concat ([RAM, HD])

g = sns.lineplot(x="Year", y="value", hue="type", data=mem, legend="full")
g.set_yscale('log')
g.set_title("Memory in Megabytes (log scale)")
g.figure.savefig("hardware_mem.png")

plt.clf()
g2 = sns.lineplot(x="Year", y="Pixels", data=hardware)
g2.set_title("Supported screen sizes in Megapixels")
# g2.set_yscale('log')
g2.figure.savefig("hardware_pix.png")



# Pixels = hardware.copy()
# del Pixels['RAM']
# del Pixels['HD']
# Pixels.columns=['Year', 'value']
# Pixels['type'] = 'Pixels'

# RAM = hardware.copy()
# del RAM['Pixels']
# del RAM['HD']
# RAM.columns=['Year', 'value']
# RAM['type'] = 'RAM'

# HD = hardware.copy()
# del HD['Pixels']
# del HD['RAM']
# # Remove unknown row
# HD=HD[HD.Year != 2008]
# HD.columns=['Year', 'value']
# HD['type'] = 'HD'

# # flat = pd.concat ([Pixels, RAM, HD])
# flat = pd.concat ([RAM, HD])
# print(flat)

# g = sns.lineplot(x="Year", y="value", hue="type", data=flat)
# g.set_yscale('log')
# g.figure.savefig("hardware.png")