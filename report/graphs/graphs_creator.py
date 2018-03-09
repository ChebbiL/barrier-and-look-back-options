import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.DataFrame.from_csv('core_task.csv', sep=';')
data.tail()



data["delta lr err"] = np.abs(data["delta lr"] - data["delta th"])

fig, ax = plt.subplots()
ax.semilogx([data[data.index==i]["delta lr err"].mean() for i in data.index.values], label='Error mean')
ax.semilogx([data[data.index==i]["delta lr err"].var() for i in data.index.values], label='Error variance')
ax.grid()
plt.title('Analysis of delta pw error for the European call option')
plt.xlabel('log(iterations)')
plt.ylabel('value')
plt.legend()
plt.show()
