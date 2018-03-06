import matplotlib.pyplot as plt
import pandas as pd

data = pd.DataFrame.from_csv('core_task.csv', sep=';')
data.head()

fig, ax = plt.subplots()
#ax.semilogx(data['delta lr'])
ax.semilogx(data['delta pw'])
ax.grid()
plt.show()
