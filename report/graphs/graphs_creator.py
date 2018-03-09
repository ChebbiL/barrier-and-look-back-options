import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.DataFrame.from_csv('core_task.csv', sep=';')
data.tail()



data["delta lr err"] = np.abs(data["delta lr"] - data["delta th"])


mytype = "euro"
greeks = [('delta',['lr', 'pw']), ('gamma',['pwlr', 'lrpw', 'lrlr']), ('vega',['lr', 'pw'])]

greeks_files_temp = [[greek[0]+'_'+greek[1][j] for j in range(len(greek[1]))] for greek in greeks]
greeks_files = []
for element in greeks_files_temp:
    greeks_files += element

greeks_names_temp = [[greek[0]+' '+greek[1][j] for j in range(len(greek[1]))] for greek in greeks]
greeks_names = []
for element in greeks_names_temp:
    greeks_names += element


for element in greeks_names:
    data[str(element + " err")] = np.abs(data[element] - data[str(element.rsplit(' ', 1)[0] + " th")])
    my_dpi = 72
    plt.figure(figsize=(1000/my_dpi, 600/my_dpi), dpi=my_dpi)
    fig, ax = plt.subplots()
    ax.semilogx([data[data.index==i][element + " err"].mean() for i in data.index.values], label='Error mean')
    ax.semilogx([data[data.index==i][element + " err"].var() for i in data.index.values], label='Error variance')
    ax.grid()
    plt.title('Analysis of ' + element + " error for the " + mytype + ' option')
    plt.xlabel('log(iterations)')
    plt.ylabel('value')
    plt.legend()
    plt.savefig(mytype + element + '.png', dpi=my_dpi)
