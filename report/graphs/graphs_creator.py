import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.DataFrame.from_csv('core_task.csv', sep=';')
#data.tail()
print('Data loaded.')


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
    ax.semilogx(data.index, [data[data.index==i][element + " err"].mean() for i in data.index.values], label='Error mean')
    ax.semilogx(data.index,[data[data.index==i][element + " err"].var() for i in data.index.values], label='Error variance')
    ax.grid()
    plt.title('Analysis of ' + element + " error for the " + mytype + ' option')
    plt.xlabel('Monte-Carlo Iterations')
    plt.ylabel('value')
    plt.legend()
    plt.savefig(mytype + element.replace(' ', '') + '.png', dpi=my_dpi)
    print('Created: ' + mytype + element.replace(' ', '') + '.png')
    plt.close()

for element in greeks_names:
    data[str(element + " err")] = np.abs(data[element] - data[str(element.rsplit(' ', 1)[0] + " th")])
    my_dpi = 72
    plt.figure(figsize=(1000/my_dpi, 600/my_dpi), dpi=my_dpi)
    fig, ax = plt.subplots()
    runtime = [data[data.index==i][element + " time"].mean() for i in data.index.values]
    ax.semilogx(runtime, [data[data.index==i][element + " err"].mean() for i in data.index.values], label='Error mean')
    ax.semilogx(runtime,[data[data.index==i][element + " err"].var() for i in data.index.values], label='Error variance')
    ax.grid()
    plt.title('Analysis of ' + element + " error for the " + mytype + ' option')
    plt.xlabel('Computing time')
    plt.ylabel('value')
    plt.legend()
    plt.savefig(mytype + element.replace(' ', '') + 'time.png', dpi=my_dpi)
    print('Created: ' + mytype + element.replace(' ', '') + 'time.png')
    plt.close()

for element in greeks_names:
    print(element + ' error mean: ' + str(data[data.index==100000][element + " err"].mean()))

for element in greeks_names:
    print(element + ' error var: ' + str(data[data.index==100000][element + " err"].var()))

for element in greeks_names:
    print(element + ' time: ' + str(data[data.index==100000][element + " time"].mean()))
