import matplotlib.pyplot as plt
import numpy as np
import csv
import os, os.path

x = []
y = []
DIR = '/tmp'
i =  len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])
print i



firstline = True
with open('example.csv','r') as csvfile:

    plots = csv.reader(csvfile, delimiter=',')

    # skip first line
    next(csvfile)

    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))


x1 = -1*np.array(x)
plt.plot(x,y,x1,y, marker=".",linestyle='None')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement')
plt.legend()

plt.xlim([-40,40])
plt.ylim([-140,77])
plt.gca().set_aspect('equal', adjustable='box')
plt.show()



#f = open("/home/stephen/Documents/Meshless/build/bin/preform/History/Strain/strain_1.txt", "r")
#print(f.read()) 



