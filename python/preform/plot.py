import matplotlib.pyplot as plt
import numpy as np
import csv


x = []
y = []

firstline = True
with open('example.csv','r') as csvfile:

    plots = csv.reader(csvfile, delimiter=',')

    # skip first line
    next(csvfile)

    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))



plt.plot(x,y, marker=".",linestyle='None')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement')
plt.legend()

plt.xlim([0,40])
plt.ylim([-140,77])
plt.gca().set_aspect('equal', adjustable='box')


plt.plot(x,y, marker=".",linestyle='None')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement')
plt.legend()

plt.xlim([0,40])
plt.ylim([-140,77])
plt.gca().set_aspect('equal', adjustable='box')
plt.show()


f = open("/home/stephen/Documents/Meshless/build/bin/preform/History/Strain/strain_1.txt", "r")
print(f.read()) 



