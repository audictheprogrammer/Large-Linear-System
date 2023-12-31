import matplotlib.pyplot as plt
import math
import sys

input = "Graph1.txt"
output = "Graph1.png"
title = "Error analysis with MAX_R = 15 and Mat_N = 300 (log_2 Scale)"

if len(sys.argv) == 2:
    if sys.argv[1] == "2":
        input = "Graph2.txt"
        output = "Graph2.png"
        title = "Etude d'erreur avec MAX_R = 15 et Mat_N = 300 en log_2"
        print("Plotting the Graph 2.")
    else:
        print("Plotting the Graph 1.")
else:
    print("Plotting the Graph 1.")


with open(input, "r") as f:
    X = []
    Y = []
    for line in f:
        x, y = map(float, line.split())
        x = math.log(x, 2)
        y = math.log(y, 2)
        X.append(x)
        Y.append(y)

print("X = ", X)
print("Y = ", Y)

plt.plot(X, Y, label="Log_2 || B-Br ||")
plt.title(title)
plt.xlabel("r")
plt.ylabel("Erreur")
plt.legend()
plt.savefig(output)
plt.clf()
