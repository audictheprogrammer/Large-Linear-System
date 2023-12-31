import matplotlib.pyplot as plt
import math
import sys

input = "Graph1.txt"
output = "Graph1.png"
title = "Error analysis with MAX_R = 15 and Mat_N = 300 (log_2 Scale)"

""" Examples.
Default: Graph 1.
$ python3 script.py 1
Optional: Graph2.
$ python3 script.py 2
"""

if len(sys.argv) == 2:
    if sys.argv[1] == "2":
        input = "Graph2.txt"
        output = "Graph2.png"
        title = "Error analysis with MAX_R = 300 and Mat_N = 300 (log_2 Scale)"
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
        # We can choose between log_2 and natural log.
        # x = math.log(x)
        # y = math.log(y)
        x = math.log(x, 2)
        y = math.log(y, 2)
        X.append(x)
        Y.append(y)

print("X = ", X)
print("Y = ", Y)

plt.plot(X, Y, label="log(r) -> log_2 || B-Br ||")
plt.title(title)
plt.xlabel("log(r)")
plt.ylabel("log_2 || B-Br ||")
plt.legend()
plt.savefig(output)
plt.clf()
