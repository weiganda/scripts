#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("lcurve.out", names=True)
for name in data.dtype.names[1:-1]:
    plt.plot(data["step"], data[name], label=name)
plt.legend()
plt.xlabel("Step")
plt.ylabel("Loss")
plt.grid()
plt.show()

# Create a logarithmic scale plot
plt.figure(figsize=(10, 5))
for name in data.dtype.names[1:-1]:
    plt.plot(data["step"], data[name], label=name)
plt.yscale("log")
plt.legend()
plt.xlabel("Step")
plt.ylabel("Loss")
plt.title("Logarithmic Scale")
plt.grid()
plt.show()
