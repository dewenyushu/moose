import numpy as np
import matplotlib.pyplot as plt

# Data for processor counts and execution times for three numerical experiments
processors = np.array([1, 2, 4, 8])
pre_no_reset = np.array([11971.215, 3708.090, 1750.639, 1860.622])
pre_reset = np.array([5489.454, 3551.250, 3007.742, 1787.323])
hash = np.array([3355.417, 2126.166, 1556.862, 1289.094])

# Calculate the ideal execution times based on the first experiment's single processor time
ideal_time = pre_no_reset[0] / processors

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(processors, pre_no_reset, marker='o', label='Preallocation, no reset')
plt.plot(processors, pre_reset, marker='s', label='Preallocation, reset')
plt.plot(processors, hash, marker='^', label='Hash')
plt.plot(processors, ideal_time, linestyle='--', color='k', label='Ideal')

# Adding labels and title
plt.xlabel('Number of Processors')
plt.ylabel('Execution Time (s)')
plt.title('Scaling')
plt.legend()
plt.grid(True)
plt.xscale('log', base=2)
plt.yscale('log', base=10)

# Save the plot to a PNG file
plt.savefig('scaling_study.png')
