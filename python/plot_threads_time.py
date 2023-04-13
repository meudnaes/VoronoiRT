import numpy as np
import matplotlib.pyplot as plt

# Number of workers, 1 thread per core
num_threads = np.array([1,
                        2,
                        3,
                        4,
                        6,
                        8,
                        10])

# time usage in seconds
time_usage = np.array([6175,
                       4600,
                       3675,
                       3130,
                       2655,
                       2434,
                       2330])

fig, ax = plt.subplots()

ax.plot(num_threads,
        1/time_usage*3600)

ax.set_xlabel("# threads")
ax.set_ylabel("speed (1/h)")

plt.show()

# plt.savefig("../img/threads.pdf")
