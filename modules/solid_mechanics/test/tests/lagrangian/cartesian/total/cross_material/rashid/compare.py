import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    gn = np.loadtxt("new_out.csv", delimiter=",", skiprows=1)
    ra = np.loadtxt("old_out.csv", delimiter=",", skiprows=1)

    for i in range(3):
        (l,) = plt.plot(gn[:, 0], gn[:, i + 1], label=f"{i}")
        plt.plot(ra[:, 0], ra[:, i + 1], ls="--", color=l.get_color(), lw = 5, alpha = 0.75)

    plt.show()
