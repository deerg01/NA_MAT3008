import matplotlib.pyplot as plt
import numpy as np

n = [100, 1000, 10000, 100000]

def uniform_hist():
    global index
    index = 1
    for i in n:
        uniform = open("u_" + str(i) + ".txt")
        uniform_dis = np.array([], np.float64)
        while True:
            line = uniform.readline()
            if not line: break
            uniform_dis = np.append(uniform_dis, np.float64(line.rstrip("\n")));
        plt.subplot(2, len(n), index)
        plt.title("Uniform: " + str(i) + " samples")
        plt.hist(uniform_dis, 100, density=True, color='indigo', rwidth=0.6)
        index += 1

def gauss_hist():
    global index
    for i in n:
        gauss = open("g_" + str(i) + ".txt")
        gauss_dis = np.array([], np.float64)
        while True:
            line = gauss.readline()
            if not line: break
            gauss_dis = np.append(gauss_dis, np.float64(line.rstrip("\n")));
        plt.subplot(2, len(n), index)
        plt.title("Gauss: m" + str(i) + " samples")
        plt.hist(gauss_dis, 100, density=True, color='seagreen', rwidth=0.6)
        # plt.show()
        index += 1

plt.figure(1, figsize=(12, 6))
uniform_hist()
plt.figure(1, figsize=(12, 6))
gauss_hist()
plt.tight_layout()
plt.show()