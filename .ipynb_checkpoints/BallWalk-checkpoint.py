import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
def oracle(x):
    return ((np.sum(x) <= 2) and np.sum(x < 0)) == 0 
X = np.array([0, 0]).reshape(2,1)
print(X)
calls = 0
pos = X
out_pos = X

delta = .989
for j in range(1,100):
    for i in range(1, 100):
        tmpx = X + np.random.uniform(0, delta, (1,2))
        if(oracle(tmpx)):
            calls = calls + 1
            X = tmpx
            print(X)
            pos = np.hstack((pos, X))
            plt.plot(pos[0, :], pos[1,:])
    X = np.array([0, 0]).reshape(2,1)
    pos = X
    out_pos = X

plt.show()
print(calls)

print(pos)
print(out_pos)
