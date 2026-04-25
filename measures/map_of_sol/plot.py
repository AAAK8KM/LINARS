import numpy as np 
import matplotlib.pyplot as plt


A=np.array([
    [10,3,6],
    [3,5,1],
    [6,1,8]
])

b=np.array([
    [1],
    [2],
    [1]
])

x=np.linalg.inv(A) @ b

print(x)

e,v=np.linalg.eig(A)

vx=v[0]
vn=v[1]

print(x,vx)
mx=(x.T@vx)[0]
mn=(x.T@vn)[0]
r=x.reshape(3,)-vx*mx-mn*vn

print(r)

delta = 0.001
xv = np.arange(-0.05, mx+0.1, delta)
yv = np.arange(mn-0.1, 0.05, delta)
X, Y = np.meshgrid(xv, yv)


# AI generated thing for fast point computation 
V = np.stack((vx[:, None, None] * X + vn[:, None, None] * Y + r[:, None, None]*Y/Y), axis=0)
Z = (0.5 * np.einsum('ikl, ikl -> kl', V, AV) 
     - np.einsum('k,kij->ij', b.squeeze(), V))

fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z)
ax.clabel(CS, fontsize=10)
ax.set_title('Solution finding')

import numpy as np

def load_vectors(filepath):
    # np.loadtxt automatically ignores blank lines and extra whitespace
    flat_data = np.loadtxt(filepath)
    
    if flat_data.size % 3 != 0:
        raise ValueError("File must contain groups of exactly 3 numbers.")
        
    return flat_data.reshape(-1, 3)

print(load_vectors("CGDV.csv"))

cg=load_vectors("CGDV.csv")
ax.plot(vx.T@cg.T,vn.T@cg.T, '--b', label="CGD")

sg=load_vectors("SGDV.csv")
print(cg.shape,sg.shape)
ax.plot(vx.T@sg.T,vn.T@sg.T, '-.r', label="SteepreGD")

s=load_vectors("SimpleV.csv")
ax.plot(vx.T@s.T,vn.T@s.T, ':g', label="Simple")

ax.legend()

ax.set_xlabel("$\lambda_{max}$ coordinate")

ax.set_ylabel("$\lambda_{min}$ coordinate")

fig.savefig("1.png")

#print(e,v)
#print(v[0].T@v[2])

#matplotlib.pypylot.countour
#пстроить по плоскотb на lambdamin, lamdamax проекции шагов