import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


c = 1.0

n = 100
dx = 0
dt = 1e-3
t_max = 0.5

x = None
u = None
u_new = None
nu = 0.0

if __name__ == "__main__":

    # 初期化
    dx = 1/float(n)
    x = np.linspace(0, 1, n+1)
    u = np.array([0.0]*(n+1))
    u_new = np.array([0.0]*(n+1))
    nu = c*dt*float(n)

    fig = plt.figure()
    ims = []

    #初期条件
    for i, x_i in enumerate(x):
        if x_i <= 0.5 and x_i > 0:
            u[i] = 1

    print("Courant number :{}".format(nu))

    for j in range(int(np.round(t_max / dt))):
        #差分方程式を計算
        for i in range(1,n):
            u_new[i] = u[i]-1/2.0*nu*(u[i+1]-u[i-1])
        
        #境界条件を適用
        u_new[0] = 0
        u_new[-1] = 0


        im = plt.plot(x,u,color='blue')
        ims.append(im)

        #更新
        u = u_new
    
    print("caluculation finished.")

    ani =  animation.ArtistAnimation(fig, ims, interval=33)
    ani.save("Advection1D\\output.gif", writer="imagemagick")

