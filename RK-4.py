from skimage.restoration import denoise_nl_means, estimate_sigma
from scipy.io import loadmat
from pathlib import Path
from cv2 import cv2
import scipy.io
import numpy as np
import skimage.color
import skimage.io
from matplotlib import pyplot as plt
import math
import scipy
from scipy import ndimage
import copy

#For Scenario 1
def f(t,v):
    g=9.81
    m=85
    rho=1.21
    C=1
    A=0.7
    y = g-(C*rho*A*(v**2))/(2*m)
    return y

#For Scenario 1
def g(t,v):
    g1=9.81
    m=85
    rho=1.21
    C=0.7
    A=0.18
    y = g1-(C*rho*A*(v**2))/(2*m)
    return y
#



def A2Q3aScnA():
    h=0.2           #Step-size
    i=0
    t=np.zeros(500)
    v=np.zeros(500)
    err=np.zeros(500)
    t[1]=10
    v[1]=1000
    x=0
    err[x]=100
    # v[x]=0.00000000000000001
    v[x] = 10 ** (-15)
    v_analytical = np.zeros(500)
    err_Num_Formula= np.zeros(500)
    g1 = 9.81
    m = 85
    rho = 1.21
    C = 1
    A = 0.7

    vss = ((2 * m * g1) / (C * rho * A)) ** 0.5
    tau = ((2 * m) / (C * rho * A * g1)) ** 0.5

    # print(tau)
    while( np.abs(err[x]/v[x]) > 10**(-6)   ):
        k1=f(t[i],v[i])
        k2=f(t[i] +h/2,v[i]+ h*(k1)/2  )
        k3=f(t[i] +h/2 ,v[i] + h*(k2)/2 )
        k4=f(t[i] +h , v[i] +h*k3)

        v[i+1]=v[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4 )
        t[i+1]=t[i] + h

        v_analytical[i] = vss * (math.exp(2 * t[i] / tau) - 1) / (math.exp(2 * t[i] / tau) + 1)
        # v_analytical[i]=vss*np.tanh(t[i]/tau)
        err_Num_Formula[i]= np.abs(v_analytical[i]-v[i])

        err[i]=v[i+1]-v[i]
        x=i
        i+=1
    print(f'The velocity is {v[i]}')
    print(f'The time required to reach is {t[i]}')
    print(f'The no. of iterations is {i}')
    v1=[]

    for j in range(1,len(v)):
        if(v[j]!=0):
            v1.append(v[j])
        j+=1

    t1=[]
    for j in range(1,len(t)):
        if(t[j]!=0):
            t1.append(t[j])
        j+=1

    err_Num_Formula = err_Num_Formula[0:len(v1)]
    # print(err_Num_Formula)
    err_Num_Formulax10=(10**(7))*err_Num_Formula

    plt.title('Velocity vs Time Scenario A ')
    plt.plot(t1, v1)
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    plt.show()

    # plt.plot(t1, v1)
    plt.plot(t1, err_Num_Formulax10)
    plt.xlabel('Time')
    plt.ylabel('b part err_Num_Formulax10^(-7))')
    plt.title('Error analytical numerical Scenario A ')
    plt.show()



def A2Q3aScnB():
    h = 0.2  # Step-size
    i = 0
    t = np.zeros(500)
    v = np.zeros(500)
    err = np.zeros(500)
    #t[1] = 10
    v[1] = 1000
    x = 0
    err[x] = 100
    # v[x] = 0.0000000000000001
    v[x]=10**(-15)
    v_analytical = np.zeros(500)
    err_Num_Formula= np.zeros(500)
    g1 = 9.81
    m = 85
    rho = 1.21
    C = 0.7
    A = 0.18

    vss = ((2 * m * g1) / (C * rho * A)) ** 0.5
    tau = ((2 * m) / (C * rho * A * g1)) ** 0.5

    while (np.abs(err[x] / v[x]) > 10 ** (-6)):
        k1 = g(t[i], v[i])
        k2 = g(t[i] + h / 2, v[i] + h * (k1) / 2)
        k3 = g(t[i] + h / 2, v[i] + h * (k2) / 2)
        k4 = g(t[i] + h, v[i] + h * k3)

        v[i + 1] = v[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        t[i + 1] = t[i] + h

        v_analytical[i] = vss * (math.exp(2 * t[i] / tau) - 1) / (math.exp(2 * t[i] / tau) + 1)

        err_Num_Formula[i] = np.abs(v_analytical[i] - v[i])

        err[i] = v[i + 1] - v[i]
        x = i
        i += 1
    print(f'The velocity is {v[i]}')
    print(f'The time required to reach is {t[i]}')
    print(f'The no. of iterations is {i}')

    v1=[]
    for j in range(1, len(v)):
        if (v[j] != 0):
            v1.append(v[j])
        j += 1

    t1 = []
    for j in range(1, len(t)):
        if (t[j] != 0):
            t1.append(t[j])
        j += 1


    err_Num_Formula = err_Num_Formula[0:len(v1)]

    err_Num_Formulax10=(10**(7))*err_Num_Formula

    plt.title('Velocity vs Time Scenario B ')
    plt.plot(t1, v1, )
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    plt.show()

    plt.plot(t1, err_Num_Formulax10)
    plt.xlabel('Time')
    plt.ylabel('b part err_Num_Formulax10^(-7)')
    plt.title('Error analytical numerical Scenario B ')
    plt.show()



def A2Q3_CPART_ScnA(h):
              #Step-size
    i=0
    t=np.zeros(500)
    v=np.zeros(500)
    err=np.zeros(500)
    t[1]=10
    v[1]=1000
    x=0
    err[x]=100
    # v[x]=0.000000000000001
    v[x]= 10**(-15)

    v_analytical = np.zeros(500)
    err_Num_Formula= np.zeros(500)
    g1 = 9.81
    m = 85
    rho = 1.21
    C = 1
    A = 0.7

    vss = ((2 * m * g1) / (C * rho * A)) ** 0.5
    tau = ((2 * m) / (C * rho * A * g1)) ** 0.5


    while( np.abs(err[x]/v[x]) > 10**(-6)   ):
        k1=f(t[i],v[i])
        k2=f(t[i] +h/2,v[i]+ h*(k1)/2  )
        k3=f(t[i] +h/2 ,v[i] + h*(k2)/2 )
        k4=f(t[i] +h , v[i] +h*k3)

        v[i+1]=v[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4 )
        t[i+1]=t[i] + h

        v_analytical[i] = vss * (math.exp(2 * t[i] / tau) - 1) / (math.exp(2 * t[i] / tau) + 1)

        err_Num_Formula[i]= np.abs(v_analytical[i] -v[i])

        err[i]=v[i+1]-v[i]
        x=i
        i+=1
    # print(f'The velocity is {v[i]}')
    # print(f'The time required to reach is {t[i]}')
    # print(f'The no. of iterations is {i}')
    v1=[]

    for j in range(1,len(v)):
        if(v[j]!=0):
            v1.append(v[j])
        j+=1

    t1=[]
    for j in range(1,len(t)):
        if(t[j]!=0):
            t1.append(t[j])
        j+=1

    err_Num_Formula = err_Num_Formula[0:len(v1)]

    # err_Num_Formulax10=(10**(7))*err_Num_Formula

    # print(len(err_Num_Formula))

    sse_root = np.linalg.norm(err_Num_Formula)
    n_root=(len(err_Num_Formula))**0.5

    norm_fn_of_h=sse_root/n_root
    return norm_fn_of_h



def A2Q3_CPART_ScnB(h):
              #Step-size
    i=0
    t=np.zeros(1000)
    v=np.zeros(1000)
    err=np.zeros(1000)
    t[1]=10
    v[1]=1000
    x=0
    err[x]=100
    # v[x]=0.000000000000001
    v[x]= 10**(-15)

    v_analytical = np.zeros(1000)
    err_Num_Formula= np.zeros(1000)
    g1 = 9.81
    m = 85
    rho = 1.21
    C = 0.7
    A = 0.18

    vss = ((2 * m * g1) / (C * rho * A)) ** 0.5
    tau = ((2 * m) / (C * rho * A * g1)) ** 0.5


    while( np.abs(err[x]/v[x]) > 10**(-6)   ):
        k1=g(t[i],v[i])
        k2=g(t[i] +h/2,v[i]+ h*(k1)/2  )
        k3=g(t[i] +h/2 ,v[i] + h*(k2)/2 )
        k4=g(t[i] +h , v[i] +h*k3)

        v[i+1]= v[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4 )
        t[i+1]=t[i] + h

        #v_analytical[i] = vss * (math.exp(2 * t[i] / tau) - 1) / (math.exp(2 * t[i] / tau) + 1)
        v_analytical[i]= vss*np.tanh(t[i]/tau)
        err_Num_Formula[i]= np.abs(v_analytical[i] -v[i])

        err[i]=v[i+1]-v[i]
        x=i
        i+=1
    # print(f'The velocity is {v[i]}')
    # print(f'The time required to reach is {t[i]}')
    # print(f'The no. of iterations is {i}')
    v1=[]

    for j in range(1,len(v)):
        if(v[j]!=0):
            v1.append(v[j])
        j+=1

    t1=[]
    for j in range(1,len(t)):
        if(t[j]!=0):
            t1.append(t[j])
        j+=1

    err_Num_Formula = err_Num_Formula[0:len(v1)]

    # err_Num_Formulax10=(10**(7))*err_Num_Formula

    # print(len(err_Num_Formula))

    sse_root = np.linalg.norm(err_Num_Formula)
    n_root=(len(err_Num_Formula))**0.5

    norm_fn_of_h=sse_root/n_root
    return norm_fn_of_h




def A2Q3_I():

    print('The answer to Scenario A is \n')
    A2Q3aScnA()
    print('\n')

    print('The answer to Scenario B is \n')
    A2Q3aScnB()
    print('\n')

    variation_A = np.zeros(4)
    variation_B = np.zeros(4)

    variation_A[0] = A2Q3_CPART_ScnA(0.1)
    variation_A[1] = A2Q3_CPART_ScnA(0.2)
    variation_A[2] = A2Q3_CPART_ScnA(0.5)
    variation_A[3] = A2Q3_CPART_ScnA(1)

    variation_B[0] = A2Q3_CPART_ScnB(0.1)
    variation_B[1] = A2Q3_CPART_ScnB(0.2)
    variation_B[2] = A2Q3_CPART_ScnB(0.5)
    variation_B[3] = A2Q3_CPART_ScnB(1)

    variation_A=variation_A*(10**(5))
    variation_B=variation_B*(10**(5))

    h=[0.1,0.2,0.5,1]

    plt.plot(h, variation_A, 'ro')
    plt.xlabel('Step sizes')
    plt.ylabel('RMSE valuesA x 10^(-5)')
    plt.title('RMS errors Scenario A ')
    plt.show()

    plt.plot(h, variation_B, 'bo')
    plt.xlabel('Step sizes')
    plt.ylabel('RMSE valuesB x 10^(-5)')
    plt.title('RMS errors Scenario B ')
    plt.show()

    # v_num,t_num =A2Q3aScnA()
    # print(v_num)
    # print(t_num)

A2Q3_I()

