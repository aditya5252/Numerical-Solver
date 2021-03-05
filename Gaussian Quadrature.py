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



def A2Q3IIScnA():
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
    distance_analytical=y=np.empty([1, 1])
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
        terminal_condition=i
    print(f'The velocity is {v[i]}')
    print(f'The time required to reach is {t[i]}')
    print(f'The no. of iterations is {i}')
    v1=[0]

    for j in range(1,len(v)):
        if(v[j]!=0):
            v1.append(v[j])
        j+=1

    t1=[0]
    for j in range(1,len(t)):
        if(t[j]!=0):
            t1.append(t[j])
        j+=1

    err_Num_Formula = err_Num_Formula[0:len(v1)]
    # print(err_Num_Formula)
    err_Num_Formulax10=(10**(7))*err_Num_Formula

    # plt.title('Velocity vs Time Scenario A ')
    # plt.plot(t1, v1)
    # plt.xlabel('Time')
    # plt.ylabel('Velocity')
    # # plt.show()
    #
    # # plt.plot(t1, v1)
    # plt.plot(t1, err_Num_Formulax10)
    # plt.xlabel('Time')
    # plt.ylabel('b part err_Num_Formulax10^(-7))')
    # plt.title('Error analytical numerical Scenario A ')
    # # plt.show()



    # fx_interval_j=[]
    distance_per_interval=[0]
    for j in range(0,len(v1)-1):
        # print(v1[j])
        a=t1[j]
        b=t1[j+1]
        z=np.sqrt(3/5)
        x1=(b+a)/2 +0.5*(b-a)*(-z)
        x2=(b+a)/2
        x3=(b+a)/2 +0.5*(b-a)*(z)

        y1= v1[j]+ (v1[j+1]-v1[j])*(x1-a)/(b-a)
        y2= v1[j]+ (v1[j+1]-v1[j])*(x2-a)/(b-a)
        y3= v1[j]+ (v1[j+1]-v1[j])*(x3-a)/(b-a)

        d= 0.5*(b-a)*(  (5/9)*y1 + (8/9)*y2 + (5/9)*y3      )

        distance_per_interval.append(d)

    distance_analytical =  []
    for j in range(0, len(t1)):
        d=tau*vss*np.log( np.cosh(  t1[j]/tau  )  )
        # distance_analytical=np.append(distance_analytical,[d])
        distance_analytical.append(d)

    distance_analytical=np.array(distance_analytical)
    # print(distance_per_interval)
    x=0
    tot_distance=[]
    for i in range(0,len(distance_per_interval)):
        x+=distance_per_interval[i]
        tot_distance.append(x)
    # print(tot_distance)
    tot_distance_m = np.array(tot_distance)


    #Total distance from the top point
    tot_distance_ft=tot_distance_m*3.280839895
    # print(tot_distance_ft)
    distance_analytical_ft=distance_analytical*3.280839895
    # print(distance_analytical_ft)

    #Error between analytical and Numerical distances

    err_dist=np.abs(distance_analytical_ft - tot_distance_ft)

    # t_terminal=t1[len(t1)-1]
    t_terminal=t1[terminal_condition]
    v_terminal=v1[terminal_condition]

    #a) part answer
    #Time at which parachute open
    t_parachute=t_terminal + (13500 -5000 - tot_distance_ft[len(tot_distance_ft)-1])/(v_terminal*3.280839895)

    print(f'\nThe time (in sec) at which Parachute is opened in Scenario A is {t_parachute} \n')

    print(f'The total distance covered from the top (in ft) before reaching terminal velocity is {tot_distance_ft[len(tot_distance_ft) - 1]}\n')

    print('The time during which the diver experiences terminal velocity')
    print(f' till the time he opens the parachute is {t_parachute-t_terminal} sec\n')
    print('The maximum velocity reached is the terminal velocity after which')
    print(f' it stays the same which is given by {v_terminal*18/5} kmph')

    #
    # print(len(err_dist))
    # print(len(t1))
    #
    # print(err_dist)
    # print(t1)

    # err_dist=err_dist.tolist()
    # print(err_dist)

    plt.plot(t1, err_dist)
    plt.xlabel('Time')
    plt.ylabel('|d_analytical - d_numerical | in feet')
    plt.title('Distance error plots Scenario A ')
    plt.show()

A2Q3IIScnA()





def A2Q3IIScnB():
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
    C = 0.7
    A = 0.18

    vss = ((2 * m * g1) / (C * rho * A)) ** 0.5
    tau = ((2 * m) / (C * rho * A * g1)) ** 0.5

    # print(tau)
    distance_analytical=y=np.empty([1, 1])
    while( np.abs(err[x]/v[x]) > 10**(-6)   ):
        k1=g(t[i],v[i])
        k2=g(t[i] +h/2,v[i]+ h*(k1)/2  )
        k3=g(t[i] +h/2 ,v[i] + h*(k2)/2 )
        k4=g(t[i] +h , v[i] +h*k3)

        v[i+1]=v[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4 )
        t[i+1]=t[i] + h

        v_analytical[i] = vss * (math.exp(2 * t[i] / tau) - 1) / (math.exp(2 * t[i] / tau) + 1)
        # v_analytical[i]=vss*np.tanh(t[i]/tau)
        err_Num_Formula[i]= np.abs(v_analytical[i]-v[i])

        err[i]=v[i+1]-v[i]
        x=i
        i+=1
        terminal_condition=i
    print(f'The velocity is {v[i]}')
    print(f'The time required to reach is {t[i]}')
    print(f'The no. of iterations is {i}')
    v1=[0]

    for j in range(1,len(v)):
        if(v[j]!=0):
            v1.append(v[j])
        j+=1

    t1=[0]
    for j in range(1,len(t)):
        if(t[j]!=0):
            t1.append(t[j])
        j+=1

    err_Num_Formula = err_Num_Formula[0:len(v1)]
    # print(err_Num_Formula)
    err_Num_Formulax10=(10**(7))*err_Num_Formula

    # plt.title('Velocity vs Time Scenario A ')
    # plt.plot(t1, v1)
    # plt.xlabel('Time')
    # plt.ylabel('Velocity')
    # # plt.show()
    #
    # # plt.plot(t1, v1)
    # plt.plot(t1, err_Num_Formulax10)
    # plt.xlabel('Time')
    # plt.ylabel('b part err_Num_Formulax10^(-7))')
    # plt.title('Error analytical numerical Scenario A ')
    # # plt.show()



    # fx_interval_j=[]
    distance_per_interval=[0]
    for j in range(0,len(v1)-1):
        # print(v1[j])
        a=t1[j]
        b=t1[j+1]
        z=np.sqrt(3/5)
        x1=(b+a)/2 +0.5*(b-a)*(-z)
        x2=(b+a)/2
        x3=(b+a)/2 +0.5*(b-a)*(z)

        y1= v1[j]+ (v1[j+1]-v1[j])*(x1-a)/(b-a)
        y2= v1[j]+ (v1[j+1]-v1[j])*(x2-a)/(b-a)
        y3= v1[j]+ (v1[j+1]-v1[j])*(x3-a)/(b-a)

        d= 0.5*(b-a)*(  (5/9)*y1 + (8/9)*y2 + (5/9)*y3      )

        distance_per_interval.append(d)

    distance_analytical =  []
    for j in range(0, len(t1)):
        d=tau*vss*np.log( np.cosh(  t1[j]/tau  )  )
        # distance_analytical=np.append(distance_analytical,[d])
        distance_analytical.append(d)

    distance_analytical=np.array(distance_analytical)
    # print(distance_per_interval)
    x=0
    tot_distance=[]
    for i in range(0,len(distance_per_interval)):
        x+=distance_per_interval[i]
        tot_distance.append(x)
    # print(tot_distance)
    tot_distance_m = np.array(tot_distance)


    #Total distance from the top point
    tot_distance_ft=tot_distance_m*3.280839895
    # print(tot_distance_ft)
    distance_analytical_ft=distance_analytical*3.280839895
    # print(distance_analytical_ft)

    #Error between analytical and Numerical distances

    err_dist=np.abs(distance_analytical_ft - tot_distance_ft)

    # t_terminal=t1[len(t1)-1]
    t_terminal=t1[terminal_condition]
    v_terminal=v1[terminal_condition]

    #a) part answer
    #Time at which parachute open
    t_parachute=t_terminal + (13500 -5000 - tot_distance_ft[len(tot_distance_ft)-1])/(v_terminal*3.280839895)

    print(f'\nThe time (in sec) at which Parachute is opened in Scenario B is {t_parachute} \n')

    print(f'The total distance covered from the top (in ft) before reaching terminal velocity is {tot_distance_ft[len(tot_distance_ft) - 1]}\n')

    print('The time during which the diver experiences terminal velocity')
    print(f' till the time he opens the parachute is {t_parachute-t_terminal} sec\n')
    print('The maximum velocity reached is the terminal velocity after which')
    print(f' it stays the same which is given by {v_terminal*18/5} kmph')

    #
    # print(len(err_dist))
    # print(len(t1))
    #
    # print(err_dist)
    # print(t1)

    # err_dist=err_dist.tolist()
    # print(err_dist)

    plt.plot(t1, err_dist)
    plt.xlabel('Time')
    plt.ylabel('|d_analytical - d_numerical| in feet')
    plt.title('Distance error plots Scenario B ')
    plt.show()

    print('\nWe are getting these absurd results because the diver would never reach terminal velocity  ')
    print('For safe USPA guidelines he will have to open at 5000 ft (abv ground altitude) ')
    print('Or 8500 ft from the starting point ')
    for i in range(0,len(tot_distance_ft)):
        if (tot_distance_ft[i]> 8490 and tot_distance_ft[i]< 8550):
            parachute_release_time=t1[i]
            break
    print(f'\nThe parachute release time in Scenario B is {parachute_release_time} sec')


    print('Since the terminal velocity is never reached:')
    print('The distance covered in FREE FALL before reaching terminal velocity is 8500 ft ')
    print('The time during which the diver experiences terminal velocity is 0 sec since he never experiences it ')
    # print('The maximum velocity reached is the terminal velocity after which')
    # print(f' it stays the same which is given by {v_terminal*18/5} kmph')
    # print(tot_distance_ft)

    print('Maximum velocity reached is Velocity before deploying parachute since parachute slows him down')
    print(f'Maximum Velocity is {v1[i]*18/5}')


A2Q3IIScnB()

