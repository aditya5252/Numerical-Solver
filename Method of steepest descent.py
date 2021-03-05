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

def func3d(x, y):
    A = np.array([[3, 2], [2, 6]])
    b = np.array([[2],[-8]])
    x_arr=np.array([[x],[y]])
    x_arr_tr=np.transpose(x_arr)

    y= 0.5*np.dot(np.dot(x_arr_tr,A),x_arr) - np.dot(x_arr_tr,b)
    y1=y[0][0]
    return y1

def A2Q1():
    count= 0

    A = np.array([[3, 2], [2, 6]])
    b = np.array([[2],[-8]])

    x0 = np.array([[-2],[-2]])

    alpha=x0
    beta=np.array([[-1000],[-1000]])
    temp=x0
    soln=[]
    err=[]
    itr=[]
    residual=[]
    while( (np.linalg.norm( beta- alpha))**2 >= 0.001):

        alpha=temp

        nk= np.dot(A,alpha) - b
        nkt=np.transpose(nk)

        beta = alpha - (np.dot(nkt,nk))/( np.dot(np.dot(nkt, A),nk)   )*nk

        soln.append(beta)
        err.append((np.linalg.norm( beta- alpha)))
        residual.append(np.linalg.norm( b-  np.dot(A,beta) ))
        temp=beta

        count+=1

        itr.append(count)
    # print(len(residual))
    print(f" \n The total no. of iterations is  {count} ")
    # print(f" \n The final value of alpha is   {alpha} ")
    print(f" \n The final solution of the equation Ax=b is   {beta} ")


    plt.plot(itr,err)
    plt.xlabel('No. of the iterations')
    plt.ylabel('Error ||x(k+1)-x(k)|| ')
    plt.title('Error ||x(k+1)-x(k)|| with no. of iterations')
    plt.show()


    plt.plot(itr,residual)
    plt.xlabel('No. of the iterations')
    plt.ylabel('Residual ||b-A*x(k)|| ')
    plt.title('Residual ||b-A*x(k)|| with no. of iterations')
    plt.show()


    x_interval = (0, 4)
    y_interval = (-4, 0)
    # print(x_interval[1])

    x_points = np.linspace(x_interval[0], x_interval[1], 100)
    y_points = np.linspace(y_interval[0], y_interval[1], 100)

    X, Y = np.meshgrid(x_points, y_points)
    func3d_vectorized = np.vectorize(func3d)
    Z = func3d_vectorized(X, Y)

    plt.figure(figsize=(20, 10))
    ax = plt.axes(projection="3d")

    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='terrain', edgecolor=None)

    ax.set(xlabel='x', ylabel='y', zlabel='f(x, y)', title="Objective  Function")
    plt.show()

    fig, ax = plt.subplots()
    CS = ax.contour(X, Y, Z)
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_title('Contour plot Objective function')
    plt.xlabel('x axis')
    plt.ylabel('y axis')
    plt.show()

    x_interval = (-4, 4)
    y_interval = (-4, 4)
    # print(x_interval[1])

    x_points = np.linspace(x_interval[0], x_interval[1], 100)
    y_points = np.linspace(y_interval[0], y_interval[1], 100)

    X, Y = np.meshgrid(x_points, y_points)
    func3d_vectorized = np.vectorize(func3d)
    Z = func3d_vectorized(X, Y)

    fig, ax = plt.subplots()
    CS = ax.contour(X, Y, Z)
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_title('Contour plot with arrows Objective function')
    plt.xlabel('x axis')
    plt.ylabel('y axis')
    soln_new = []
    list = []
    list_new_x=[]
    list_new_y=[]
    for i in range(len(soln) - 1):
        soln_new.append(np.transpose(soln[i]))
        list.append(soln_new[i][0])
        list_new_x.append(list[i][0])
        list_new_y.append(list[i][1])
    c = np.transpose(x0)
    plt.plot([-2,list_new_x[0]], [-2,list_new_y[0]], '-k')
    for i in range(len(list)-1):
        plt.plot(list_new_x, list_new_y, '-k')
    plt.show()

#print(b)


A2Q1()