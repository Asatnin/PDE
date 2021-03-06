# coding=utf8
import algo
import matrix_utils
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np


# решение методом Гаусса
def first_method():
    step = 0.08
    x = np.arange(0, 1, step)
    z = np.arange(0, 1, step)
    X, Z = np.meshgrid(x, z)
    # X, Z = np.meshgrid(x[1:len(x) - 1], z[1:len(z) - 1])

    z_len = len(z)
    x_len = len(x)

    lambdas = [[0.0 for x in xrange(x_len + 2)] for z in xrange(z_len + 2)]
    ans = [20.0 for i in xrange(z_len * x_len)]
    ans_new = [20.0 for i in xrange(z_len * x_len)]
    lambdas = matrix_utils.recalc_lambda_coef(lambdas, ans_new)

    flag = True
    it = 0
    while flag or not algo.is_conv(ans, ans_new):
        flag = False  # :((
        it += 1

        # algo.is_conv(ans, ans_new)
        ans = ans_new

        lambdas = matrix_utils.recalc_lambda_coef(lambdas, ans)
        matr = matrix_utils.gen_matr(z_len, x_len, step, lambdas, ans)
        matr, ans_new = algo.gauss(matr)

    print "With " + str(it) + " iterations"

    # dirty hack
    # res = []
    # for i in xrange(1, z_len + 1):
    #     for j in xrange(1, x_len + 1):
    #         if i != 1 and i != z_len and j != 1 and j != x_len:
    #             res.append(ans_new[matrix_utils.number(i, j, z_len)])
    # ans_new[0] = ans_new[(x_len - 1) * z_len]
    # ans_new[x_len - 1] = ans_new[len(ans_new) - 1]
    # ans_new[x_len - 2] = ans_new[len(ans_new) - 2]
    # ans_new[len(ans_new) - 1] = -500

    # y = np.array(res)
    y = np.array(ans_new)
    Y = y.reshape(X.shape)
    draw_plot(X, Z, Y)


def draw_plot(X, Z, Y):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Z, Y, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Z Label')
    ax.set_zlabel('Y Label')
    plt.show()


# решение локально-одномерным методом
def second_method():
    t_step = 0.05  # шаг по времени
    step = 0.05
    x = np.arange(0, 1, step)
    z = np.arange(0, 1, step)
    X, Z = np.meshgrid(x, z)

    z_len = len(z)
    x_len = len(x)

    lambdas = [[0.0 for x in xrange(x_len + 2)] for z in xrange(z_len + 2)]
    ans = [20.0 for i in xrange(z_len * x_len)]
    ans_prom = [20.0 for i in xrange(z_len * x_len)]
    ans_new = [20.0 for i in xrange(z_len * x_len)]

    lambdas = matrix_utils.recalc_lambda_coef(lambdas, ans_new)

    flag = True
    it = 0
    while flag or not algo.is_conv(ans, ans_new):
        flag = False  # :((
        it += 1

        algo.is_conv(ans, ans_new)
        ans = list(ans_new)

        lambdas = matrix_utils.recalc_lambda_coef(lambdas, ans)

        for i in xrange(1, z_len + 1):
            pod, gl, nad, f = matrix_utils.gen_pr_matr_prom(i, z_len, x_len, step, t_step, lambdas, ans)
            res = algo.tdma(pod, nad, gl, f)
            for j in xrange(1, x_len + 1):
                ans_prom[matrix_utils.number(i, j, z_len)] = res[j - 1]

        lambdas = matrix_utils.recalc_lambda_coef(lambdas, ans_prom)

        for j in xrange(1, x_len + 1):
            pod, gl, nad, f = matrix_utils.gen_pr_matr_new(j, z_len, x_len, step, t_step, lambdas, ans_prom)
            res = algo.tdma(pod, nad, gl, f)
            for i in xrange(1, z_len + 1):
                ans_new[matrix_utils.number(i, j, z_len)] = res[i - 1]

    print "With " + str(it) + " iterations"

    y = np.array(ans_new)
    Y = y.reshape(X.shape)
    draw_plot(X, Z, Y)


first_method()
# second_method()

#
# y = np.array(ans)
# Y = y.reshape(X.shape)
# #print Y
#
# # plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")
# ax.plot_surface(X, Z, Y, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# ax.set_xlabel('X Label')
# ax.set_ylabel('Z Label')
# ax.set_zlabel('Y Label')
# plt.show()

# print algo.tdma([0, 1, 1, 1], [1, -5, 2], [2, 10, -5, 4], [5, -18, -40, -27])
# print algo.TDMA1([0, 1, 1, 1], [1, -5, 2], [2, 10, -5, 4], [5, -18, -40, -27])
# print
# print algo.tdma([0, -3, -5, -6, -5], [-1, -1, 2, -4], [2, 8, 12, 18, 10], [-25, 72, -69, -156, 20])
# print algo.TDMA1([0, -3, -5, -6, -5], [-1, -1, 2, -4], [2, 8, 12, 18, 10], [-25, 72, -69, -156, 20])
# print
# print algo.tdma([0, 3, 1, 1], [3, 1, -2], [5, 6, 4, -3], [8, 10, 3, -2])
# print algo.TDMA1([0, 3, 1, 1], [3, 1, -2], [5, 6, 4, -3], [8, 10, 3, -2])
# print
# print algo.tdma([0, 2, 1, 1], [1, -1, 3], [2, 3, -1, -1], [4, 9, 12, -4])
# print algo.TDMA1([0, 2, 1, 1], [1, -1, 3], [2, 3, -1, -1], [4, 9, 12, -4])


# a = [[3.0, 2.0, -5.0, -1.0],
#      [2.0, -1.0, 3.0, 13.0],
#      [1.0, 2.0, -1.0, 9.0]]
# a = [[4.0, 2.0, -1.0, 1.0],
#      [5.0, 3.0, -2.0, 2.0],
#      [3.0, 2.0, -3.0, 0.0]]
# a = [[8.0, 7.0, 3.0, 18.0],
#      [-7.0, -4.0, -4.0, -11.0],
#      [-6.0, 5.0, -4.0, -15.0]]
a = [[1.0, 2.0, 3.0, 1.0],
     [2.0, -1.0, 2.0, 6.0],
     [1.0, 1.0, 5.0, -1.0]]
# a = [[2.0, 5.0, 4.0, 1.0, 20.0],
#      [1.0, 3.0, 2.0, 1.0, 11.0],
#      [2.0, 10.0, 9.0, 7.0, 40.0],
#      [3.0, 8.0, 9.0, 2.0, 37.0]]
# a = [[4.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
#      [-1.0, 4.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
#      [0.0, -1.0, 4.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
#      [-1.0, 0.0, 0.0, 4.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0],
#      [0.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 0.0, 0.0],
#      [0.0, 0.0, -1.0, 0.0, -1.0, 4.0, 0.0, 0.0, -1.0, 0.0],
#      [0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 4.0, -1.0, 0.0, 0.0],
#      [0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0],
#      [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 4.0, 0.0]]
# a = [[4.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
#      [-1.0, 4.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2.0],
#      [0.0, -1.0, 4.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 3.0],
#      [-1.0, 0.0, 0.0, 4.0, -1.0, 0.0, -1.0, 0.0, 0.0, 4.0],
#      [0.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 0.0, 5.0],
#      [0.0, 0.0, -1.0, 0.0, -1.0, 4.0, 0.0, 0.0, -1.0, 6.0],
#      [0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 4.0, -1.0, 0.0, 7.0],
#      [0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 4.0, -1.0, 8.0],
#      [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 4.0, 9.0]]


#a, ans = algo.gauss(a)

# for i in xrange(len(a)):
#     for j in xrange(len(a[0])):
#         print a[i][j],
#     print

# print ans
