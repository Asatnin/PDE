# coding=utf8
import algo
import matrix_utils
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np

# решение методом Гаусса
def first_method():
     step = 0.05
     x = np.arange(0, 1, step)
     z = np.arange(0, 1, step)
     X, Z = np.meshgrid(x, z)

     z_len = len(z)
     x_len = len(x)

     lambdas = [[0.0 for z in xrange(z_len + 2)] for x in xrange(x_len + 2)]
     ans = ans_new = [20.0 for i in xrange(z_len * x_len)]
     lambdas = matrix_utils.recalc_lambda_coef(lambdas, ans_new)

     flag = 1
     while flag < 2: #or not algo.is_conv(ans, ans_new):
          flag += 1  # :((

          algo.is_conv(ans, ans_new)
          ans = ans_new

          lambdas = matrix_utils.recalc_lambda_coef(lambdas, ans)
          matr = matrix_utils.gen_matr(z_len, x_len, step, lambdas, ans)
          matr, ans_new = algo.gauss(matr)

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


first_method()

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
