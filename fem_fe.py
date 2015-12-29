#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#           Реализация функций форм конечных элементов
###################################################################

import math
from abc import abstractmethod
from numpy.linalg import solve
from numpy import array
from numpy import zeros
from numpy.linalg import det
from numpy.linalg import inv
from fem_error import TFEMException
from fem_defs import eps


# Абстрактный базовый класс, описывающий конечный элемент (КЭ)
class TFE:
    def __init__(self):
        self.size = 0
        self.density = self.damping = 0.0
        self.e = []  # Модуль Юнга и коэффициент Пуассона
        self.m = []
        self.x = []  # Координаты вершин КЭ
        self.y = []
        self.z = []
        self.vx = []  # Компоненты вектора объемной нагрузки
        self.vy = []
        self.vz = []
        self.K = [[]]  # Матрица жесткости
        self.M = [[]]  # ... масс
        self.D = [[]]  # ... демпфирования
        self.c = [[]]  # Коэффициенты функций форм

    # Задание координат
    def set_coord(self, *args):
        if len(args) == 1:
            self.x = args[0]
        elif len(args) == 2:
            self.x = args[0]
            self.y = args[1]
        elif len(args) == 3:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
        self.__create__()

    # Задание объемной нагрузки
    def set_volume_load(self, *args):
        if len(args) == 1:
            self.vx = args[0]
        elif len(args) == 2:
            self.vx = args[0]
            self.vy = args[1]
        elif len(args) == 3:
            self.vx = args[0]
            self.vy = args[1]
            self.vz = args[2]

    # Задание параметров упругости
    def set_elasticity(self, p1, p2):
        self.e = p1
        self.m = p2

    # Задание плотности
    def set_density(self, d):
        self.density = d

    # Задание параметра демпфирования
    def set_damping(self, d):
        self.damping = d

    @abstractmethod
    # Вычисление функций форм КЭ
    def __create__(self):
        pass

    # Формирование матриц жесткости, масс и демпфирования
    def generate(self, is_static=True):
        pass

    # Вычисление объема (площади) КЭ
    def volume(self):
        pass

    # Вычисления стандартных результатов КЭ
    def calc(self, u, index):
        pass


# Линейный (двухузловой) одномерный КЭ
class TFE1D2(TFE):
    def __init__(self):
        super().__init__()
        self.size = 2
        self.K = [
            [0, 0, 0],
            [0, 0, 0]
        ]
        self.c = [
            [0, 0],
            [0, 0]
        ]

    def volume(self):
        return math.fabs(self.x[1] - self.x[0])

    def __create__(self):
        if self.volume() == 0.0:
            raise TFEMException('incorrect_fe_err')
        vol = self.x[1] - self.x[0]
        self.c[0][0] = self.x[1]/vol
        self.c[0][1] = 1.0/vol
        self.c[1][0] = self.x[0]/vol
        self.c[1][1] = 1.0/vol

    def calc(self, u, index):
        u0 = u[0]
        u1 = u[1]
        if index == 0:
            res = u0*self.c[0][1] + u1*self.c[1][1]
        else:
            res = self.e[0]*(u0*self.c[0][1] + u1*self.c[1][1])
        u[0] = u[1] = res

    def generate(self, is_static=True):
        self.K[0][0] = self.volume()*2.0*self.e[0]*self.c[0][1]*self.c[0][1]
        self.K[0][1] = self.volume()*2.0*self.e[0]*self.c[0][1]*self.c[1][1]
        self.K[1][1] = 2.0*self.e[0]*self.c[1][1]*self.c[1][1]

        # Вычисление интеграла для объемных сил
        if len(self.vx):
            self.K[0][2] += self.vx[0]*(self.c[0][1]*(self.x[1]*self.x[1] - self.x[0]*self.x[0])*0.5 +
                                        self.c[0][0]*(self.x[1] - self.x[0]))
            self.K[1][2] += self.vx[1]*(self.c[1][1]*(self.x[1]*self.x[1] - self.x[0]*self.x[0])*0.5 +
                                        self.c[1][0]*(self.x[1] - self.x[0]))
        if not is_static:
            # Формирование матрицы массы
            k00 = self.volume()*(-2.0/3.0*self.c[0][1]*self.c[0][1]*self.x[0]*self.x[0]*self.x[0] +
                                 2.0/3.0*self.c[0][1]*self.c[0][1]*self.x[1]*self.x[1]*self.x[1] -
                                 2.0*self.c[0][0]*self.c[0][0]*self.x[0] + 2.0*self.c[0][0]*self.c[0][0]*self.x[1] -
                                 2.0*self.c[0][0]*self.c[0][1]*self.x[0]*self.x[0] +
                                 2.0*self.c[0][0]*self.c[0][1]*self.x[1]*self.x[1])
            k01 = self.volume()*(-2.0/3.0*self.c[0][1]*self.c[1][1]*self.x[0]*self.x[0]*self.x[0] +
                                 2.0/3.0*self.c[0][1]*self.c[1][1]*self.x[1]*self.x[1]*self.x[1] -
                                 self.c[0][1]*self.c[1][0]*self.x[0]*self.x[0] +
                                 self.c[0][1]*self.c[1][0]*self.x[1]*self.x[1] -
                                 self.c[0][0]*self.c[1][1]*self.x[0]*self.x[0] +
                                 self.c[0][0]*self.c[1][1]*self.x[1]*self.x[1] -
                                 2.0*self.c[0][0]*self.c[1][0]*self.x[0] + 2.0*self.c[0][0]*self.c[1][0]*self.x[1])
            k11 = self.volume()*(2.0*self.c[1][0]*self.c[1][1]*self.x[1]*self.x[1] -
                                 2.0*self.c[1][0]*self.c[1][1]*self.x[0]*self.x[0] +
                                 2.0/3.0*self.c[1][1]*self.c[1][1]*self.x[1]*self.x[1]*self.x[1] -
                                 2.0/3.0*self.c[1][1]*self.c[1][1]*self.x[0]*self.x[0]*self.x[0] +
                                 2.0*self.c[1][0]*self.c[1][0]*self.x[1] - 2.0*self.c[1][0]*self.c[1][0]*self.x[0])

            self.M = [[0, 0], [0, 0]]
            self.M[0][0] = self.density*k00
            self.M[0][1] = self.density*k01
            self.M[1][0] = self.density*k01
            self.M[1][1] = self.density*k11

            # Формирование матрицы демпфирования
            self.D = [[0, 0], [0, 0]]
            self.D[0][0] = self.damping*k00
            self.D[0][1] = self.damping*k01
            self.D[1][0] = self.damping*k01
            self.D[1][1] = self.damping*k11


# Линейный (трехузловой) треугольный КЭ
class TFE2D3(TFE):
    def __init__(self):
        super().__init__()
        self.size = 3
        self.K = [
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0]
        ]
        self.c = [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ]

    def volume(self):
        a = math.sqrt((self.x[0] - self.x[1])*(self.x[0] - self.x[1]) + (self.y[0] - self.y[1])*(self.y[0] - self.y[1]))
        b = math.sqrt((self.x[0] - self.x[2])*(self.x[0] - self.x[2]) + (self.y[0] - self.y[2])*(self.y[0] - self.y[2]))
        c = math.sqrt((self.x[2] - self.x[1])*(self.x[2] - self.x[1]) + (self.y[2] - self.y[1])*(self.y[2] - self.y[1]))
        p = 0.5*(a + b + c)
        return math.sqrt(p*(p - a)*(p - b)*(p - c))

    def __create__(self):
        det = self.y[2]*self.x[1] - self.y[2]*self.x[0] - self.y[0]*self.x[1] - self.y[1]*self.x[2] + \
              self.y[1]*self.x[0] + self.y[0]*self.x[2]
        det0 = self.y[2]*self.x[1] - self.y[1]*self.x[2]
        det1 = self.y[1] - self.y[2]
        det2 = self.x[2] - self.x[1]

        if math.fabs(det) < eps:
            raise TFEMException('incorrect_fe_err')
        self.c[0][0] = det0/det
        self.c[0][1] = det1/det
        self.c[0][2] = det2/det

        det0 = -self.y[2]*self.x[0] + self.y[0]*self.x[2]
        det1 = self.y[2] - self.y[0]
        det2 = -self.x[2] + self.x[0]

        self.c[1][0] = det0/det
        self.c[1][1] = det1/det
        self.c[1][2] = det2/det

        det0 = -self.y[0]*self.x[1] + self.y[1]*self.x[0]
        det1 = -self.y[1] + self.y[0]
        det2 = self.x[1] - self.x[0]

        self.c[2][0] = det0/det
        self.c[2][1] = det1/det
        self.c[2][2] = det2/det

    def calc(self, u, index):
        u0 = u[0]
        v0 = u[1]
        u1 = u[2]
        v1 = u[3]
        u2 = u[4]
        v2 = u[5]
        res = 0
        c01 = self.c[0][1]
        c02 = self.c[0][2]
        c11 = self.c[1][1]
        c12 = self.c[1][2]
        c21 = self.c[2][1]
        c22 = self.c[2][2]
        e = self.e[0]
        m = self.m[0]
        if index == 0:      # Exx
            res = u0*c01 + u1*c11 + u2*c21
        elif index == 1:    # Eyy
            res = v0*c02 + v1*c12 + v2*c22
        elif index == 2:    # Exy
            res = u0*c02 + u1*c12 + u2*c22 + v0*c01 + v1*c11 + v2*c21
        elif index == 3:    # Sxx
            res = e/(1.0 - m*m)*((u0*c01 + u1*c11 + u2*c21) + self.m[0]*(v0*c02 + v1*c12 + v2*c22))
        elif index == 4:    # Syy
            res = e/(1.0 - m*m)*(self.m[0]*(u0*c01 + u1*c11 + u2*c21) + (v0*c02 + v1*c12 + v2*c22))
        elif index == 5:    # Sxy
            res = e/(2.0 + 2.0*m)*(u0*c02 + u1*c12 + u2*c22 + v0*c01 + v1*c11 + v2*c21)
        u[0] = u[1] = u[2] = res

    def generate(self, is_static=True):
        k = self.e[0]/(1.0 - self.m[0]*self.m[0])
        g = self.e[0]/(2.0 + 2.0*self.m[0])
        vol = self.volume()

        self.K[0][0] = vol*(k*self.c[0][1]*self.c[0][1] + g*self.c[0][2]*self.c[0][2])
        self.K[0][1] = vol*(k*self.c[0][1]*self.m[0]*self.c[0][2] + g*self.c[0][2]*self.c[0][1])
        self.K[0][2] = vol*(k*self.c[0][1]*self.c[1][1] + g*self.c[0][2]*self.c[1][2])
        self.K[0][3] = vol*(k*self.c[0][1]*self.m[0]*self.c[1][2] + g*self.c[0][2]*self.c[1][1])
        self.K[0][4] = vol*(k*self.c[0][1]*self.c[2][1] + g*self.c[0][2]*self.c[2][2])
        self.K[0][5] = vol*(g*self.c[0][2]*self.c[2][1] + k*self.c[0][1]*self.m[0]*self.c[2][2])

        self.K[1][1] = vol*(k*self.c[0][2]*self.c[0][2] + g*self.c[0][1]*self.c[0][1])
        self.K[1][2] = vol*(g*self.c[1][2]*self.c[0][1] + k*self.c[1][1]*self.m[0]*self.c[0][2])
        self.K[1][3] = vol*(g*self.c[0][1]*self.c[1][1] + k*self.c[0][2]*self.c[1][2])
        self.K[1][4] = vol*(k*self.c[2][1]*self.m[0]*self.c[0][2] + g*self.c[2][2]*self.c[0][1])
        self.K[1][5] = vol*(g*self.c[0][1]*self.c[2][1] + k*self.c[0][2]*self.c[2][2])

        self.K[2][2] = vol*(k*self.c[1][1]*self.c[1][1] + g*self.c[1][2]*self.c[1][2])
        self.K[2][3] = vol*(g*self.c[1][2]*self.c[1][1] + k*self.c[1][1]*self.m[0]*self.c[1][2])
        self.K[2][4] = vol*(g*self.c[1][2]*self.c[2][2] + k*self.c[1][1]*self.c[2][1])
        self.K[2][5] = vol*(g*self.c[1][2]*self.c[2][1] + k*self.c[1][1]*self.m[0]*self.c[2][2])

        self.K[3][3] = vol*(g*self.c[1][1]*self.c[1][1] + k*self.c[1][2]*self.c[1][2])
        self.K[3][4] = vol*(g*self.c[2][2]*self.c[1][1] + k*self.c[2][1]*self.m[0]*self.c[1][2])
        self.K[3][5] = vol*(k*self.c[1][2]*self.c[2][2] + g*self.c[1][1]*self.c[2][1])

        self.K[4][4] = vol*(g*self.c[2][2]*self.c[2][2] + k*self.c[2][1]*self.c[2][1])
        self.K[4][5] = vol*(g*self.c[2][2]*self.c[2][1] + k*self.c[2][1]*self.m[0]*self.c[2][2])

        self.K[5][5] = vol*(g*self.c[2][1]*self.c[2][1] + k*self.c[2][2]*self.c[2][2])

        # Вычисление интеграла для объемных сил
        if len(self.vx) or len(self.vy):
            self.K[0][6] += self.vx[0]*vol/6.0
            self.K[2][6] += self.vx[1]*vol/6.0
            self.K[4][6] += self.vx[2]*vol/6.0
            self.K[1][6] += self.vy[0]*vol/6.0
            self.K[3][6] += self.vy[1]*vol/6.0
            self.K[5][6] += self.vy[2]*vol/6.0

        if not is_static:
            # Формирование матриц массы и демпфирования
            k00 = vol/12.0
            k01 = vol/6.0
            self.M = [
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0]
            ]
            self.M[0][0] = k00
            self.M[0][2] = k01
            self.M[0][4] = k01
            self.M[1][1] = k00
            self.M[1][3] = k01
            self.M[1][5] = k01
            self.M[2][0] = k01
            self.M[2][2] = k00
            self.M[2][4] = k01
            self.M[3][1] = k01
            self.M[3][3] = k00
            self.M[3][5] = k01
            self.M[4][0] = k01
            self.M[4][2] = k01
            self.M[4][4] = k00
            self.M[5][1] = k01
            self.M[5][3] = k01
            self.M[5][5] = k00

            self.D = list(self.M)
            for i in range(0, len(self.M)):
                for j in range(0, len(self.M)):
                    self.M[i][j] *= self.density
                    self.D[i][j] *= self.damping


# Линейный (четырехузловой) тетраэдральный КЭ
class TFE3D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4
        self.K = [[0]*13]*12
        self.c = [[0]*4]*4

    def volume(self):
        a = (self.x[1] - self.x[0])*(self.y[2] - self.y[0])*(self.z[3] - self.z[0]) + \
            (self.x[3] - self.x[0])*(self.y[1] - self.y[0])*(self.z[2] - self.z[0]) + \
            (self.x[2] - self.x[0])*(self.y[3] - self.y[0])*(self.z[1] - self.z[0])
        b = (self.x[3] - self.x[0])*(self.y[2] - self.y[0])*(self.z[1] - self.z[0]) + \
            (self.x[2] - self.x[0])*(self.y[1] - self.y[0])*(self.z[3] - self.z[0]) + \
            (self.x[1] - self.x[0])*(self.y[3] - self.y[0])*(self.z[2] - self.z[0])
        return math.fabs(a - b)/6.0

    def __create__(self):
        if self.volume() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a = array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i]
                a[i][2] = self.y[i]
                a[i][3] = self.z[i]
            b[j] = 1.0
            x = solve(a, b)
            self.c[j] = list(x)

    def calc(self, u, index):
        e = self.e[0]
        m = self.m[0]
        g = e/(2.0 + 2.0*m)
        l = 2.0*m*g/(1.0 - 2.0*m)
        u0 = u[0]
        v0 = u[1]
        w0 = u[2]
        u1 = u[3]
        v1 = u[4]
        w1 = u[5]
        u2 = u[6]
        v2 = u[7]
        w2 = u[8]
        u3 = u[9]
        v3 = u[10]
        w3 = u[11]
        c01 = self.c[0][1]
        c02 = self.c[0][2]
        c03 = self.c[0][3]
        c11 = self.c[1][1]
        c12 = self.c[1][2]
        c13 = self.c[1][3]
        c21 = self.c[2][1]
        c22 = self.c[2][2]
        c23 = self.c[2][3]
        c31 = self.c[3][1]
        c32 = self.c[3][2]
        c33 = self.c[3][3]
        res = 0
        if index == 0:      # Exx
            res = u0*c01 + u1*c11 + u2*c21 + u3*c31
        elif index == 1:    # Eyy
            res = v0*c02 + v1*c12 + v2*c22 + v3*c32
        elif index == 2:    # Ezz
            res = w0*c03 + w1*c13 + w2*c23 + w3*c33
        elif index == 3:    # Exy
            res = u0*c02 + u1*c12 + u2*c22 + u3*c32 + v0*c01 + v1*c11 + v2*c21 + v3*c31
        elif index == 4:    # Exz
            res = u0*c03 + u1*c13 + u2*c23 + u3*c33 + w0*c01 + w1*c11 + w2*c21 + w3*c31
        elif index == 5:    # Eyz
            res = v0*c03 + v1*c13 + v2*c23 + v3*c33 + w0*c02 + w1*c12 + w2*c22 + w3*c32
        elif index == 6:    # Sxx
            res = 2.0*g*u0*c01 + 2.0*g*u1*c11 + 2.0*g*u2*c21 + 2.0*g*u3*c31 + l*u0*c01 + l*u1*c11 + l*u2*c21 + \
                  l*u3*c31 + l*v0*c02 + l*v1*c12 + l*v2*c22 + l*v3*c32 + l*w0*c03 + l*w1*c13 + l*w2*c23 + l*w3*c33
        elif index == 7:    # Syy
            res = 2.0*g*v0*c02 + 2.0*g*v1*c12 + 2.0*g*v2*c22 + 2.0*g*v3*c32 + l*u0*c01 + l*u1*c11 + l*u2*c21 + \
                  l*u3*c31 + l*v0*c02 + l*v1*c12 + l*v2*c22 + l*v3*c32 + l*w0*c03 + l*w1*c13 + l*w2*c23 + l*w3*c33
        elif index == 8:    # Szz
            res = 2.0*g*w0*c03 + 2.0*g*w1*c13 + 2.0*g*w2*c23 + 2.0*g*w3*c33 + l*u0*c01 + l*u1*c11 + l*u2*c21 + \
                  l*u3*c31 + l*v0*c02 + l*v1*c12 + l*v2*c22 + l*v3*c32 + l*w0*c03 + l*w1*c13 + l*w2*c23 + l*w3*c33
        elif index == 9:    # Sxy
            res = g*(u0*c02 + u1*c12 + u2*c22 + u3*c32 + v0*c01 + v1*c11 + v2*c21 + v3*c31)
        elif index == 10:   # Sxz
            res = g*(u0*c03 + u1*c13 + u2*c23 + u3*c33 + w0*c01 + w1*c11 + w2*c21 + w3*c31)
        elif index == 11:   # Syz
            res = g*(v0*c03 + v1*c13 + v2*c23 + v3*c33 + w0*c02 + w1*c12 + w2*c22 + w3*c32)
        u[0] = u[1] = u[2] = u[3] = u[4] = res

    def generate(self, is_static=True):
        vol = self.volume()
        g = self.e[0]/(2.0 + 2.0*self.m[0])
        l = 2.0*self.m[0]*g/(1.0 - 2.0*self.m[0])
        c01 = self.c[0][1]
        c02 = self.c[0][2]
        c03 = self.c[0][3]
        c11 = self.c[1][1]
        c12 = self.c[1][2]
        c13 = self.c[1][3]
        c21 = self.c[2][1]
        c22 = self.c[2][2]
        c23 = self.c[2][3]
        c31 = self.c[3][1]
        c32 = self.c[3][2]
        c33 = self.c[3][3]

        self.K[0][0] = vol*0.5*(2.0*(2.0*g*c01 + l*c01)*c01 + 2.0*g*c02*c02 + 2.0*g*c03*c03)
        self.K[0][1] = vol*0.5*(2.0*l*c02*c01 + 2.0*g*c01*c02)
        self.K[0][2] = vol*0.5*(2.0*l*c01*c03 + 2.0*g*c01*c03)
        self.K[0][3] = vol*0.5*((2.0*g*c01 + l*c01)*c11 + (2.0*g*c11 + l*c11)*c01 + 2.0*g*c12*c02 + 2.0*g*c13*c03)
        self.K[0][4] = vol*0.5*(2.0*g*c11*c02 + 2.0*l*c12*c01)
        self.K[0][5] = vol*0.5*(2.0*l*c13*c01 + 2.0*g*c11*c03)
        self.K[0][6] = vol*0.5*((2.0*g*c01 + l*c01)*c21 + 2.0*g*c23*c03 + (2.0*g*c21 + l*c21)*c01 + 2.0*g*c22*c02)
        self.K[0][7] = vol*0.5*(2.0*l*c22*c01 + 2.0*g*c21*c02)
        self.K[0][8] = vol*0.5*(2.0*l*c23*c01 + 2.0*g*c21*c03)
        self.K[0][9] = vol*0.5*(2.0*g*c32*c02 + (2.0*g*c31 + l*c31)*c01 + 2.0*g*c33*c03 + (2.0*g*c01 + l*c01)*c31)
        self.K[0][10] = vol*0.5*(2.0*l*c32*c01 + 2.0*g*c31*c02)
        self.K[0][11] = vol*0.5*(2.0*l*c33*c01 + 2.0*g*c31*c03)

        self.K[1][1] = vol*0.5*(2.0*(2.0*g*c02 + l*c02)*c02 + 2.0*g*c03*c03 + 2.0*g*c01*c01)
        self.K[1][2] = vol*0.5*(2.0*l*c03*c02 + 2.0*g*c03*c02)
        self.K[1][3] = vol*0.5*(2.0*l*c02*c11 + 2.0*g*c12*c01)
        self.K[1][4] = vol*0.5*(2.0*g*c13*c03 + (2.0*g*c12 + l*c12)*c02 + (2.0*g*c02 + l*c02)*c12 + 2.0*g*c01*c11)
        self.K[1][5] = vol*0.5*(2.0*g*c03*c12 + 2.0*l*c13*c02)
        self.K[1][6] = vol*0.5*(2.0*l*c02*c21 + 2.0*g*c22*c01)
        self.K[1][7] = vol*0.5*((2.0*g*c02 + l*c02)*c22 + (2.0*g*c22 + l*c22)*c02 + 2.0*g*c01*c21 + 2.0*g*c23*c03)
        self.K[1][8] = vol*0.5*(2.0*l*c23*c02 + 2.0*g*c03*c22)
        self.K[1][9] = vol*0.5*(2.0*l*c02*c31 + 2.0*g*c32*c01)
        self.K[1][10] = vol*0.5*((2.0*g*c02 + l*c02)*c32 + (2.0*g*c32 + l*c32)*c02 + 2.0*g*c33*c03 + 2.0*g*c01*c31)
        self.K[1][11] = vol*0.5*(2.0*g*c03*c32 + 2.0*l*c33*c02)
    
        self.K[2][2] = vol*0.5*(2.0*g*c01*c01 + 2*(2.0*g*c03 + l*c03)*c03 + 2.0*g*c02*c02)
        self.K[2][3] = vol*0.5*(2.0*g*c13*c01 + 2.0*l*c03*c11)
        self.K[2][4] = vol*0.5*(2.0*l*c03*c12 + 2.0*g*c13*c02)
        self.K[2][5] = vol*0.5*((2.0*g*c13 + l*c13)*c03 + (2.0*g*c03 + l*c03)*c13 + 2.0*g*c01*c11 + 2.0*g*c12*c02)
        self.K[2][6] = vol*0.5*(2.0*l*c03*c21 + 2.0*g*c23*c01)
        self.K[2][7] = vol*0.5*(2.0*l*c03*c22 + 2.0*g*c23*c02)
        self.K[2][8] = vol*0.5*((2.0*g*c03 + l*c03)*c23 + 2.0*g*c01*c21 + 2.0*g*c22*c02 + (2.0*g*c23 + l*c23)*c03)
        self.K[2][9] = vol*0.5*(2.0*g*c33*c01 + 2.0*l*c03*c31)
        self.K[2][10] = vol*0.5*(2.0*g*c33*c02 + 2.0*l*c03*c32)
        self.K[2][11] = vol*0.5*((2.0*g*c03 + l*c03)*c33 + 2.0*g*c01*c31 + 2.0*g*c32*c02 + (2.0*g*c33 + l*c33)*c03)
    
        self.K[3][3] = vol*0.5*(2.0*(2.0*g*c11 + l*c11)*c11 + 2.0*g*c12*c12 + 2.0*g*c13*c13)
        self.K[3][4] = vol*0.5*(2.0*l*c12*c11 + 2.0*g*c12*c11)
        self.K[3][5] = vol*0.5*(2.0*g*c13*c11 + 2.0*l*c13*c11)
        self.K[3][6] = vol*0.5*((2.0*g*c11 + l*c11)*c21 + 2.0*g*c13*c23 + (2.0*g*c21 + l*c21)*c11 + 2.0*g*c12*c22)
        self.K[3][7] = vol*0.5*(2.0*g*c12*c21 + 2.0*l*c22*c11)
        self.K[3][8] = vol*0.5*(2.0*l*c23*c11 + 2.0*g*c13*c21)
        self.K[3][9] = vol*0.5*(2.0*g*c12*c32 + (2.0*g*c31 + l*c31)*c11 + 2.0*g*c13*c33 + (2.0*g*c11 + l*c11)*c31)
        self.K[3][10] = vol*0.5*(2.0*l*c32*c11 + 2.0*g*c12*c31)
        self.K[3][11] = vol*0.5*(2.0*l*c33*c11 + 2.0*g*c13*c31)
    
        self.K[4][4] = vol*0.5*(2.0*g*c13*c13 + 2*(2.0*g*c12 + l*c12)*c12 + 2.0*g*c11*c11)
        self.K[4][5] = vol*0.5*(2.0*g*c13*c12 + 2.0*l*c13*c12)
        self.K[4][6] = vol*0.5*(2.0*l*c12*c21 + 2.0*g*c22*c11)
        self.K[4][7] = vol*0.5*((2.0*g*c12 + l*c12)*c22 + (2.0*g*c22 + l*c22)*c12 + 2.0*g*c11*c21 + 2.0*g*c13*c23)
        self.K[4][8] = vol*0.5*(2.0*l*c23*c12 + 2.0*g*c13*c22)
        self.K[4][9] = vol*0.5*(2.0*g*c32*c11 + 2.0*l*c12*c31)
        self.K[4][10] = vol*0.5*((2.0*g*c12 + l*c12)*c32 + (2.0*g*c32 + l*c32)*c12 + 2.0*g*c13*c33 + 2.0*g*c11*c31)
        self.K[4][11] = vol*0.5*(2.0*g*c13*c32 + 2.0*l*c33*c12)
    
        self.K[5][5] = vol*0.5*(2.0*(2.0*g*c13 + l*c13)*c13 + 2.0*g*c11*c11 + 2.0*g*c12*c12)
        self.K[5][6] = vol*0.5*(2.0*g*c23*c11 + 2.0*l*c13*c21)
        self.K[5][7] = vol*0.5*(2.0*l*c13*c22 + 2.0*g*c23*c12)
        self.K[5][8] = vol*0.5*((2.0*g*c13 + l*c13)*c23 + 2.0*g*c11*c21 + 2.0*g*c12*c22 + (2.0*g*c23 + l*c23)*c13)
        self.K[5][9] = vol*0.5*(2.0*l*c13*c31 + 2.0*g*c33*c11)
        self.K[5][10] = vol*0.5*(2.0*l*c13*c32 + 2.0*g*c33*c12)
        self.K[5][11] = vol*0.5*((2.0*g*c13 + l*c13)*c33 + 2.0*g*c11*c31 + 2.0*g*c12*c32 + (2.0*g*c33 + l*c33)*c13)
    
        self.K[6][6] = vol*0.5*(2.0*(2.0*g*c21 + l*c21)*c21 + 2.0*g*c23*c23 + 2.0*g*c22*c22)
        self.K[6][7] = vol*0.5*(2.0*l*c22*c21 + 2.0*g*c22*c21)
        self.K[6][8] = vol*0.5*(2.0*l*c23*c21 + 2.0*g*c23*c21)
        self.K[6][9] = vol*0.5*(2.0*g*c22*c32 + (2.0*g*c31 + l*c31)*c21 + 2.0*g*c23*c33 + (2.0*g*c21 + l*c21)*c31)
        self.K[6][10] = vol*0.5*(2.0*g*c22*c31 + 2.0*l*c32*c21)
        self.K[6][11] = vol*0.5*(2.0*g*c23*c31 + 2.0*l*c33*c21)
    
        self.K[7][7] = vol*0.5*(2.0*(2.0*g*c22 + l*c22)*c22 + 2.0*g*c21*c21 + 2.0*g*c23*c23)
        self.K[7][8] = vol*0.5*(2.0*l*c23*c22 + 2.0*g*c23*c22)
        self.K[7][9] = vol*0.5*(2.0*g*c32*c21 + 2.0*l*c22*c31)
        self.K[7][10] = vol*0.5*((2.0*g*c22 + l*c22)*c32 + (2.0*g*c32 + l*c32)*c22 + 2.0*g*c23*c33 + 2.0*g*c21*c31)
        self.K[7][11] = vol*0.5*(2.0*g*c23*c32 + 2.0*l*c33*c22)
    
        self.K[8][8] = vol*0.5*(2.0*(2.0*g*c23 + l*c23)*c23 + 2.0*g*c21*c21 + 2.0*g*c22*c22) 
        self.K[8][9] = vol*0.5*(2.0*l*c23*c31 + 2.0*g*c33*c21)
        self.K[8][10] = vol*0.5*(2.0*l*c23*c32 + 2.0*g*c33*c22)
        self.K[8][11] = vol*0.5*((2.0*g*c23 + l*c23)*c33 + 2.0*g*c21*c31 + 2.0*g*c22*c32 + (2.0*g*c33 + l*c33)*c23)
    
        self.K[9][9] = vol*0.5*(2.0*g*c32*c32 + 2*(2.0*g*c31 + l*c31)*c31 + 2.0*g*c33*c33)
        self.K[9][10] = vol*0.5*(2.0*l*c32*c31 + 2.0*g*c32*c31)
        self.K[9][11] = vol*0.5*(2.0*l*c33*c31 + 2.0*g*c33*c31)
    
        self.K[10][10] = vol*0.5*(2.0*(2.0*g*c32 + l*c32)*c32 + 2.0*g*c33*c33 + 2.0*g*c31*c31) 
        self.K[10][11] = vol*0.5*(2.0*g*c33*c32 + 2.0*l*c33*c32)
    
        self.K[11][11] = vol*0.5*(2.0*(2.0*g*c33 + l*c33)*c33 + 2.0*g*c31*c31 + 2.0*g*c32*c32)

        if len(self.vx) or len(self.vy) or len(self.vz):
            self.K[0][12] += 0.25*self.vx[0]*vol
            self.K[1][12] += 0.25*self.vy[0]*vol
            self.K[2][12] += 0.25*self.vz[0]*vol
            self.K[3][12] += 0.25*self.vx[1]*vol
            self.K[4][12] += 0.25*self.vy[1]*vol
            self.K[5][12] += 0.25*self.vz[1]*vol
            self.K[6][12] += 0.25*self.vx[2]*vol
            self.K[7][12] += 0.25*self.vy[2]*vol
            self.K[8][12] += 0.25*self.vz[2]*vol
            self.K[9][12] += 0.25*self.vx[3]*vol
            self.K[10][12] += 0.25*self.vy[3]*vol
            self.K[11][12] += 0.25*self.vz[3]*vol
        if not is_static:
            self.M = [
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            ]

            k1 = 0.1*vol
            k2 = 0.05*vol

            self.M[0][0] = k1
            self.M[0][3] = k2
            self.M[0][6] = k2
            self.M[0][9] = k2
    
            self.M[1][1] = k1
            self.M[1][4] = k2
            self.M[1][7] = k2
            self.M[1][10] = k2
    
            self.M[2][2] = k1
            self.M[2][5] = k2
            self.M[2][8] = k2
            self.M[2][11] = k2
    
            self.M[3][0] = k2
            self.M[3][3] = k1
            self.M[3][6] = k2
            self.M[3][9] = k2
    
            self.M[4][1] = k2
            self.M[4][4] = k1
            self.M[4][7] = k2
            self.M[4][10] = k2
    
            self.M[5][2] = k2
            self.M[5][5] = k1
            self.M[5][8] = k2
            self.M[5][11] = k2
    
            self.M[6][0] = k2
            self.M[6][3] = k2
            self.M[6][6] = k1
            self.M[6][9] = k2
    
            self.M[7][1] = k2
            self.M[7][4] = k2
            self.M[7][7] = k1
            self.M[7][10] = k2
    
            self.M[8][2] = k2
            self.M[8][5] = k2
            self.M[8][8] = k1
            self.M[8][11] = k2
    
            self.M[9][0] = k2
            self.M[9][3] = k2
            self.M[9][6] = k2
            self.M[9][9] = k1
    
            self.M[10][1] = k2
            self.M[10][4] = k2
            self.M[10][7] = k2
            self.M[10][10] = k1
    
            self.M[11][2] = k2
            self.M[11][5] = k2
            self.M[11][8] = k2
            self.M[11][11] = k1
            
            self.D = list(self.M)

            for i in range(0, len(self.M)):
                for j in range(0, len(self.M)):
                    self.M[i][j] *= self.density
                    self.D[i][j] *= self.damping


# Четырехузловой двумерный КЭ
class TFE2D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4
#        self.K = [[0]*9]*8
#        self.c = [[0]*4]*4
        self.K = [
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]
        self.c = [
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ]

    def volume(self):
        return math.sqrt((self.x[0] - self.x[1])**2 + (self.y[0] - self.y[1])**2)

    def __create__(self):
        if self.volume() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a = array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i]
                a[i][2] = self.y[i]
                a[i][3] = self.x[i]*self.y[i]
            b[j] = 1.0
            x = solve(a, b)
            self.c[j] = list(x)

    def generate(self, is_static=True):
        # Параметры квадратур Гаусса
        xi = [-0.57735027, -0.57735027,  0.57735027,  0.57735027]
        eta = [-0.57735027,  0.57735027, -0.57735027,  0.57735027]
        w = [1.0,  1.0,  1.0,  1.0]
        # Матрица упругих свойст
        d = array([
            [1.0 - self.m[0], self.m[0], 0.0],
            [self.m[0], 1.0 - self.m[0], 0.0],
            [0.0, 0.0, (1.0 - 2.0*self.m[0])/2.0]
            ])*self.e[0]/((1.0 + self.m[0])*(1.0 - 2.0*self.m[0]))
        # Формирование локальных матриц жесткости, масс и демпфирования
        local_k = zeros((8, 8))
        volume_load = zeros(8)
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(w)):
            for j in range(len(w)):
                # Изопараметрические функции формы и их производные
                shape = array([
                    (1.0 - xi[i])*(1.0 - eta[j])/4.0,
                    (1.0 + xi[i])*(1.0 - eta[j])/4.0,
                    (1.0 + xi[i])*(1.0 + eta[j])/4.0,
                    (1.0 - xi[i])*(1.0 + eta[j])/4.0
                    ])
                shape_dxi = array([
                    -(1.0 - eta[j])/4.0,
                    (1.0 - eta[j])/4.0,
                    (1.0 + eta[j])/4.0,
                    -(1.0 + eta[j])/4.0
                    ])
                shape_deta = array([
                    -(1.0 - xi[i])/4.0,
                    -(1.0 + xi[i])/4.0,
                     (1.0 + xi[i])/4.0,
                     (1.0 - xi[i])/4.0
                    ])
                # Матрица Якоби
                jacobi = array([
                    [sum(shape_dxi*self.x), sum(shape_dxi*self.y)],
                    [sum(shape_deta*self.x), sum(shape_deta*self.y)]
                    ])
                # Якобиан
                jacobian = det(jacobi)
                inverted_jacobi = inv(jacobi)
                shape_dx = inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta
                shape_dy = inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta
                # Матрица градиентов
                b = array([
                    [shape_dx[0], 0.0, shape_dx[1], 0.0, shape_dx[2], 0.0, shape_dx[3], 0.0],
                    [0.0, shape_dy[0], 0.0, shape_dy[1], 0.0, shape_dy[2], 0.0, shape_dy[3]],
                    [shape_dy[0], shape_dx[0], shape_dy[1], shape_dx[1], shape_dy[2], shape_dx[2], shape_dy[3], shape_dx[3]]
                    ])
                bt = b.conj().transpose()
                local_k += bt.dot(d).dot(b)*jacobian*w[i]
                # Учет объемной нагрузки
                if len(self.vx) or len(self.vy):
                    for k in range (0, 4):
                        volume_load[2*k] = self.vx[k]*shape[k]
                        volume_load[2*k + 1] = self.vy[k]*shape[k]

        for i in range(0, 8):
            for j in range(i, 8):
                self.K[i][j] = local_k[i][j]
            self.K[i][8] = volume_load[i]
#        import sys
#        print('*******************************')
#        for i in range(0, len(self.K)):
#            for j in range(0, len(self.K[0])):
#                sys.stdout.write('%f' % self.K[i][j])
#                sys.stdout.write(' ')
#            sys.stdout.write('\n')
#        print('*******************************')

    def calc(self, u, index):
        l = self.e[0]/(1.0 - self.m[0]*self.m[0])
        g = self.e[0]/(2.0 + 2.0*self.m[0])
        u0 = u[0]
        v0 = u[1]
        u1 = u[2]
        v1 = u[3]
        u2 = u[4]
        v2 = u[5]
        u3 = u[6]
        v3 = u[7]

        for i in range(0, self.size):
            if index == 0:      # Exx
                u[i] = u0*self.c[0][1] + u0*self.c[0][3]*self.y[i] + u1*self.c[1][1] + u1*self.c[1][3]*self.y[i] + \
                       u2*self.c[2][1] + u2*self.c[2][3]*self.y[i] + u3*self.c[3][1] + u3*self.c[3][3]*self.y[i]
            elif index == 1:    # Eyy
                u[i] = v0*self.c[0][2] + v0*self.c[0][3]*self.x[i] + v1*self.c[1][2] + v1*self.c[1][3]*self.x[i] + \
                       v2*self.c[2][2] + v2*self.c[2][3]*self.x[i] + v3*self.c[3][2] + v3*self.c[3][3]*self.x[i]
            elif index == 2:    # Exy
                u[i] = u0*self.c[0][2] + u0*self.c[0][3]*self.x[i] + u1*self.c[1][2] + u1*self.c[1][3]*self.x[i] + \
                       u2*self.c[2][2] + u2*self.c[2][3]*self.x[i] + u3*self.c[3][2] + u3*self.c[3][3]*self.x[i] + \
                       v0*self.c[0][1] + v0*self.c[0][3]*self.y[i] + v1*self.c[1][1] + v1*self.c[1][3]*self.y[i] + \
                       v2*self.c[2][1] + v2*self.c[2][3]*self.y[i] + v3*self.c[3][1] + v3*self.c[3][3]*self.y[i]
            elif index == 3:    # Sxx
                u[i] = l*(u0*(self.c[0][1] + self.c[0][3]*self.y[i]) + u1*(self.c[1][1] + self.c[1][3]*self.y[i]) +
                          u2*(self.c[2][1] + self.c[2][3]*self.y[i]) + u3*(self.c[3][1] + self.c[3][3]*self.y[i]) +
                          self.m[0]*(v0*(self.c[0][2] + self.c[0][3]*self.x[i]) + v1*(self.c[1][2] +
                          self.c[1][3]*self.x[i]) + v2*(self.c[2][2] + self.c[2][3]*self.x[i]) + v3*(self.c[3][2] +
                          self.c[3][3]*self.x[i])))
            elif index == 4:    # Syy
                u[i] = l*(self.m[0]*(u0*(self.c[0][1] + self.c[0][3]*self.y[i]) + u1*(self.c[1][1] +
                                         self.c[1][3]*self.y[i]) + u2*(self.c[2][1] + self.c[2][3]*self.y[i]) +
                                     u3*(self.c[3][1] + self.c[3][3]*self.y[i])) +
                          v0*(self.c[0][2] + self.c[0][3]*self.x[i]) + v1*(self.c[1][2] + self.c[1][3]*self.x[i]) +
                          v2*(self.c[2][2] + self.c[2][3]*self.x[i]) + v3*(self.c[3][2] + self.c[3][3]*self.x[i]))
            elif index == 5:    # Sxy
                u[i] = g*(u0*(self.c[0][2] + self.c[0][3]*self.x[i]) + u1*(self.c[1][2] + self.c[1][3]*self.x[i]) +
                          u2*(self.c[2][2] + self.c[2][3]*self.x[i]) + u3*(self.c[3][2] + self.c[3][3]*self.x[i]) +
                          v0*(self.c[0][1] + self.c[0][3]*self.y[i]) + v1*(self.c[1][1] + self.c[1][3]*self.y[i]) +
                          v2*(self.c[2][1] + self.c[2][3]*self.y[i]) + v3*(self.c[3][1] + self.c[3][3]*self.y[i]))
