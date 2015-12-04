#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#           Реализация функций форм конечных элементов
###################################################################

import math
from abc import abstractmethod
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
    def set_volume_force(self, *args):
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
        self.K[0][2] = self.K[1][2] = 0.0
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
        if index == 0:      # Exx
            res = u0*self.c[0][1] + u1*self.c[1][1] + u2*self.c[2][1]
        elif index == 1:    # Eyy
            res = v0*self.c[0][2] + v1*self.c[1][2] + v2*self.c[2][2]
        elif index == 2:    # Exy
            res = u0*self.c[0][2] + u1*self.c[1][2] + u2*self.c[2][2] + v0*self.c[0][1] + v1*self.c[1][1] + \
                  v2*self.c[2][1]
        elif index == 3:    # Sxx
            res = self.e[0]/(1.0 - self.m[0]*self.m[0])*((u0*self.c[0][1] + u1*self.c[1][1] + u2*self.c[2][1]) +
                                   self.m[0]*(v0*self.c[0][2] + v1*self.c[1][2] + v2*self.c[2][2]))
        elif index == 4:    # Syy
            res = self.e[0]/(1.0 - self.m[0]*self.m[0])*(self.m[0]*(u0*self.c[0][1] + u1*self.c[1][1] +
                                u2*self.c[2][1]) + (v0*self.c[0][2] + v1*self.c[1][2] + v2*self.c[2][2]))
        elif index == 5:    # Sxy
            res = self.e[0]/(2.0 + 2.0*self.m[0])*(u0*self.c[0][2] + u1*self.c[1][2] + u2*self.c[2][2] + v0*self.c[0][1] + v1*self.c[1][1] + v2*self.c[2][1])
        u[0] = u[1] = u[2] = res

    def generate(self, is_static=True):
        k = self.e[0]/(1.0 - self.m[0]*self.m[0])
        g = self.e[0]/(2.0 + 2.0*self.m[0])
        volume = self.volume()

        self.K[0][0] = volume*(k*self.c[0][1]*self.c[0][1] + g*self.c[0][2]*self.c[0][2])
        self.K[0][1] = volume*(k*self.c[0][1]*self.m[0]*self.c[0][2] + g*self.c[0][2]*self.c[0][1])
        self.K[0][2] = volume*(k*self.c[0][1]*self.c[1][1] + g*self.c[0][2]*self.c[1][2])
        self.K[0][3] = volume*(k*self.c[0][1]*self.m[0]*self.c[1][2] + g*self.c[0][2]*self.c[1][1])
        self.K[0][4] = volume*(k*self.c[0][1]*self.c[2][1] + g*self.c[0][2]*self.c[2][2])
        self.K[0][5] = volume*(g*self.c[0][2]*self.c[2][1] + k*self.c[0][1]*self.m[0]*self.c[2][2])

        self.K[1][1] = volume*(k*self.c[0][2]*self.c[0][2] + g*self.c[0][1]*self.c[0][1])
        self.K[1][2] = volume*(g*self.c[1][2]*self.c[0][1] + k*self.c[1][1]*self.m[0]*self.c[0][2])
        self.K[1][3] = volume*(g*self.c[0][1]*self.c[1][1] + k*self.c[0][2]*self.c[1][2])
        self.K[1][4] = volume*(k*self.c[2][1]*self.m[0]*self.c[0][2] + g*self.c[2][2]*self.c[0][1])
        self.K[1][5] = volume*(g*self.c[0][1]*self.c[2][1] + k*self.c[0][2]*self.c[2][2])

        self.K[2][2] = volume*(k*self.c[1][1]*self.c[1][1] + g*self.c[1][2]*self.c[1][2])
        self.K[2][3] = volume*(g*self.c[1][2]*self.c[1][1] + k*self.c[1][1]*self.m[0]*self.c[1][2])
        self.K[2][4] = volume*(g*self.c[1][2]*self.c[2][2] + k*self.c[1][1]*self.c[2][1])
        self.K[2][5] = volume*(g*self.c[1][2]*self.c[2][1] + k*self.c[1][1]*self.m[0]*self.c[2][2])

        self.K[3][3] = volume*(g*self.c[1][1]*self.c[1][1] + k*self.c[1][2]*self.c[1][2])
        self.K[3][4] = volume*(g*self.c[2][2]*self.c[1][1] + k*self.c[2][1]*self.m[0]*self.c[1][2])
        self.K[3][5] = volume*(k*self.c[1][2]*self.c[2][2] + g*self.c[1][1]*self.c[2][1])

        self.K[4][4] = volume*(g*self.c[2][2]*self.c[2][2] + k*self.c[2][1]*self.c[2][1])
        self.K[4][5] = volume*(g*self.c[2][2]*self.c[2][1] + k*self.c[2][1]*self.m[0]*self.c[2][2])

        self.K[5][5] = volume*(g*self.c[2][1]*self.c[2][1] + k*self.c[2][2]*self.c[2][2])

        # Вычисление интеграла для объемных сил
        self.K[0][6] = self.K[1][6] = self.K[2][6] = self.K[3][6] = self.K[4][6] = self.K[5][6] = 0
        if len(self.vx) or len(self.vy):
            self.K[0][6] += self.vx[0]*volume/6.0
            self.K[2][6] += self.vx[1]*volume/6.0
            self.K[4][6] += self.vx[2]*volume/6.0
            self.K[1][6] += self.vy[0]*volume/6.0
            self.K[3][6] += self.vy[1]*volume/6.0
            self.K[5][6] += self.vy[2]*volume/6.0

        if not is_static:
            # Формирование матриц массы и демпфирования
            k00 = volume/12.0
            k01 = volume/6.0
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

            self.D = self.M
            for i in range(0, len(self.M)):
                for j in range(0, len(self.M)):
                    self.M[i][j] *= self.density
                    self.D[i][j] *= self.damping
