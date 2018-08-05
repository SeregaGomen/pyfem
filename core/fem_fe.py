#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#           Реализация функций форм конечных элементов
###################################################################

import math
from abc import abstractmethod
from numpy.linalg import solve, LinAlgError
from numpy import array
from numpy import zeros
from numpy.linalg import det
from numpy.linalg import inv
from core.fem_error import TFEMException
from core.fem_defs import eps


# Абстрактный базовый класс, описывающий конечный элемент (КЭ)
class TFE:
    def __init__(self):
        self.size = 0
        self.h = 1                          # Толщина (для оболочек и пластин)
        self.density = 0                    # Плотность
        self.damping = [0, 0]               # Параметры демпфирования
        self.e = []                         # Модуль Юнга
        self.m = []                         # Коэффициент Пуассона
        self.x = []                         # Координаты вершин КЭ
        self.y = []
        self.z = []
        self.K = [[]]                       # Локальная матрица жесткости
        self.M = [[]]                       # ... масс
        self.C = [[]]                       # ... демпфирования
        self.a = [[]]                       # Коэффициенты функций форм

    # Задание толщины
    def set_h(self, h):
        self.h = h

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

    # Формирование матриц жесткости, масс и демпфирования
    def generate(self, is_static=True):
        self.K = self._generate_stiffness_matrix_()
        if not is_static:
            self.M, self.C = self._generate_mass_damping_matrix_()

    @abstractmethod
    # Вычисление функций форм КЭ
    def __create__(self):
        raise NotImplementedError('Method TFE.__create__ is pure virtual')

    # Вычисления стандартных результатов КЭ
    @abstractmethod
    def calc(self, u):
        return [[]]

    # Вычисление матрицы жесткости
    @abstractmethod
    def _generate_stiffness_matrix_(self):
        return [[]]

    # Вычисление матриц масс и демпфирования
    @abstractmethod
    def _generate_mass_damping_matrix_(self):
        return [[]], [[]]


# Линейный (двухузловой) одномерный КЭ
class TFE1D2(TFE):
    def __init__(self):
        super().__init__()
        self.size = 2

    def __length__(self):
        return math.fabs(self.x[1] - self.x[0])

    def __create__(self):
        if self.__length__() == 0.0:
            raise TFEMException('incorrect_fe_err')
        self.a = zeros((2, 2))
        self.a[0][0] = self.x[1] / (self.x[1] - self.x[0])
        self.a[0][1] = -1.0 / (self.x[1] - self.x[0])
        self.a[1][0] = self.x[0] / (self.x[0] - self.x[1])
        self.a[1][1] = -1.0 / (self.x[0] - self.x[1])

    def calc(self, u):
        res = zeros((self.size, self.size))
        res[0][0] = res[0][1] = u[0] * self.a[0][1] + u[1] * self.a[1][1]
        res[1][0] = res[1][1] = self.e[0]*(u[0] * self.a[0][1] + u[1] * self.a[1][1])
        return res

    def _generate_stiffness_matrix_(self):
        return array([
                    [1.0, -1.0],
                    [-1.0, .0]
                ])*self.e[0]/self.__length__()

    # Формирование матрицы демпфирования по Релею
    def _generate_mass_damping_matrix_(self):
        a = array([
                [2.0, 1.0],
                [1.0, 2.0]
            ])
        m = a*self.__length__()/6.0*self.density
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Линейный (трехузловой) треугольный КЭ
class TFE2D3(TFE):
    def __init__(self):
        super().__init__()
        self.size = 3

    def __square__(self):
        a = math.sqrt((self.x[0] - self.x[1])*(self.x[0] - self.x[1]) + (self.y[0] - self.y[1])*(self.y[0] - self.y[1]))
        b = math.sqrt((self.x[0] - self.x[2])*(self.x[0] - self.x[2]) + (self.y[0] - self.y[2])*(self.y[0] - self.y[2]))
        c = math.sqrt((self.x[2] - self.x[1])*(self.x[2] - self.x[1]) + (self.y[2] - self.y[1])*(self.y[2] - self.y[1]))
        p = 0.5*(a + b + c)
        return math.sqrt(p*(p - a)*(p - b)*(p - c))

    def __create__(self):
        det0 = self.y[2]*self.x[1] - self.y[2]*self.x[0] - self.y[0]*self.x[1] - self.y[1]*self.x[2] + \
               self.y[1]*self.x[0] + self.y[0]*self.x[2]
        if math.fabs(det0) < eps:
            raise TFEMException('incorrect_fe_err')
        index = [[2, 1], [0, 2], [1, 0]]
        self.a = zeros((self.size, self.size))
        for i in range(0, 3):
            det1 = self.y[index[i][0]]*self.x[index[i][1]] - self.y[index[i][1]]*self.x[index[i][0]]
            det2 = self.y[index[i][1]] - self.y[index[i][0]]
            det3 = self.x[index[i][0]] - self.x[index[i][1]]
            self.a[i][0] = det1 / det0
            self.a[i][1] = det2 / det0
            self.a[i][2] = det3 / det0

    def calc(self, u):
        m = self.m[0]
        g = self.e[0]/(2.0 + 2.0*m)
        k = self.e[0]/(1.0 - m**2)
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        res = zeros((6, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = self.a[j][1]
                dy[i][j] = self.a[j][2]
                res[0][i] += u[2*j]*dx[i][j]
                res[1][i] += u[2*j + 1]*dy[i][j]
                res[2][i] += u[2*j]*dy[i][j] + u[2*j + 1]*dx[i][j]
                res[3][i] += k*(u[2*j]*dx[i][j] + m*u[2*j + 1]*dy[i][j])
                res[4][i] += k*(u[2*j + 1]*dy[i][j] + m*u[2*j]*dx[i][j])
                res[5][i] += g*(u[2*j]*dy[i][j] + u[2*j + 1]*dx[i][j])
        return res

    # Формирование локальной матриц жесткости
    def _generate_stiffness_matrix_(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
            ])*self.e[0]/(1.0 - self.m[0]**2)
        # Производные функций формы
        shape_dx = array([self.a[0][1], self.a[1][1], self.a[2][1]])
        shape_dy = array([self.a[0][2], self.a[1][2], self.a[2][2]])
        # Матрица градиентов
        b = array([
            [shape_dx[0], 0.0, shape_dx[1], 0.0, shape_dx[2], 0.0],
            [0.0, shape_dy[0], 0.0, shape_dy[1], 0.0, shape_dy[2]],
            [shape_dy[0], shape_dx[0], shape_dy[1], shape_dx[1], shape_dy[2], shape_dx[2]]
        ])
        return b.conj().transpose().dot(d).dot(b)*self.__square__()

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        a = array([
                [0.5, 0.0, 0.25, 0.0, 0.25, 0.0],
                [0.0, 0.5, 0.0, 0.25, 0.0, 0.25],
                [0.25, 0.0, 0.5, 0.0, 0.25, 0.0],
                [0.0, 0.25, 0.0, 0.5, 0.0, 0.25],
                [0.25, 0.0, 0.25, 0.0, 0.5, 0.0],
                [0.0, 0.25, 0.0, 0.25, 0.0, 0.5]
        ])
        m = a*self.density*self.__square__()
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Линейный (четырехузловой) тетраэдральный КЭ
class TFE3D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4

    def __volume__(self):
        a = (self.x[1] - self.x[0])*(self.y[2] - self.y[0])*(self.z[3] - self.z[0]) + \
            (self.x[3] - self.x[0])*(self.y[1] - self.y[0])*(self.z[2] - self.z[0]) + \
            (self.x[2] - self.x[0])*(self.y[3] - self.y[0])*(self.z[1] - self.z[0])
        b = (self.x[3] - self.x[0])*(self.y[2] - self.y[0])*(self.z[1] - self.z[0]) + \
            (self.x[2] - self.x[0])*(self.y[1] - self.y[0])*(self.z[3] - self.z[0]) + \
            (self.x[1] - self.x[0])*(self.y[3] - self.y[0])*(self.z[2] - self.z[0])
        return math.fabs(a - b)/6.0

    def __create__(self):
        if self.__volume__() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a, self.a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i]
                a[i][2] = self.y[i]
                a[i][3] = self.z[i]
            b[j] = 1.0
            x = solve(a, b)
            self.a[j] = list(x)

    def calc(self, u):
        g = self.e[0]/(2.0 + 2.0*self.m[0])
        l = 2.0*self.m[0]*g/(1.0 - 2.0*self.m[0])
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        dz = zeros((self.size, self.size))
        res = zeros((12, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = self.a[j][1]
                dy[i][j] = self.a[j][2]
                dz[i][j] = self.a[j][3]
                res[0][i] += u[3*j]*dx[i][j]
                res[1][i] += u[3*j + 1]*dy[i][j]
                res[2][i] += u[3*j + 2]*dz[i][j]
                res[3][i] += u[3*j]*dy[i][j] + u[3*j + 1]*dx[i][j]
                res[4][i] += u[3*j]*dz[i][j] + u[3*j + 2]*dx[i][j]
                res[5][i] += u[3*j + 1]*dz[i][j] + u[3*j + 2]*dy[i][j]
                res[6][i] += 2.0*g*u[3*j]*dx[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[7][i] += 2.0*g*u[3*j + 1]*dy[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[8][i] += 2.0*g*u[3*j + 2]*dz[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[9][i] += g*(u[3*j]*dy[i][j] + u[3*j + 1]*dx[i][j])
                res[10][i] += g*(u[3*j]*dz[i][j] + u[3*j + 2]*dx[i][j])
                res[11][i] += g*(u[3*j + 1]*dz[i][j] + u[3*j + 2]*dy[i][j])
        return res

    # Формирование локальной матрицы жесткости
    def _generate_stiffness_matrix_(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), 1.0, self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0])],
            ])*self.e[0]*(1.0 - self.m[0])/(1.0 + self.m[0])/(1.0 - 2.0*self.m[0])
        # Производные функций формы
        shape_dx = array([self.a[0][1], self.a[1][1], self.a[2][1], self.a[3][1]])
        shape_dy = array([self.a[0][2], self.a[1][2], self.a[2][2], self.a[3][2]])
        shape_dz = array([self.a[0][3], self.a[1][3], self.a[2][3], self.a[3][3]])
        # Матрица градиентов
        b = array([
            [shape_dx[0], 0.0, 0.0, shape_dx[1], 0.0, 0.0, shape_dx[2], 0.0, 0.0, shape_dx[3], 0.0, 0.0],
            [0.0, shape_dy[0], 0.0, 0.0, shape_dy[1], 0.0, 0.0, shape_dy[2], 0.0, 0.0, shape_dy[3], 0.0],
            [0.0, 0.0, shape_dz[0], 0.0, 0.0, shape_dz[1], 0.0, 0.0, shape_dz[2], 0.0, 0.0, shape_dz[3]],
            [shape_dy[0], shape_dx[0], 0.0, shape_dy[1], shape_dx[1], 0.0, shape_dy[2], shape_dx[2], 0.0, shape_dy[3],
             shape_dx[3], 0.0],
            [0.0, shape_dz[0], shape_dy[0], 0.0, shape_dz[1], shape_dy[1], 0.0, shape_dz[2], shape_dy[2], 0.0,
             shape_dz[3], shape_dy[3]],
            [shape_dz[0], 0.0, shape_dx[0], shape_dz[1], 0.0, shape_dx[1], shape_dz[2], 0.0, shape_dx[2], shape_dz[3],
             0.0, shape_dx[3]]
        ])
        return b.conj().transpose().dot(d).dot(b)*self.__volume__()

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        a = array([
                [1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0],
                [0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5],
                [0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0],
                [0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0],
                [0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5],
                [0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0],
                [0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0],
                [0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5],
                [0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0]
        ])*0.1
        m = a*self.density*self.__volume__()
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Билинейный четырехузловой двумерный КЭ
class TFE2D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4
        # Параметры квадратур Гаусса
        self.__xi__ = [-0.57735027, -0.57735027, 0.57735027, 0.57735027]
        self.__eta__ = [-0.57735027, 0.57735027, -0.57735027, 0.57735027]
        self.__w__ = [1.0, 1.0, 1.0, 1.0]

    def __square__(self):
        return math.sqrt((self.x[0] - self.x[1])**2 + (self.y[0] - self.y[1])**2)

    def __create__(self):
        if self.__square__() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a, self.a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i]
                a[i][2] = self.y[i]
                a[i][3] = self.x[i]*self.y[i]
            b[j] = 1.0
            x = solve(a, b)
            self.a[j] = list(x)

    def calc(self, u):
        m = self.m[0]
        g = self.e[0]/(2.0 + 2.0*m)
        k = self.e[0]/(1.0 - m**2)
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        res = zeros((6, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = self.a[j][1] + self.a[j][3] * self.y[i]
                dy[i][j] = self.a[j][2] + self.a[j][3] * self.x[i]
                res[0][i] += u[2*j]*dx[i][j]
                res[1][i] += u[2*j + 1]*dy[i][j]
                res[2][i] += u[2*j]*dy[i][j] + u[2*j + 1]*dx[i][j]
                res[3][i] += k*(u[2*j]*dx[i][j] + m*u[2*j + 1]*dy[i][j])
                res[4][i] += k*(u[2*j + 1]*dy[i][j] + m*u[2*j]*dx[i][j])
                res[5][i] += g*(u[2*j]*dy[i][j] + u[2*j + 1]*dx[i][j])
        return res

    # Изопараметрические функции формы и их производные
    def __shape__(self, i):
        return array([
                0.25*(1.0 - self.__xi__[i])*(1.0 - self.__eta__[i]),
                0.25*(1.0 + self.__xi__[i])*(1.0 - self.__eta__[i]),
                0.25*(1.0 + self.__xi__[i])*(1.0 + self.__eta__[i]),
                0.25*(1.0 - self.__xi__[i])*(1.0 + self.__eta__[i])
                ])

    def __shape_dxi__(self, i):
        return array([
                -0.25*(1.0 - self.__eta__[i]),
                0.25*(1.0 - self.__eta__[i]),
                0.25*(1.0 + self.__eta__[i]),
                -0.25*(1.0 + self.__eta__[i])
                ])
    
    def __shape_deta__(self, i):
        return array([
                -0.25*(1.0 - self.__xi__[i]),
                -0.25*(1.0 + self.__xi__[i]),
                0.25*(1.0 + self.__xi__[i]),
                0.25*(1.0 - self.__xi__[i])
                ])

    # Формирование локальной матрицы жесткости
    def _generate_stiffness_matrix_(self):
        local_k = zeros((8, 8))
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
            ])*self.e[0]/(1.0 - self.m[0]**2)
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
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
            local_k += b.conj().transpose().dot(d).dot(b)*jacobian*self.__w__[i]
        return local_k

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        local_m = zeros((8, 8))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x), sum(shape_dxi*self.y)],
                [sum(shape_deta*self.x), sum(shape_deta*self.y)]
                ])
            # Якобиан
            jacobian = det(jacobi)
            # Вспомогательная матрица для построения матриц масс и демпфирования
            c = array([
                [shape[0], 0.0, shape[1], 0.0, shape[2], 0.0, shape[3], 0.0],
                [0.0, shape[0], 0.0, shape[1], 0.0, shape[2], 0.0, shape[3]]
                ])
            local_m += c.conj().transpose().dot(c)*jacobian*self.__w__[i]
        m = self.density*local_m
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Восьмиузловой призматический КЭ
class TFE3D8(TFE):
    def __init__(self):
        super().__init__()
        self.size = 8
        # Параметры квадратур Гаусса
        self.__xi__ = [-0.57735027, -0.57735027, -0.57735027, -0.57735027,
                       0.57735027, 0.57735027, 0.57735027, 0.57735027]
        self.__eta__ = [-0.57735027, -0.57735027, 0.57735027, 0.57735027,
                        -0.57735027, -0.57735027, 0.57735027, 0.57735027]
        self.__psi__ = [-0.57735027, 0.57735027, -0.57735027, 0.57735027,
                        -0.57735027, 0.57735027, -0.57735027, 0.57735027]
        self.__w__ = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    def __create__(self):
        self.a, a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0]*self.size)
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i]
                a[i][2] = self.y[i]
                a[i][3] = self.z[i]
                a[i][4] = self.x[i]*self.y[i]
                a[i][5] = self.x[i]*self.z[i]
                a[i][6] = self.y[i]*self.z[i]
                a[i][7] = self.x[i]*self.y[i]*self.z[i]
            b[j] = 1.0
            try:
                x = solve(a, b)
            except LinAlgError:
                raise TFEMException('incorrect_fe_err')
            self.a[j] = list(x)

    def calc(self, u):
        g = self.e[0]/(2.0 + 2.0*self.m[0])
        l = 2.0*self.m[0]*g/(1.0 - 2.0*self.m[0])
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        dz = zeros((self.size, self.size))
        res = zeros((12, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = (self.a[j][1] + self.a[j][4] * self.y[i]) + \
                           (self.a[j][5] * self.z[i] + self.a[j][7] * self.y[i] * self.z[i])
                dy[i][j] = (self.a[j][2] + self.a[j][4] * self.x[i]) + \
                           (self.a[j][6] * self.z[i] + self.a[j][7] * self.x[i] * self.z[i])
                dz[i][j] = (self.a[j][3] + self.a[j][5] * self.x[i]) + \
                           (self.a[j][6] * self.y[i] + self.a[j][7] * self.x[i] * self.y[i])
                res[0][i] += u[3*j]*dx[i][j]
                res[1][i] += u[3*j + 1]*dy[i][j]
                res[2][i] += u[3*j + 2]*dz[i][j]
                res[3][i] += u[3*j]*dy[i][j] + u[3*j + 1]*dx[i][j]
                res[4][i] += u[3*j]*dz[i][j] + u[3*j + 2]*dx[i][j]
                res[5][i] += u[3*j + 1]*dz[i][j] + u[3*j + 2]*dy[i][j]
                res[6][i] += 2.0*g*u[3*j]*dx[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[7][i] += 2.0*g*u[3*j + 1]*dy[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[8][i] += 2.0*g*u[3*j + 2]*dz[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[9][i] += g*(u[3*j]*dy[i][j] + u[3*j + 1]*dx[i][j])
                res[10][i] += g*(u[3*j]*dz[i][j] + u[3*j + 2]*dx[i][j])
                res[11][i] += g*(u[3*j + 1]*dz[i][j] + u[3*j + 2]*dy[i][j])
        return res

    # Изопараметрические функции формы и их производные
    def __shape__(self, i): 
        return array([
            0.125 * (1.0 - self.__xi__[i]) * (1.0 - self.__eta__[i]) * (1.0 - self.__psi__[i]),
            0.125 * (1.0 + self.__xi__[i]) * (1.0 - self.__eta__[i]) * (1.0 - self.__psi__[i]),
            0.125 * (1.0 + self.__xi__[i]) * (1.0 + self.__eta__[i]) * (1.0 - self.__psi__[i]),
            0.125 * (1.0 - self.__xi__[i]) * (1.0 + self.__eta__[i]) * (1.0 - self.__psi__[i]),
            0.125 * (1.0 - self.__xi__[i]) * (1.0 - self.__eta__[i]) * (1.0 + self.__psi__[i]),
            0.125 * (1.0 + self.__xi__[i]) * (1.0 - self.__eta__[i]) * (1.0 + self.__psi__[i]),
            0.125 * (1.0 + self.__xi__[i]) * (1.0 + self.__eta__[i]) * (1.0 + self.__psi__[i]),
            0.125 * (1.0 - self.__xi__[i]) * (1.0 + self.__eta__[i]) * (1.0 + self.__psi__[i])
        ])
    
    def __shape_dxi__(self, i):
        return array([
            -0.125 * (1.0 - self.__eta__[i]) * (1.0 - self.__psi__[i]),
            0.125 * (1.0 - self.__eta__[i]) * (1.0 - self.__psi__[i]),
            0.125 * (1.0 + self.__eta__[i]) * (1.0 - self.__psi__[i]),
            -0.125 * (1.0 + self.__eta__[i]) * (1.0 - self.__psi__[i]),
            -0.125 * (1.0 - self.__eta__[i]) * (1.0 + self.__psi__[i]),
            0.125 * (1.0 - self.__eta__[i]) * (1.0 + self.__psi__[i]),
            0.125 * (1.0 + self.__eta__[i]) * (1.0 + self.__psi__[i]),
            -0.125 * (1.0 + self.__eta__[i]) * (1.0 + self.__psi__[i])
        ])
    
    def __shape_deta__(self, i): 
        return array([
            -0.125 * (1.0 - self.__xi__[i]) * (1.0 - self.__psi__[i]),
            -0.125 * (1.0 + self.__xi__[i]) * (1.0 - self.__psi__[i]),
            0.125 * (1.0 + self.__xi__[i]) * (1.0 - self.__psi__[i]),
            0.125 * (1.0 - self.__xi__[i]) * (1.0 - self.__psi__[i]),
            -0.125 * (1.0 - self.__xi__[i]) * (1.0 + self.__psi__[i]),
            -0.125 * (1.0 + self.__xi__[i]) * (1.0 + self.__psi__[i]),
            0.125 * (1.0 + self.__xi__[i]) * (1.0 + self.__psi__[i]),
            0.125 * (1.0 - self.__xi__[i]) * (1.0 + self.__psi__[i])
        ])
    
    def __shape_dpsi__(self, i): 
        return array([
                -0.125 * (1.0 - self.__xi__[i]) * (1.0 - self.__eta__[i]),
                -0.125 * (1.0 + self.__xi__[i]) * (1.0 - self.__eta__[i]),
                -0.125 * (1.0 + self.__xi__[i]) * (1.0 + self.__eta__[i]),
                -0.125 * (1.0 - self.__xi__[i]) * (1.0 + self.__eta__[i]),
                0.125 * (1.0 - self.__xi__[i]) * (1.0 - self.__eta__[i]),
                0.125 * (1.0 + self.__xi__[i]) * (1.0 - self.__eta__[i]),
                0.125 * (1.0 + self.__xi__[i]) * (1.0 + self.__eta__[i]),
                0.125 * (1.0 - self.__xi__[i]) * (1.0 + self.__eta__[i])
            ])

    # Формирование локальной ыматрицы жесткости
    def _generate_stiffness_matrix_(self):
        local_k = zeros((3*self.size, 3*self.size))
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), 1.0, self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0])],
            ])*self.e[0]*(1.0 - self.m[0])/(1.0 + self.m[0])/(1.0 - 2.0*self.m[0])
        # Формирование локальных матриц жесткости, масс и демпфирования
        # Интегрирование по кубу [-1; 1] x [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
            shape_dpsi = self.__shape_dpsi__(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x), sum(shape_dxi*self.y), sum(shape_dxi*self.z)],
                [sum(shape_deta*self.x), sum(shape_deta*self.y), sum(shape_deta*self.z)],
                [sum(shape_dpsi*self.x), sum(shape_dpsi*self.y), sum(shape_dpsi*self.z)]
                ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = (inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta) + \
                       (inverted_jacobi[0, 2]*shape_dpsi)
            shape_dy = (inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta) + \
                       (inverted_jacobi[1, 2]*shape_dpsi)
            shape_dz = (inverted_jacobi[2, 0]*shape_dxi + inverted_jacobi[2, 1]*shape_deta) + \
                       (inverted_jacobi[2, 2]*shape_dpsi)
            # Матрица градиентов
            b = array([
                [shape_dx[0], 0.0, 0.0, shape_dx[1], 0.0, 0.0, shape_dx[2], 0.0, 0.0, shape_dx[3], 0.0, 0.0,
                 shape_dx[4], 0.0, 0.0, shape_dx[5], 0.0, 0.0, shape_dx[6], 0.0, 0.0, shape_dx[7], 0.0, 0.0],
                [0.0, shape_dy[0], 0.0, 0.0, shape_dy[1], 0.0, 0.0, shape_dy[2], 0.0, 0.0, shape_dy[3], 0.0, 0.0,
                 shape_dy[4], 0.0, 0.0, shape_dy[5], 0.0, 0.0, shape_dy[6], 0.0, 0.0, shape_dy[7], 0.0],
                [0.0, 0.0, shape_dz[0], 0.0, 0.0, shape_dz[1], 0.0, 0.0, shape_dz[2], 0.0, 0.0, shape_dz[3], 0.0, 0.0,
                 shape_dz[4], 0.0, 0.0, shape_dz[5], 0.0, 0.0, shape_dz[6], 0.0, 0.0, shape_dz[7]],
                [shape_dy[0], shape_dx[0], 0.0, shape_dy[1], shape_dx[1], 0.0, shape_dy[2], shape_dx[2], 0.0,
                 shape_dy[3], shape_dx[3], 0.0, shape_dy[4], shape_dx[4], 0.0, shape_dy[5], shape_dx[5], 0.0,
                 shape_dy[6], shape_dx[6], 0.0, shape_dy[7], shape_dx[7], 0.0],
                [0.0, shape_dz[0], shape_dy[0], 0.0, shape_dz[1], shape_dy[1], 0.0, shape_dz[2], shape_dy[2], 0.0,
                 shape_dz[3], shape_dy[3], 0.0, shape_dz[4], shape_dy[4], 0.0, shape_dz[5], shape_dy[5], 0.0,
                 shape_dz[6], shape_dy[6], 0.0, shape_dz[7], shape_dy[7]],
                [shape_dz[0], 0.0, shape_dx[0], shape_dz[1], 0.0, shape_dx[1], shape_dz[2], 0.0, shape_dx[2],
                 shape_dz[3], 0.0, shape_dx[3], shape_dz[4], 0.0, shape_dx[4], shape_dz[5], 0.0, shape_dx[5],
                 shape_dz[6], 0.0, shape_dx[6], shape_dz[7], 0.0, shape_dx[7]],
                ])
            local_k += b.conj().transpose().dot(d).dot(b)*jacobian*self.__w__[i]
        return local_k

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        local_m = zeros((3 * self.size, 3 * self.size))
        # Интегрирование по кубу [-1; 1] x [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
            shape_dpsi = self.__shape_dpsi__(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi * self.x), sum(shape_dxi * self.y), sum(shape_dxi * self.z)],
                [sum(shape_deta * self.x), sum(shape_deta * self.y), sum(shape_deta * self.z)],
                [sum(shape_dpsi * self.x), sum(shape_dpsi * self.y), sum(shape_dpsi * self.z)]
            ])
            # Якобиан
            jacobian = det(jacobi)
            # Вспомогательная матрица для построения матриц масс и демпфирования
            c = array([
                    [shape[0], 0.0, 0.0, shape[1], 0.0, 0.0, shape[2], 0.0, 0.0, shape[2], 0.0, 0.0, shape[3], 0.0, 0.0,
                     shape[4], 0.0, 0.0, shape[5], 0.0, 0.0, shape[6], 0.0, 0.0, shape[7], 0.0, 0.0],
                    [0.0, shape[0], 0.0, 0.0, shape[1], 0.0, 0.0, shape[2], 0.0, 0.0, shape[2], 0.0, 0.0, shape[3], 0.0,
                     0.0, shape[4], 0.0, 0.0, shape[5], 0.0, 0.0, shape[6], 0.0, 0.0, shape[7], 0.0],
                    [0.0, 0.0, shape[0], 0.0, 0.0, shape[1], 0.0, 0.0, shape[2], 0.0, 0.0, shape[2], 0.0, 0.0, shape[3],
                     0.0, 0.0, shape[4], 0.0, 0.0, shape[5], 0.0, 0.0, shape[6], 0.0, 0.0, shape[7]]
                ])
            local_m += c.conj().transpose().dot(c)*jacobian*self.__w__[i]
        m = local_m*self.density
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Четырехугольный КЭ плиты
class TFE2D4P(TFE2D4):
    def __init__(self):
        super().__init__()

    def calc(self, u):
        d = self.e[0]*self.h**3/(1 - self.m[0]**2)/12.0
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        res = zeros((6, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = self.a[j][1] + self.a[j][3] * self.y[i]
                dy[i][j] = self.a[j][2] + self.a[j][3] * self.x[i]
                res[0][i] += u[2*j + 1]*dx[i][j]
                res[1][i] += u[2*j + 2]*dy[i][j]
                res[2][i] += u[2*j + 1]*dy[i][j] + u[2*j + 2]*dx[i][j]
                res[3][i] += d*(u[2*j + 1]*dx[i][j] + self.m[0]*u[2*j + 2]*dy[i][j])
                res[4][i] += d*(u[2*j + 2]*dy[i][j] + self.m[0]*u[2*j + 1]*dx[i][j])
                res[5][i] += 2.0*(1.0 - self.m[0])*d*(u[2*j + 1]*dy[i][j] + u[2*j + 2]*dx[i][j])
        return res

    # Формирование локальной матрицы жесткости
    def _generate_stiffness_matrix_(self):
        local_k = zeros((3*self.size, 3*self.size))
        # Матрицы упругих свойст
        k = 5.0/6.0
        cb = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
            ])*(self.e[0]*self.h**3)/(1.0 - self.m[0]**2)/12.0
        cs = array([
            [1.0, 0.0],
            [0.0, 1.0],
            ])*self.e[0]*self.h*k/(2.0 + 2.0*self.m[0])
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
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
            # Матрицы градиентов
            b1 = array([
                [0, 0, -shape_dx[0], 0, 0, -shape_dx[1], 0, 0, -shape_dx[2], 0, 0, -shape_dx[3]],
                [0, shape_dy[0], 0, 0, shape_dy[1], 0, 0, shape_dy[2], 0, 0, shape_dy[3], 0],
                [0, shape_dx[0], -shape_dy[0], 0, shape_dx[1], -shape_dy[1], 0, shape_dx[2], -shape_dy[2], 0,
                 shape_dx[3], -shape_dy[3]]
                ])
            b2 = array([
                [shape_dx[0], 0, shape[0], shape_dx[1], 0, shape[1], shape_dx[2], 0, shape[2], shape_dx[3], 0,
                 shape[3]],
                [shape_dy[0], -shape[0], 0, shape_dy[1], -shape[1], 0, shape_dy[2], -shape[2], 0, shape_dy[3],
                 -shape[3], 0]
                ])
            local_k += (b1.conj().transpose().dot(cb).dot(b1) +
                        b2.conj().transpose().dot(cs).dot(b2))*jacobian*self.__w__[i]*0.5
        return local_k

    # Формирование локальных матриц жесткости, масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        local_m = zeros((3*self.size, 3*self.size))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x), sum(shape_dxi*self.y)],
                [sum(shape_deta*self.x), sum(shape_deta*self.y)]
                ])
            # Якобиан
            jacobian = det(jacobi)
            # Вспомогательная матрица для построения матриц масс и демпфирования
            c = array([
                    [shape[0], 0, 0, shape[1], 0, 0, shape[2], 0, 0, shape[3], 0, 0],
                    [0, shape[0], 0, 0, shape[1], 0, 0, shape[2], 0, 0, shape[3], 0],
                    [0, 0, shape[0], 0, 0, shape[1], 0, 0, shape[2], 0, 0, shape[3]],
                ])
            mi = array([
                    [self.density*self.h, 0, 0],
                    [0, self.density*self.h**3/12.0, 0],
                    [0, 0, self.density*self.h**3/12.0],
                ])
            local_m += c.conj().transpose().dot(mi).dot(c)*jacobian*self.__w__[i]
        m = local_m*self.density
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Треугольный КЭ пластины
class TFE2D3P(TFE2D3):
    def __init__(self):
        super().__init__()
        # Параметры квадратур Гаусса
        self.__xi__ = [0.57735027, 0.57735027]
        self.__eta__ = [0.57735027, 0.57735027]
        self.__w__ = [1.0, 1.0]

    def calc(self, u):
        d = self.e[0]*self.h**3/(1 - self.m[0]**2)/12.0
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        res = zeros((6, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = self.a[j][1]
                dy[i][j] = self.a[j][2]
                res[0][i] += u[2*j + 1]*dx[i][j]
                res[1][i] += u[2*j + 2]*dy[i][j]
                res[2][i] += u[2*j + 1]*dy[i][j] + u[2*j + 2]*dx[i][j]
                res[3][i] += d*(u[2*j + 1]*dx[i][j] + self.m[0]*u[2*j + 2]*dy[i][j])
                res[4][i] += d*(u[2*j + 2]*dy[i][j] + self.m[0]*u[2*j + 1]*dx[i][j])
                res[5][i] += 2.0*d*(1.0 - self.m[0])*(u[2*j + 1]*dy[i][j] + u[2*j + 2]*dx[i][j])
        return res

    def __shape__(self, i):
        return array([1.0 - self.__xi__[i], self.__xi__[i] - self.__eta__[i], self.__eta__[i]])

    def _generate_stiffness_matrix_(self):
        local_k = zeros((3*self.size, 3*self.size))
        k = 5.0/6.0
        # Матрицы упругих свойст
        cb = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
            ])*(self.e[0]*self.h**3)/(1.0 - self.m[0]**2)/12.0
        cs = array([
            [1.0, 0.0],
            [0.0, 1.0],
            ])*self.e[0]*self.h*k/(2.0 + 2.0*self.m[0])
        # Интегрирование по треугольнику [0,0]-[1,0]-[1,1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = array([-1.0, 1.0, 0.0])
            shape_deta = array([0.0, -1.0, 1.0])
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
            # Матрицы градиентов
            b1 = array([
                [0, 0, -shape_dx[0], 0, 0, -shape_dx[1], 0, 0, -shape_dx[2]],
                [0, shape_dy[0], 0, 0, shape_dy[1], 0, 0, shape_dy[2], 0],
                [0, shape_dx[0], -shape_dy[0], 0, shape_dx[1], -shape_dy[1], 0, shape_dx[2], -shape_dy[2]]
                ])
            b2 = array([
                [shape_dx[0], 0, shape[0], shape_dx[1], 0, shape[1], shape_dx[2], 0, shape[2]],
                [shape_dy[0], -shape[0], 0, shape_dy[1], -shape[1], 0, shape_dy[2], -shape[2], 0]
                ])
            # Вычисление компонент локальной матрицы жесткости
            local_k += (b1.conj().transpose().dot(cb).dot(b1) +
                        b2.conj().transpose().dot(cs).dot(b2))*jacobian*self.__w__[i]/4.0
        return local_k

    def _generate_mass_damping_matrix_(self):
        local_m = zeros((3*self.size, 3*self.size))
        # Интегрирование по треугольнику [0,0]-[1,0]-[1,1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = array([-1.0, 1.0, 0.0])
            shape_deta = array([0.0, -1.0, 1.0])
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x), sum(shape_dxi*self.y)],
                [sum(shape_deta*self.x), sum(shape_deta*self.y)]
                ])
            # Якобиан
            jacobian = det(jacobi)
            # Вспомогательная матрица для построения матриц масс и демпфирования
            c = array([
                    [shape[0], 0, 0, shape[1], 0, 0, shape[2], 0, 0],
                    [0, shape[0], 0, 0, shape[1], 0, 0, shape[2], 0],
                    [0, 0, shape[0], 0, 0, shape[1], 0, 0, shape[2]]
                ])
            mi = array([
                    [self.h, 0, 0],
                    [0, self.h**3/12.0, 0],
                    [0, 0, self.h**3/12.0],
                ])
            local_m += c.conj().transpose().dot(mi).dot(c)*jacobian*self.__w__[i]
        m = local_m*self.density
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Треугольный КЭ оболочки
class TFE2D3S(TFE2D3P, TFE2D3):
    def __init__(self):
        super().__init__()

    def calc(self, u):
        d = self.e[0]*self.h**3/(1 - self.m[0]**2)/12.0
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        res = zeros((6, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = self.a[j][1]
                dy[i][j] = self.a[j][2]
                res[0][i] += u[2*j + 1]*dx[i][j]
                res[1][i] += u[2*j + 2]*dy[i][j]
                res[2][i] += u[2*j + 1]*dy[i][j] + u[2*j + 2]*dx[i][j]
                res[3][i] += d*(u[2*j + 1]*dx[i][j] + self.m[0]*u[2*j + 2]*dy[i][j])
                res[4][i] += d*(u[2*j + 2]*dy[i][j] + self.m[0]*u[2*j + 1]*dx[i][j])
                res[5][i] += 2.0*d*(1.0 - self.m[0])*(u[2*j + 1]*dy[i][j] + u[2*j + 2]*dx[i][j])
        return res

    def generate(self, is_static=True):
        local_k, local_m = zeros((6*self.size, 6*self.size)), zeros((6*self.size, 6*self.size))
        k = 5.0/6.0
        # Матрицы упругих свойст
        cb = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
            ])*(self.e[0]*self.h**3)/(1.0 - self.m[0]**2)/12.0
        cs = array([
            [1.0, 0.0],
            [0.0, 1.0],
            ])*self.e[0]*self.h*k/(2.0 + 2.0*self.m[0])
        # Параметры квадратур Гаусса
        xi = [0.57735027, 0.57735027]
        eta = [0.57735027, 0.57735027]
        w = [1.0, 1.0]
        # Интегрирование по треугольнику [0,0]-[1,0]-[1,1] (по формуле Гаусса)
        for i in range(len(w)):
            # Изопараметрические функции формы и их производные
            shape = array([1.0 - xi[i], xi[i] - eta[i], eta[i]])
            shape_dxi = array([-1.0, 1.0, 0.0])
            shape_deta = array([0.0, -1.0, 1.0])
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
            # Матрицы градиентов
            b1 = array([
                [0, 0, -shape_dx[0], 0, 0, -shape_dx[1], 0, 0, -shape_dx[2]],
                [0, shape_dy[0], 0, 0, shape_dy[1], 0, 0, shape_dy[2], 0],
                [0, shape_dx[0], -shape_dy[0], 0, shape_dx[1], -shape_dy[1], 0, shape_dx[2], -shape_dy[2]]
                ])
            b2 = array([
                [shape_dx[0], 0, shape[0], shape_dx[1], 0, shape[1], shape_dx[2], 0, shape[2]],
                [shape_dy[0], -shape[0], 0, shape_dy[1], -shape[1], 0, shape_dy[2], -shape[2], 0]
                ])
            # Вычисление компонент локальной матрицы жесткости
            local_k += (b1.conj().transpose().dot(cb).dot(b1) + b2.conj().transpose().dot(cs).dot(b2))*jacobian*w[i]/4.0
            if not is_static:
                # Вспомогательная матрица для построения матриц масс и демпфирования
                c = array([
                    [shape[0], 0, 0, shape[1], 0, 0, shape[2], 0, 0],
                    [0, shape[0], 0, 0, shape[1], 0, 0, shape[2], 0],
                    [0, 0, shape[0], 0, 0, shape[1], 0, 0, shape[2]]
                ])
                mi = array([
                    [self.h, 0, 0],
                    [0, self.h**3/12.0, 0],
                    [0, 0, self.h**3/12.0],
                ])
                local_m += c.conj().transpose().dot(mi).dot(c)*jacobian*w[i]
        # Заполнение локальных матриц жесткости, масс и демпфирования
        self.K = local_k
        if not is_static:
            self.M, self.D = self.density*local_m, self.damping*local_m
