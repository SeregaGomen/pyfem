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
from numpy import identity
from numpy.linalg import det
from numpy.linalg import inv
from numpy.linalg import norm
from numpy import cross
from core.fem_error import TFEMException
from core.fem_defs import eps


# Абстрактный базовый класс, описывающий конечный элемент (КЭ)
class TFE:
    def __init__(self):
        self.size = 0           # Размерность КЭ
        self.freedom = 0        # Количество степеней свободы
        self.h = 1              # Толщина (для оболочек и пластин)
        self.density = 0        # Плотность
        self.damping = [0, 0]   # Параметры демпфирования
        self.e = []             # Модуль Юнга
        self.m = []             # Коэффициент Пуассона
        self.x = []             # Координаты вершин КЭ
        self.K = []             # Локальная матрица жесткости
        self.M = []             # ... масс
        self.C = []             # ... демпфирования
        self.a = []             # Коэффициенты функций форм

    # Задание толщины
    def set_h(self, h):
        self.h = h

    # Задание координат
    def set_coord(self, x):
        self.x = array(x)
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

    # Построение вектора для заданной стороны элемента
    def __vector__(self, i, j):
        v = array(self.x[j]) - array(self.x[i])
        return v/norm(v)

    # Матрица преобразования в локальную систему координат
    def __create_transform_matrix__(self):
        v_x = self.__vector__(1, 0)
        v_z = self.__cross_product__(self.__vector__(1, 0), self.__vector__(2, 0))
        v_y = self.__cross_product__(v_z, v_x)
        return array([v_x, v_y, v_z])

    # Векторное произведение a и b
    @staticmethod
    def __cross_product__(a, b):
        v = cross(a, b)
        return v/norm(v)

    # Вычисление функций форм КЭ
    @abstractmethod
    def __create__(self):
        raise NotImplementedError('Method TFE.__create__() is pure virtual')

    # Матрица упругих свойств
    @abstractmethod
    def __elastic_matrix__(self):
        raise NotImplementedError('Method TFE.__elastic_matrix__() is pure virtual')

    # Вычисления стандартных результатов КЭ
    @abstractmethod
    def calc(self, u):
        raise NotImplementedError('Method TFE.calc() is pure virtual')

    # Вычисление матрицы жесткости
    @abstractmethod
    def _generate_stiffness_matrix_(self):
        raise NotImplementedError('Method TFE._generate_stiffness_matrix_() is pure virtual')

    # Вычисление матриц масс и демпфирования
    @abstractmethod
    def _generate_mass_damping_matrix_(self):
        raise NotImplementedError('Method TFE._generate_mass_damping_matrix_() is pure virtual')


# Линейный (двухузловой) одномерный КЭ
class TFE1D2(TFE):
    def __init__(self):
        super().__init__()
        self.size = 2
        self.freedom = 1

    def __length__(self):
        return math.fabs(self.x[1][0] - self.x[0][0])

    def __create__(self):
        if self.__length__() == 0.0:
            raise TFEMException('incorrect_fe_err')
        self.a = zeros((self.size, self.size))
        self.a[0][0] = self.x[1][0]/(self.x[1][0] - self.x[0][0])
        self.a[0][1] = -1.0/(self.x[1][0] - self.x[0][0])
        self.a[1][0] = self.x[0][0]/(self.x[0][0] - self.x[1][0])
        self.a[1][1] = -1.0/(self.x[0][0] - self.x[1][0])

    def calc(self, u):
        res = zeros((self.size, self.size))
        res[0][0] = res[0][1] = u[0]*self.a[0][1] + u[1]*self.a[1][1]
        res[1][0] = res[1][1] = self.e[0]*(u[0]*self.a[0][1] + u[1]*self.a[1][1])
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
        self.freedom = 2

    def __elastic_matrix__(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
        ])*self.e[0]/(1.0 - self.m[0]**2)
        return d

    def __gradient_matrix__(self):
        # Матрица градиентов
        b = zeros([3, self.freedom*self.size])
        for i in range(0, self.size):
            # Производные функций формы
            dx = self.a[i][1]
            dy = self.a[i][2]
            b[0][self.freedom*i + 0] = dx
            b[1][self.freedom*i + 1] = dy
            b[2][self.freedom*i + 0] = dy
            b[2][self.freedom*i + 1] = dx
        return b

    def __square__(self):
        a = math.sqrt((self.x[0][0] - self.x[1][0])**2 + (self.x[0][1] - self.x[1][1])**2)
        b = math.sqrt((self.x[0][0] - self.x[2][0])**2 + (self.x[0][1] - self.x[2][1])**2)
        c = math.sqrt((self.x[2][0] - self.x[1][0])**2 + (self.x[2][1] - self.x[1][1])**2)
        p = 0.5*(a + b + c)
        return math.sqrt(p*(p - a)*(p - b)*(p - c))

    def __create__(self):
        det0 = self.x[2][1]*self.x[1][0] - self.x[2][1]*self.x[0][0] - self.x[0][1]*self.x[1][0] - \
               self.x[1][1]*self.x[2][0] + self.x[1][1]*self.x[0][0] + self.x[0][1]*self.x[2][0]
        if math.fabs(det0) < eps:
            raise TFEMException('incorrect_fe_err')
        index = [[2, 1], [0, 2], [1, 0]]
        self.a = zeros((self.size, self.size))
        for i in range(0, self.size):
            det1 = self.x[index[i][0]][1]*self.x[index[i][1]][0] - self.x[index[i][1]][1]*self.x[index[i][0]][0]
            det2 = self.x[index[i][1]][1] - self.x[index[i][0]][1]
            det3 = self.x[index[i][0]][0] - self.x[index[i][1]][0]
            self.a[i][0] = det1/det0
            self.a[i][1] = det2/det0
            self.a[i][2] = det3/det0

    def calc(self, u):
        e = self.__gradient_matrix__().dot(u)
        s = self.__elastic_matrix__().dot(e)
        res = zeros((6, self.size))
        for i in range(0, 3):
            for j in range(0, self.size):
                res[i][j] = e[i]
                res[i + 3][j] = s[i]
        return res

    # Формирование локальной матриц жесткости
    def _generate_stiffness_matrix_(self):
        return self.__gradient_matrix__().conj().transpose().dot(self.__elastic_matrix__()).\
                   dot(self.__gradient_matrix__())*self.__square__()*self.h

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
        self.freedom = 3

    def __elastic_matrix__(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), 1.0, self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0])],
        ])*self.e[0]*(1.0 - self.m[0])/(1.0 + self.m[0])/(1.0 - 2.0*self.m[0])
        return d

    def __gradient_matrix__(self):
        # Матрица градиентов
        b = zeros([6, self.size*self.freedom])
        for i in range(0, self.size):
            # Производные функций формы
            dx = self.a[i][1]
            dy = self.a[i][2]
            dz = self.a[i][3]

            b[0][self.freedom*i + 0] = dx
            b[1][self.freedom*i + 1] = dy
            b[2][self.freedom*i + 2] = dz
            b[3][self.freedom*i + 0] = dy
            b[3][self.freedom*i + 1] = dx
            b[4][self.freedom*i + 1] = dz
            b[4][self.freedom*i + 2] = dy
            b[5][self.freedom*i + 0] = dz
            b[5][self.freedom*i + 2] = dx
        return b

    def __volume__(self):
        a = (self.x[1][0] - self.x[0][0])*(self.x[2][1] - self.x[0][1])*(self.x[3][2] - self.x[0][2]) + \
            (self.x[3][0] - self.x[0][0])*(self.x[1][1] - self.x[0][1])*(self.x[2][2] - self.x[0][2]) + \
            (self.x[2][0] - self.x[0][0])*(self.x[3][1] - self.x[0][1])*(self.x[1][2] - self.x[0][2])
        b = (self.x[3][0] - self.x[0][0])*(self.x[2][1] - self.x[0][1])*(self.x[1][2] - self.x[0][2]) + \
            (self.x[2][0] - self.x[0][0])*(self.x[1][1] - self.x[0][1])*(self.x[3][2] - self.x[0][2]) + \
            (self.x[1][0] - self.x[0][0])*(self.x[3][1] - self.x[0][1])*(self.x[2][2] - self.x[0][2])
        return math.fabs(a - b)/6.0

    def __create__(self):
        if self.__volume__() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a, self.a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][2]
            b[j] = 1.0
            x = solve(a, b)
            self.a[j] = list(x)

    def calc(self, u):
        e = self.__gradient_matrix__().dot(u)
        s = self.__elastic_matrix__().dot(e)
        res = zeros((12, self.size))
        for i in range(0, 6):
            for j in range(0, self.size):
                res[i][j] = e[i]
                res[i + 6][j] = s[i]
        return res

    # Формирование локальной матрицы жесткости
    def _generate_stiffness_matrix_(self):
        return self.__gradient_matrix__().conj().transpose().dot(self.__elastic_matrix__()).\
                   dot(self.__gradient_matrix__())*self.__volume__()

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
        self.freedom = 2
        # Параметры квадратур Гаусса
        self.__xi__ = [-0.57735027, -0.57735027, 0.57735027, 0.57735027]
        self.__eta__ = [-0.57735027, 0.57735027, -0.57735027, 0.57735027]
        self.__w__ = [1.0, 1.0, 1.0, 1.0]

    def __square__(self):
        return math.sqrt((self.x[0][0] - self.x[1][0])**2 + (self.x[0][1] - self.x[1][1])**2)

    def __create__(self):
        if self.__square__() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a, self.a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][0]*self.x[i][1]
            b[j] = 1.0
            x = solve(a, b)
            self.a[j] = list(x)

    def __elastic_matrix__(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
        ])*self.e[0]/(1.0 - self.m[0]**2)
        return d

    def calc(self, u):
        res = zeros((6, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([3, self.size*self.freedom])
            for j in range(0, self.size):
                dx = self.a[j][1] + self.a[j][3] * self.x[i][1]
                dy = self.a[j][2] + self.a[j][3] * self.x[i][0]
                b[0][j*self.freedom + 0] = dx
                b[1][j*self.freedom + 1] = dy
                b[2][j*self.freedom + 0] = dy
                b[2][j*self.freedom + 1] = dx
            e = b.dot(u)
            s = self.__elastic_matrix__().dot(e)
            for j in range(0, 3):
                res[j][i] += e[j]
                res[j + 3][i] += s[j]
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
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta
            shape_dy = inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta
            # Изопараметрическвая матрица градиентов
            b = zeros([3, self.freedom*self.size])
            for j in range(0, self.size):
                b[0][j*self.freedom + 0] = shape_dx[j]
                b[1][j*self.freedom + 1] = shape_dy[j]
                b[2][j*self.freedom + 0] = shape_dy[j]
                b[2][j*self.freedom + 1] = shape_dx[j]
            local_k += b.conj().transpose().dot(self.__elastic_matrix__()).dot(b)*jacobian*self.__w__[i]*self.h
        return local_k

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        local_m = zeros((self.freedom*self.size, self.freedom*self.size))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
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
        self.freedom = 3
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
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][2]
                a[i][4] = self.x[i][0]*self.x[i][1]
                a[i][5] = self.x[i][0]*self.x[i][2]
                a[i][6] = self.x[i][1]*self.x[i][2]
                a[i][7] = self.x[i][0]*self.x[i][1]*self.x[i][2]
            b[j] = 1.0
            try:
                x = solve(a, b)
            except LinAlgError:
                raise TFEMException('incorrect_fe_err')
            self.a[j] = list(x)

    def __elastic_matrix__(self):
        # Матрица упругих свойст
        return array([
            [1.0, self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), 1.0, self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0])],
        ])*self.e[0]*(1.0 - self.m[0])/(1.0 + self.m[0])/(1.0 - 2.0*self.m[0])

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([6, self.size*self.freedom])
            for j in range(0, self.size):
                dx = self.a[j][1] + self.a[j][4] * self.x[i][1] + self.a[j][5] * self.x[i][2] + \
                     self.a[j][7] * self.x[i][1] * self.x[i][2]
                dy = self.a[j][2] + self.a[j][4] * self.x[i][0] + self.a[j][6] * self.x[i][2] + \
                     self.a[j][7] * self.x[i][0] * self.x[i][2]
                dz = self.a[j][3] + self.a[j][5] * self.x[i][0] + self.a[j][6] * self.x[i][1] + \
                     self.a[j][7] * self.x[i][0] * self.x[i][1]
                b[0][j*self.freedom + 0] = dx
                b[1][j*self.freedom + 1] = dy
                b[2][j*self.freedom + 2] = dz
                b[3][j*self.freedom + 0] = dy
                b[3][j*self.freedom + 1] = dx
                b[4][j*self.freedom + 1] = dz
                b[4][j*self.freedom + 2] = dy
                b[5][j*self.freedom + 0] = dz
                b[5][j*self.freedom + 2] = dx
            e = b.dot(u)
            s = self.__elastic_matrix__().dot(e)
            for j in range(0, 6):
                res[j][i] += e[j]
                res[j + 6][i] += s[j]
        return res

    # Изопараметрические функции формы и их производные
    def __shape__(self, i):
        return array([
            0.125*(1.0 - self.__xi__[i])*(1.0 - self.__eta__[i])*(1.0 - self.__psi__[i]),
            0.125*(1.0 + self.__xi__[i])*(1.0 - self.__eta__[i])*(1.0 - self.__psi__[i]),
            0.125*(1.0 + self.__xi__[i])*(1.0 + self.__eta__[i])*(1.0 - self.__psi__[i]),
            0.125*(1.0 - self.__xi__[i])*(1.0 + self.__eta__[i])*(1.0 - self.__psi__[i]),
            0.125*(1.0 - self.__xi__[i])*(1.0 - self.__eta__[i])*(1.0 + self.__psi__[i]),
            0.125*(1.0 + self.__xi__[i])*(1.0 - self.__eta__[i])*(1.0 + self.__psi__[i]),
            0.125*(1.0 + self.__xi__[i])*(1.0 + self.__eta__[i])*(1.0 + self.__psi__[i]),
            0.125*(1.0 - self.__xi__[i])*(1.0 + self.__eta__[i])*(1.0 + self.__psi__[i])
        ])

    def __shape_dxi__(self, i):
        return array([
            -0.125*(1.0 - self.__eta__[i])*(1.0 - self.__psi__[i]),
            0.125*(1.0 - self.__eta__[i])*(1.0 - self.__psi__[i]),
            0.125*(1.0 + self.__eta__[i])*(1.0 - self.__psi__[i]),
            -0.125*(1.0 + self.__eta__[i])*(1.0 - self.__psi__[i]),
            -0.125*(1.0 - self.__eta__[i])*(1.0 + self.__psi__[i]),
            0.125*(1.0 - self.__eta__[i])*(1.0 + self.__psi__[i]),
            0.125*(1.0 + self.__eta__[i])*(1.0 + self.__psi__[i]),
            -0.125*(1.0 + self.__eta__[i])*(1.0 + self.__psi__[i])
        ])

    def __shape_deta__(self, i):
        return array([
            -0.125*(1.0 - self.__xi__[i])*(1.0 - self.__psi__[i]),
            -0.125*(1.0 + self.__xi__[i])*(1.0 - self.__psi__[i]),
            0.125*(1.0 + self.__xi__[i])*(1.0 - self.__psi__[i]),
            0.125*(1.0 - self.__xi__[i])*(1.0 - self.__psi__[i]),
            -0.125*(1.0 - self.__xi__[i])*(1.0 + self.__psi__[i]),
            -0.125*(1.0 + self.__xi__[i])*(1.0 + self.__psi__[i]),
            0.125*(1.0 + self.__xi__[i])*(1.0 + self.__psi__[i]),
            0.125*(1.0 - self.__xi__[i])*(1.0 + self.__psi__[i])
        ])

    def __shape_dpsi__(self, i):
        return array([
            -0.125*(1.0 - self.__xi__[i])*(1.0 - self.__eta__[i]),
            -0.125*(1.0 + self.__xi__[i])*(1.0 - self.__eta__[i]),
            -0.125*(1.0 + self.__xi__[i])*(1.0 + self.__eta__[i]),
            -0.125*(1.0 - self.__xi__[i])*(1.0 + self.__eta__[i]),
            0.125*(1.0 - self.__xi__[i])*(1.0 - self.__eta__[i]),
            0.125*(1.0 + self.__xi__[i])*(1.0 - self.__eta__[i]),
            0.125*(1.0 + self.__xi__[i])*(1.0 + self.__eta__[i]),
            0.125*(1.0 - self.__xi__[i])*(1.0 + self.__eta__[i])
        ])

    # Формирование локальной ыматрицы жесткости
    def _generate_stiffness_matrix_(self):
        local_k = zeros((self.freedom*self.size, self.freedom*self.size))
        # Интегрирование по кубу [-1; 1] x [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
            shape_dpsi = self.__shape_dpsi__(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1]), sum(shape_dxi*self.x[:, 2])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1]), sum(shape_deta*self.x[:, 2])],
                [sum(shape_dpsi*self.x[:, 0]), sum(shape_dpsi*self.x[:, 1]), sum(shape_dpsi*self.x[:, 2])]
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
            # Изопараметрическая матрица градиентов
            b = zeros([6, self.freedom*self.size])
            for j in range(0, self.size):
                b[0][j*self.freedom + 0] = shape_dx[j]
                b[1][j*self.freedom + 1] = shape_dy[j]
                b[2][j*self.freedom + 2] = shape_dz[j]
                b[3][j*self.freedom + 0] = shape_dy[j]
                b[3][j*self.freedom + 1] = shape_dx[j]
                b[4][j*self.freedom + 1] = shape_dz[j]
                b[4][j*self.freedom + 2] = shape_dy[j]
                b[5][j*self.freedom + 0] = shape_dz[j]
                b[5][j*self.freedom + 2] = shape_dx[j]
            local_k += b.conj().transpose().dot(self.__elastic_matrix__()).dot(b)*jacobian*self.__w__[i]
        return local_k

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        local_m = zeros((self.freedom*self.size, self.freedom*self.size))
        # Интегрирование по кубу [-1; 1] x [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
            shape_dpsi = self.__shape_dpsi__(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1]), sum(shape_dxi*self.x[:, 2])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1]), sum(shape_deta*self.x[:, 2])],
                [sum(shape_dpsi*self.x[:, 0]), sum(shape_dpsi*self.x[:, 1]), sum(shape_dpsi*self.x[:, 2])]
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
        self.freedom = 3

    def __elastic_matrix__(self):
        d1 = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
        ])*(self.e[0])/(1.0 - self.m[0]**2)
        d2 = array([
            [1.0, 0.0],
            [0.0, 1.0],
        ])*self.e[0]/(2.0 + 2.0*self.m[0])
        return d1, d2

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b1 = zeros([3, self.size*self.freedom])
            b2 = zeros([2, self.size*self.freedom])
            for j in range(0, self.size):
                shape = self.a[j][0] + self.a[j][1]*self.x[i][0] + self.a[j][2]*self.x[i][1] + \
                        self.a[j][3]*self.x[i][0]*self.x[i][1]
                dx = self.a[j][1] + self.a[j][3]*self.x[i][1]
                dy = self.a[j][2] + self.a[j][3]*self.x[i][0]
                b1[0][self.freedom*j + 2] = -dx
                b1[1][self.freedom*j + 1] = dy
                b1[2][self.freedom*j + 1] = dx
                b1[2][self.freedom*j + 2] = -dy

                b2[0][self.freedom*j + 0] = dx
                b2[0][self.freedom*j + 2] = shape
                b2[1][self.freedom*j + 0] = dy
                b2[1][self.freedom*j + 1] = -shape

            d1, d2 = self.__elastic_matrix__()
            e1 = b1.dot(u)
            s1 = d1.dot(e1)
            e2 = b2.dot(u)
            s2 = d2.dot(e2)
            res[0][i] += e1[0]  # Exx
            res[1][i] += e1[1]  # Eyy
            res[3][i] += e1[2]  # Exy
            res[4][i] += e2[0]  # Exz
            res[5][i] += e2[1]  # Eyz

            res[6][i] += s1[0]  # Sxx
            res[7][i] += s1[1]  # Syy
            res[9][i] += s1[2]  # Sxy
            res[10][i] += s2[0]  # Sxz
            res[11][i] += s2[1]  # Syz
        return res

    # Формирование локальной матрицы жесткости
    def _generate_stiffness_matrix_(self):
        local_k = zeros((self.freedom*self.size, self.freedom*self.size))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = self.__shape_dxi__(i)
            shape_deta = self.__shape_deta__(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta
            shape_dy = inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta
            # Изопараметрические матрицы градиентов
            b1 = zeros([3, self.size*self.freedom])
            for j in range(0, self.size):
                b1[0][self.freedom*j + 2] = -shape_dx[j]
                b1[1][self.freedom*j + 1] = shape_dy[j]
                b1[2][self.freedom*j + 1] = shape_dx[j]
                b1[2][self.freedom*j + 2] = -shape_dy[j]
            b2 = zeros([2, self.size*self.freedom])
            for j in range(0, self.size):
                b2[0][self.freedom*j + 0] = shape_dx[j]
                b2[0][self.freedom*j + 2] = shape[j]
                b2[1][self.freedom*j + 0] = shape_dy[j]
                b2[1][self.freedom*j + 1] = -shape[j]
            d1, d2 = self.__elastic_matrix__()
            local_k += (b1.conj().transpose().dot(d1).dot(b1)*self.h**3/12.0 +
                        b2.conj().transpose().dot(d2).dot(b2)*self.h*5.0/6.0)*jacobian*self.__w__[i]
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
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
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
        self.freedom = 3
        # Параметры квадратур Гаусса
        self.__xi__ = [0.57735027, 0.57735027]
        self.__eta__ = [0.57735027, 0.57735027]
        self.__w__ = [1.0, 1.0]

    def __elastic_matrix__(self):
        # Матрицы упругих свойст
        d1 = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
        ])*(self.e[0])/(1.0 - self.m[0]**2)
        d2 = array([
            [1.0, 0.0],
            [0.0, 1.0],
        ])*self.e[0]/(2.0 + 2.0*self.m[0])
        return d1, d2

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b1 = zeros([3, self.size*self.freedom])
            b2 = zeros([2, self.size*self.freedom])
            for j in range(0, self.size):
                shape = self.a[j][0] + self.a[j][1]*self.x[i][0] + self.a[j][2]*self.x[i][1]
                dx = self.a[j][1]
                dy = self.a[j][2]
                b1[0][self.freedom*j + 2] = -dx
                b1[1][self.freedom*j + 1] = dy
                b1[2][self.freedom*j + 1] = dx
                b1[2][self.freedom*j + 2] = -dy

                b2[0][self.freedom*j + 0] = dx
                b2[0][self.freedom*j + 2] = shape
                b2[1][self.freedom*j + 0] = dy
                b2[1][self.freedom*j + 1] = -shape
            d1, d2 = self.__elastic_matrix__()
            e1 = b1.dot(u)
            s1 = d1.dot(e1)
            e2 = b2.dot(u)
            s2 = d2.dot(e2)
            res[0][i] += e1[0]  # Exx
            res[1][i] += e1[1]  # Eyy
            res[3][i] += e1[2]  # Exy
            res[4][i] += e2[0]  # Exz
            res[5][i] += e2[1]  # Eyz

            res[6][i] += s1[0]  # Sxx
            res[7][i] += s1[1]  # Syy
            res[9][i] += s1[2]  # Sxy
            res[10][i] += s2[0]  # Sxz
            res[11][i] += s2[1]  # Syz
        return res

    def __shape__(self, i):
        return array([1.0 - self.__xi__[i], self.__xi__[i] - self.__eta__[i], self.__eta__[i]])

    def _generate_stiffness_matrix_(self):
        local_k = zeros((self.freedom*self.size, self.freedom*self.size))
        # Интегрирование по треугольнику [0,0]-[1,0]-[1,1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = array([-1.0, 1.0, 0.0])
            shape_deta = array([0.0, -1.0, 1.0])
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta
            shape_dy = inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta
            # Матрицы градиентов
            b1 = zeros([3, self.size*self.freedom])
            for j in range(0, self.size):
                b1[0][self.freedom*j + 2] = -shape_dx[j]
                b1[1][self.freedom*j + 1] = shape_dy[j]
                b1[2][self.freedom*j + 1] = shape_dx[j]
                b1[2][self.freedom*j + 2] = -shape_dy[j]
            b2 = zeros([2, self.size*self.freedom])
            for j in range(0, self.size):
                b2[0][self.freedom*j + 0] = shape_dx[j]
                b2[0][self.freedom*j + 2] = shape[j]
                b2[1][self.freedom*j + 0] = shape_dy[j]
                b2[1][self.freedom*j + 1] = -shape[j]
            # Вычисление компонент локальной матрицы жесткости
            d1, d2 = self.__elastic_matrix__()
            local_k += (b1.conj().transpose().dot(d1).dot(b1)*self.h**3/12.0 +
                        b2.conj().transpose().dot(d2).dot(b2)*self.h*5.0/6.0)*jacobian*self.__w__[i]/4.0
        return local_k

    def _generate_mass_damping_matrix_(self):
        local_m = zeros((self.freedom*self.size, self.freedom*self.size))
        # Интегрирование по треугольнику [0,0]-[1,0]-[1,1] (по формуле Гаусса)
        for i in range(len(self.__w__)):
            # Изопараметрические функции формы и их производные
            shape = self.__shape__(i)
            shape_dxi = array([-1.0, 1.0, 0.0])
            shape_deta = array([0.0, -1.0, 1.0])
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
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
class TFE2D3S(TFE2D3P):
    def __init__(self):
        super().__init__()
        self.size = 3
        self.freedom = 6
        self.T = zeros((3, 3))  # Матрица преобразования глобальных координат в локальные
        self.global_x = zeros((3, 3))

    def __create__(self):
        self.T = self.__create_transform_matrix__()
        self.global_x = self.x
        self.x = array([self.T.dot(self.x[0, :]), self.T.dot(self.x[1, :]), self.T.dot(self.x[2, :])])
        self.x -= self.x[0, :]
        TFE2D3.__create__(self)
#        self.x = x

    def __prepare_transform_matrix__(self):
        m = zeros((self.freedom*self.size, self.freedom*self.size))
        for i in range(0, 3):
            for j in range(0, 3):
                m[i][j] = self.T[i][j]
                m[i + 3][j + 3] = self.T[i][j]
                m[i + 6][j + 6] = self.T[i][j]
                m[i + 9][j + 9] = self.T[i][j]
                m[i + 12][j + 12] = self.T[i][j]
                m[i + 15][j + 15] = self.T[i][j]
        return m

    def calc(self, u):
        # Подготовка матрицы преобразования
        m = self.__prepare_transform_matrix__()
        # Модифицируем перемещения
        lu = m.dot(u)
        # Вычисляем деформации и напряжения "мембранной" составляющей
        membrane_u = []
        for i in range(0, self.size):
            membrane_u.append(lu[self.freedom*i])      # u
            membrane_u.append(lu[self.freedom*i + 1])  # w
        membrane_res = TFE2D3.calc(self, membrane_u)
        # Вычисляем деформации и напряжения пластины
        plate_u = []
        for i in range(0, self.size):
            plate_u.append(lu[self.freedom*i + 2])     # w
            plate_u.append(lu[self.freedom*i + 3])     # Tx
            plate_u.append(lu[self.freedom*i + 4])     # Ty
        plate_res = TFE2D3P.calc(self, plate_u)
        res = zeros((12, self.size))
        for i in range(0, self.size):
            ld = array([
                [plate_res[0][i] + membrane_res[0][i], plate_res[3][i] + membrane_res[2][i], plate_res[4][i]],
                [plate_res[3][i] + membrane_res[2][i], plate_res[1][i] + membrane_res[1][i], plate_res[5][i]],
                [plate_res[4][i], plate_res[5][i], 0]
            ])
            ls = array([
                [plate_res[6][i] + membrane_res[3][i], plate_res[9][i] + membrane_res[5][i], plate_res[10][i]],
                [plate_res[9][i] + membrane_res[5][i], plate_res[7][i] + membrane_res[4][i], plate_res[11][i]],
                [plate_res[10][i], plate_res[11][i], 0]
            ])
            gd = self.T.conj().transpose().dot(ld).dot(self.T)
            gs = self.T.conj().transpose().dot(ls).dot(self.T)

            res[0][i] = gd[0][0]    # Exx
            res[1][i] = gd[1][1]    # Eyy
            res[2][i] = gd[2][2]    # Ezz
            res[3][i] = gd[0][1]    # Exy
            res[4][i] = gd[0][2]    # Exz
            res[5][i] = gd[1][2]    # Eyz
            res[6][i] = gs[0][0]    # Sxx
            res[7][i] = gs[1][1]    # Syy
            res[8][i] = gs[2][2]    # Szz
            res[9][i] = gs[0][1]    # Sxy
            res[10][i] = gs[0][2]   # Sxz
            res[11][i] = gs[1][2]   # Syz
        return res

    def _generate_stiffness_matrix_(self):
        local_k = identity(self.freedom*self.size)
        # Создание матрицы преобразования
        m = self.__prepare_transform_matrix__()
        # Локальная матрицыа жесткости плоского КЭ
        k1 = TFE2D3._generate_stiffness_matrix_(self)
        # ... КЭ пластины
        k2 = TFE2D3P._generate_stiffness_matrix_(self)
        local_freedom1 = 2
        for i in range(0, len(k1)):
            p = (i//local_freedom1)*self.freedom + i % local_freedom1
            for j in range(0, len(k1)):
                q = (j//local_freedom1)*self.freedom + j % local_freedom1
                local_k[p][q] = k1[i][j]
        local_freedom2 = 3
        for i in range(0, len(k2)):
            p = (i//local_freedom2)*self.freedom + i % local_freedom2 + local_freedom1
            for j in range(0, len(k2)):
                q = (j//local_freedom2)*self.freedom + j % local_freedom2 + local_freedom1
                local_k[p][q] = k2[i][j]
        return m.conj().transpose().dot(local_k).dot(m)


# Четырехугольный КЭ оболочки
class TFE2D4S(TFE2D4P, TFE2D4):
    def __init__(self):
        super().__init__()
        self.size = 4
        self.freedom = 6
        self.T = zeros((3, 3))  # Матрица преобразования глобальных координат в локальные
        self.global_x = zeros((4, 4))

    def calc(self, u):
        # Подготовка матрицы преобразования
        m = self.__prepare_transform_matrix__()
        # Модифицируем перемещения
        lu = m.dot(u)
        # Вычисляем деформации и напряжения "мембранной" составляющей
        membrane_u = []
        for i in range(0, self.size):
            membrane_u.append(lu[self.freedom*i])      # u
            membrane_u.append(lu[self.freedom*i + 1])  # w
        membrane_res = TFE2D4.calc(self, membrane_u)
        # Вычисляем деформации и напряжения пластины
        plate_u = []
        for i in range(0, self.size):
            plate_u.append(lu[self.freedom*i + 2])     # w
            plate_u.append(lu[self.freedom*i + 3])     # Tx
            plate_u.append(lu[self.freedom*i + 4])     # Ty
        plate_res = TFE2D4P.calc(self, plate_u)
        res = zeros((12, self.size))
        for i in range(0, self.size):
            ld = array([
                [plate_res[0][i] + membrane_res[0][i], plate_res[3][i] + membrane_res[2][i], plate_res[4][i]],
                [plate_res[3][i] + membrane_res[2][i], plate_res[1][i] + membrane_res[1][i], plate_res[5][i]],
                [plate_res[4][i], plate_res[5][i], 0]
            ])
            ls = array([
                [plate_res[6][i] + membrane_res[3][i], plate_res[9][i] + membrane_res[5][i], plate_res[10][i]],
                [plate_res[9][i] + membrane_res[5][i], plate_res[7][i] + membrane_res[4][i], plate_res[11][i]],
                [plate_res[10][i], plate_res[11][i], 0]
            ])
            gd = self.T.conj().transpose().dot(ld).dot(self.T)
            gs = self.T.conj().transpose().dot(ls).dot(self.T)

            res[0][i] = gd[0][0]    # Exx
            res[1][i] = gd[1][1]    # Eyy
            res[2][i] = gd[2][2]    # Ezz
            res[3][i] = gd[0][1]    # Exy
            res[4][i] = gd[0][2]    # Exz
            res[5][i] = gd[1][2]    # Eyz
            res[6][i] = gs[0][0]    # Sxx
            res[7][i] = gs[1][1]    # Syy
            res[8][i] = gs[2][2]    # Szz
            res[9][i] = gs[0][1]    # Sxy
            res[10][i] = gs[0][2]   # Sxz
            res[11][i] = gs[1][2]   # Syz
        return res

    def __prepare_transform_matrix__(self):
        m = zeros((self.freedom*self.size, self.freedom*self.size))
        for i in range(0, 3):
            for j in range(0, 3):
                m[i][j] = self.T[i][j]
                m[i + 3][j + 3] = self.T[i][j]
                m[i + 6][j + 6] = self.T[i][j]
                m[i + 9][j + 9] = self.T[i][j]
                m[i + 12][j + 12] = self.T[i][j]
                m[i + 15][j + 15] = self.T[i][j]
                m[i + 18][j + 18] = self.T[i][j]
                m[i + 21][j + 21] = self.T[i][j]
        return m

    # Формирование локальной матрицы жесткости
    def _generate_stiffness_matrix_(self):
        local_k = identity(self.freedom*self.size)
        # Подготовка матрицы преобразования
        m = self.__prepare_transform_matrix__()
        # Локальная матрицыа жесткости плоского КЭ
#        k1 = TFE2D4._generate_stiffness_matrix_(self)

        k1 = super(TFE2D4S, self)._generate_stiffness_matrix_()


        # ... КЭ пластины
        k2 = TFE2D4P._generate_stiffness_matrix_(self)

        local_freedom1 = 2
        for i in range(0, len(k1)):
            p = (i//local_freedom1)*self.freedom + i % local_freedom1
            for j in range(0, len(k1)):
                q = (j//local_freedom1)*self.freedom + j % local_freedom1
                local_k[p][q] = k1[i][j]
        local_freedom2 = 3
        for i in range(0, len(k2)):
            p = (i//local_freedom2)*self.freedom + i % local_freedom2 + local_freedom1
            for j in range(0, len(k2)):
                q = (j//local_freedom2)*self.freedom + j % local_freedom2 + local_freedom1
                local_k[p][q] = k2[i][j]

        # import sys
        # print('\n******************************************')
        # for i in range(0, len(local_k)):
        #     for j in range(0, len(local_k)):
        #         sys.stdout.write('%+E\t' % local_k[i][j])
        #     sys.stdout.write('\n')
        # print('******************************************')
        return m.conj().transpose().dot(local_k).dot(m)

    def __create__(self):
        self.T = self.__create_transform_matrix__()
        self.global_x = self.x

        x0 = self.T.dot(self.x[0, :].conj().transpose())
        x1 = self.T.dot(self.x[1, :].conj().transpose())
        x2 = self.T.dot(self.x[2, :].conj().transpose())
        x3 = self.T.dot(self.x[3, :].conj().transpose())

        self.x = array([x0 - x0, x1 - x0, x2 - x0, x3 - x0])
#        self.x = array([x0, x1, x2, x3])
        TFE2D4.__create__(self)
#        self.x = x
