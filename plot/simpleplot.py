#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#             Простая визуализация результатов расчета
###################################################################

import os.path
import simplejson as json
import matplotlib.pyplot as plt
import numpy as np
from math import floor
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from core.fem_result import TResult
from core.fem_object import error


class TSimplePlot:
    def __init__(self):
        self.__name__ = ''          # Название задачи
        self.__type__ = ''          # Тип задачи
        self.__date_time__ = ''     # Время и дата рачета задачи
        self.__fe_type__ = ''       # Тип КЭ
        self.__x__ = []             # Координаты узлов
        self.__fe__ = []            # Связи КЭ
        self.__be__ = []            # ... ГЭ
        self.__results__ = []       # Результаты расчета

    def set_results(self, file_name):
        if len(file_name) < 5 or file_name[len(file_name) - 5:] != '.json':
            file_name += '.json'

        if not os.path.exists(file_name):
            error('Unable open file ' + file_name)
            return False

        with open(file_name, 'r') as file:
            data = json.loads(json.load(file))
        header = data['header']
        self.__name__ = header['name']
        self.__type__ = header['type']
        self.__date_time__ = header['date_time']
        mesh = data['mesh']
        self.__fe_type__ = mesh['fe_type']
        self.__x__ = mesh['vertex']
        self.__fe__ = mesh['fe']
        self.__be__ = mesh['be']
        results = data['results']
        for i in range(0, len(results)):
            res = TResult()
            res.name = results[i]['function']
            res.t = results[i]['t']
            res.results = results[i]['results']
            self.__results__.append(res)
        return True

    # Визуализация заданной функции
    def plot(self, fun_name, t=0):
        # Поиск индекса функции в списке результатов
        index = -1
        for i in range(0, len(self.__results__)):
            if self.__results__[i].name == fun_name and self.__results__[i].t == t:
                index = i
                break
        if index == -1:
            error('Error: \'%s\' is not a recognized function name or incorrect time' % fun_name)
            return
        # Визуализация результата
        if self.__fe_type__ == 'fe_1d_2':
            self.__plot_1d_linear__(index)
        elif self.__fe_type__ == 'fe_2d_3':
            self.__plot_2d_tri__(index)
        elif self.__fe_type__ == 'fe_2d_4':
            self.__plot_2d_quad__(index)
        elif self.__fe_type__ == 'fe_3d_4':
            self.__plot_3d_tet__(index)
        elif self.__fe_type__ == 'fe_3d_8':
            self.__plot_3d_hex__(index)
        # Задание заголовка
        if self.__type__ == 'dynamic':
            fun_name += ' (t = %5.2f)' % t

        plt.gcf().canvas.set_window_title('Result image')
        plt.title(fun_name)
        plt.show()

    # Визуализация заданной функции в случае одномерного линейного КЭ
    def __plot_1d_linear__(self, index):
        plt.plot(self.__x__, self.__results__[index].results, '-', linewidth=2)

    # Визуализация заданной функции в случае плоской треугольной сетки
    def __plot_2d_tri__(self, index):
        plt.figure()
        plt.gca().set_aspect('equal')
#        c_map = cm.get_cmap(name='terrain', lut=None)
        c_map = cm.get_cmap(name='nipy_spectral', lut=None)
        plt.triplot([self.__x__[i][0] for i in range(0, len(self.__x__))],
                    [self.__x__[i][1] for i in range(0, len(self.__x__))],
                    self.__fe__, lw=0.5, color='white')
        plt.tricontourf([self.__x__[i][0] for i in range(0, len(self.__x__))],
                        [self.__x__[i][1] for i in range(0, len(self.__x__))],
                        self.__fe__, self.__results__[index].results, cmap=c_map)
        plt.colorbar()

    # Визуализация заданной функции в случае плоской четырехугольной сетки
    def __plot_2d_quad__(self, index):
        plt.figure()
        plt.gca().set_aspect('equal')
        c_map = cm.get_cmap(name='nipy_spectral', lut=None)
        tri = np.array([np.array([T[0], T[1], T[2]]) for T in self.__fe__])
        plt.tricontourf([self.__x__[i][0] for i in range(0, len(self.__x__))],
                        [self.__x__[i][1] for i in range(0, len(self.__x__))],
                        tri, self.__results__[index].results, cmap=c_map)
        tri = np.array([np.array([T[0], T[2], T[3]]) for T in self.__fe__])
        plt.tricontourf([self.__x__[i][0] for i in range(0, len(self.__x__))],
                        [self.__x__[i][1] for i in range(0, len(self.__x__))],
                        tri, self.__results__[index].results, cmap=c_map)
        plt.colorbar()

    # Визуализация заданной функции в случае КЭ в форме тетраэдра
    def __plot_3d_tet__(self, index):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        triangle_vertices = \
            np.array([np.array([[self.__x__[T[0]][0], self.__x__[T[0]][1], self.__x__[T[0]]][2],
                                [self.__x__[T[1]][0], self.__x__[T[1]][0], self.__x__[T[1]]][2],
                                [self.__x__[T[2]][0], self.__x__[T[2]][1], self.__x__[T[2]][2]]])
                      for T in self.__be__])

        c_map = cm.ScalarMappable()
        c_map.set_array([self.__results__[index].min(), self.__results__[index].max()])
        face_colors = self.get_surface_color(self.__results__[index].results)
        coll = Poly3DCollection(triangle_vertices, facecolors=face_colors, edgecolors='none')
        ax.add_collection(coll)
        ax.set_xlim(min(self.__x__[i][0] for i in range(0, len(self.__x__))),
                    max(self.__x__[i][0] for i in range(0, len(self.__x__))))
        ax.set_ylim(min(self.__x__[i][1] for i in range(0, len(self.__x__))),
                    max(self.__x__[i][1] for i in range(0, len(self.__x__))))
        ax.set_zlim(min(self.__x__[i][2] for i in range(0, len(self.__x__))),
                    max(self.__x__[i][2] for i in range(0, len(self.__x__))))
        plt.colorbar(c_map)

    # Визуализация заданной функции в случае кубического КЭ
    def __plot_3d_hex__(self, index):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_aspect("auto")
        ax.set_autoscale_on(True)
        triangle_vertices1 = \
            np.array([np.array([[self.__x__[T[0]][0], self.__x__[T[0]][1], self.__x__[T[0]][2]],
                                [self.__x__[T[1]][0], self.__x__[T[1]][1], self.__x__[T[1]][2]],
                                [self.__x__[T[2]][0], self.__x__[T[2]][1], self.__x__[T[2]][2]],
                                [self.__x__[T[0]][0], self.__x__[T[0]][1], self.__x__[T[0]][2]],
                                [self.__x__[T[2]][0], self.__x__[T[2]][1], self.__x__[T[2]][2]],
                                [self.__x__[T[3]][0], self.__x__[T[3]][1], self.__x__[T[3]][2]]])
                      for T in self.__be__])

        c_map = cm.ScalarMappable()
        c_map.set_array([self.__results__[index].min(), self.__results__[index].max()])
        face_colors = self.get_surface_color(self.__results__[index].results)
        coll = Poly3DCollection(triangle_vertices1, facecolors=face_colors, edgecolors='none')
        ax.add_collection(coll)

        ax.set_xlim(min(self.__x__[i][0] for i in range(0, len(self.__x__))),
                    max(self.__x__[i][0] for i in range(0, len(self.__x__))))
        ax.set_ylim(min(self.__x__[i][1] for i in range(0, len(self.__x__))),
                    max(self.__x__[i][1] for i in range(0, len(self.__x__))))
        ax.set_zlim(min(self.__x__[i][2] for i in range(0, len(self.__x__))),
                    max(self.__x__[i][2] for i in range(0, len(self.__x__))))
        plt.colorbar(c_map)

    # Определене цвета поверхностной грани
    def get_surface_color(self, res):
        # Градации цвета (спектр)
        colors = [  # красный - желтый
                    [1.00, 0.00, 0.00], [1.00, 0.25, 0.00], [1.00, 0.50, 0.00], [1.00, 0.75, 0.00],
                    # желтый - зеленый
                    [1.00, 1.00, 0.00], [0.75, 1.00, 0.00], [0.50, 1.00, 0.00], [0.25, 1.00, 0.00],
                    # зеленый - фиолетовый
                    [0.00, 1.00, 0.00], [0.00, 1.00, 0.25], [0.00, 1.00, 0.50], [0.00, 1.00, 0.75],
                    # фиолетовый - синий
                    [0.00, 1.00, 1.00], [0.00, 0.75, 1.00], [0.00, 0.50, 1.00], [0.00, 0.00, 1.00]]
        u_min = min(res)
        u_max = max(res)
        face_colors = []
        for i in range(0, len(self.__be__)):
            u = 0
            for j in range(0, len(self.__be__[i])):
                u += res[self.__be__[i][j]]
            u /= len(self.__be__[i])
            index = floor((u - u_min)/((u_max - u_min)/16.0))
            if index > 15:
                index = 15
            face_colors.append(colors[index])
        return face_colors
