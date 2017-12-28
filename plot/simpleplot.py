#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#             Простая визуализация результатов расчета
###################################################################

import simplejson as json
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from core.fem_result import TResult


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

    def set_file(self, file_name):
        with open(file_name, 'r') as file:
            data = json.loads(json.load(file))
        header = data['header']
        self.__name__ = header['name']
        self.__type__ = header['type']
        self.__date_time__ = header['date_time']
        mesh = data['mesh']
        self.__x__ = mesh['vertex']
        self.__fe__ = mesh['fe']
        self.__be__ = mesh['be']
        results = data['results']
        for i in range(0, len(results)):
            res = TResult()
            res.name = results['function']
            res.t = results['t']
            res.results = results['results']
            self.__results__.append(res)

    # Визуализация заданной функции
    def plot(self, fun_name, t=0):
        # Проверка корректности задания времени
        if self.__type__ == 'dynamic' and \
                ((t < self.__params__.t0 or t > self.__params__.t1) or t % self.__params__.th > eps):
            error('Error: incorrectly specified the time: %5.2f' % t)
            return
        # Поиск индекса функции в списке результатов
        index = -1
        for i in range(0, len(self.__results__)):
            if self.__results__[i].name == fun_name and self.__results__[i].t == t:
                index = i
                break
        if index == -1:
            error('Error: \'%s\' is not a recognized function name' % fun_name)
            return
        # Визуализация результата
        if self.__mesh__.fe_type == 'fe_1d_2':
            self.__plot_1d_linear__(index)
        elif self.__mesh__.fe_type == 'fe_2d_3':
            self.__plot_2d_tri__(index)
        elif self.__mesh__.fe_type == 'fe_2d_4':
            self.__plot_2d_quad__(index)
        elif self.__mesh__.fe_type == 'fe_3d_4':
            self.__plot_3d_tet__(index)
        elif self.__mesh__.fe_type == 'fe_3d_8':
            self.__plot_3d_hex__(index)
        # Задание заголовка
        if self.__params__.problem_type == 'dynamic':
            fun_name += ' (t = %5.2f)' % t

        plt.gcf().canvas.set_window_title('Result image')
        plt.title(fun_name)
        plt.show()

    # Визуализация заданной функции в случае одномерного линейного КЭ
    def __plot_1d_linear__(self, index):
        plt.plot(self.__mesh__.x, self.__results__[index].results, '-', linewidth=2)

    # Визуализация заданной функции в случае плоской треугольной сетки
    def __plot_2d_tri__(self, index):
        plt.figure()
        plt.gca().set_aspect('equal')
#        c_map = cm.get_cmap(name='terrain', lut=None)
        c_map = cm.get_cmap(name='nipy_spectral', lut=None)
        plt.triplot([self.__mesh__.x[i][0] for i in range(0, len(self.__mesh__.x))],
                    [self.__mesh__.x[i][1] for i in range(0, len(self.__mesh__.x))],
                    self.__mesh__.fe, lw=0.5, color='white')
        plt.tricontourf([self.__mesh__.x[i][0] for i in range(0, len(self.__mesh__.x))],
                        [self.__mesh__.x[i][1] for i in range(0, len(self.__mesh__.x))],
                        self.__mesh__.fe, self.__results__[index].results, cmap=c_map)
        plt.colorbar()

    # Визуализация заданной функции в случае плоской четырехугольной сетки
    def __plot_2d_quad__(self, index):
        plt.figure()
        plt.gca().set_aspect('equal')
        c_map = cm.get_cmap(name='nipy_spectral', lut=None)
        tri = np.array([np.array([T[0], T[1], T[2]]) for T in self.__mesh__.fe])
        plt.tricontourf([self.__mesh__.x[i][0] for i in range(0, len(self.__mesh__.x))],
                        [self.__mesh__.x[i][1] for i in range(0, len(self.__mesh__.x))],
                        tri, self.__results__[index].results, cmap=c_map)
        tri = np.array([np.array([T[0], T[2], T[3]]) for T in self.__mesh__.fe])
        plt.tricontourf([self.__mesh__.x[i][0] for i in range(0, len(self.__mesh__.x))],
                        [self.__mesh__.x[i][1] for i in range(0, len(self.__mesh__.x))],
                        tri, self.__results__[index].results, cmap=c_map)
        plt.colorbar()

    # Визуализация заданной функции в случае КЭ в форме тетраэдра
    def __plot_3d_tet__(self, index):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        triangle_vertices = \
            np.array([np.array([[self.__mesh__.x[T[0]][0], self.__mesh__.x[T[0]][1], self.__mesh__.x[T[0]]][2],
                                [self.__mesh__.x[T[1]][0], self.__mesh__.x[T[1]][0], self.__mesh__.x[T[1]]][2],
                                [self.__mesh__.x[T[2]][0], self.__mesh__.x[T[2]][1], self.__mesh__.x[T[2]][2]]])
                      for T in self.__mesh__.be])

        c_map = cm.ScalarMappable()
        c_map.set_array([self.__results__[index].min(), self.__results__[index].max()])
        face_colors = self.get_surface_color(self.__results__[index].results)
        coll = Poly3DCollection(triangle_vertices, facecolors=face_colors, edgecolors='none')
        ax.add_collection(coll)
        ax.set_xlim(min(self.__mesh__.x[i][0] for i in range(0, len(self.__mesh__.x))),
                    max(self.__mesh__.x[i][0] for i in range(0, len(self.__mesh__.x))))
        ax.set_ylim(min(self.__mesh__.x[i][1] for i in range(0, len(self.__mesh__.x))),
                    max(self.__mesh__.x[i][1] for i in range(0, len(self.__mesh__.x))))
        ax.set_zlim(min(self.__mesh__.x[i][2] for i in range(0, len(self.__mesh__.x))),
                    max(self.__mesh__.x[i][2] for i in range(0, len(self.__mesh__.x))))
        plt.colorbar(c_map)

    # Визуализация заданной функции в случае кубического КЭ
    def __plot_3d_hex__(self, index):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_aspect("auto")
        ax.set_autoscale_on(True)

        # for i in range(0, len(self.__mesh__.surface)):
        #     for j in range(0, len(self.__mesh__.surface[i])):
        #         ind1 = j
        #         ind2 = j + 1 if j < len(self.__mesh__.surface[i]) - 1 else 0
        #         x = [self.__mesh__.x[self.__mesh__.surface[i][ind1]], self.__mesh__.x[self.__mesh__.surface[i][ind2]]]
        #         y = [self.__mesh__.y[self.__mesh__.surface[i][ind1]], self.__mesh__.y[self.__mesh__.surface[i][ind2]]]
        #         z = [self.__mesh__.z[self.__mesh__.surface[i][ind1]], self.__mesh__.z[self.__mesh__.surface[i][ind2]]]
        #         ax.plot3D(x, y, z, color="w")

        triangle_vertices1 = \
            np.array([np.array([[self.__mesh__.x[T[0]][0], self.__mesh__.x[T[0]][1], self.__mesh__.x[T[0]][2]],
                                [self.__mesh__.x[T[1]][0], self.__mesh__.x[T[1]][1], self.__mesh__.x[T[1]][2]],
                                [self.__mesh__.x[T[2]][0], self.__mesh__.x[T[2]][1], self.__mesh__.x[T[2]][2]],
                                [self.__mesh__.x[T[0]][0], self.__mesh__.x[T[0]][1], self.__mesh__.x[T[0]][2]],
                                [self.__mesh__.x[T[2]][0], self.__mesh__.x[T[2]][1], self.__mesh__.x[T[2]][2]],
                                [self.__mesh__.x[T[3]][0], self.__mesh__.x[T[3]][1], self.__mesh__.x[T[3]][2]]])
                      for T in self.__mesh__.be])

        c_map = cm.ScalarMappable()
        c_map.set_array([self.__results__[index].min(), self.__results__[index].max()])
        face_colors = self.get_surface_color(self.__results__[index].results)
        coll = Poly3DCollection(triangle_vertices1, facecolors=face_colors, edgecolors='none')
        ax.add_collection(coll)

        ax.set_xlim(min(self.__mesh__.x[i][0] for i in range(0, len(self.__mesh__.x))),
                    max(self.__mesh__.x[i][0] for i in range(0, len(self.__mesh__.x))))
        ax.set_ylim(min(self.__mesh__.x[i][1] for i in range(0, len(self.__mesh__.x))),
                    max(self.__mesh__.x[i][1] for i in range(0, len(self.__mesh__.x))))
        ax.set_zlim(min(self.__mesh__.x[i][2] for i in range(0, len(self.__mesh__.x))),
                    max(self.__mesh__.x[i][2] for i in range(0, len(self.__mesh__.x))))
        plt.colorbar(c_map)
