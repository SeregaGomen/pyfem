#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#     Визуализация результатов расчета с использованием OpenGL
###################################################################

import os.path
import simplejson as json
from math import floor
from core.fem_result import TResult
from core.fem_object import print_error
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout
from PyQt5.QtGui import QFont, QFontMetrics
from PyQt5.QtOpenGL import QGLWidget
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *


# Базовый класс, реализующий основной функционал OpenGL
class TOpenGLWidget(QWidget):
    def __init__(self, x, fe, be, results, index):
        super(TOpenGLWidget, self).__init__()
        self.x = x
        self.fe = fe
        self.be = be
        self.results = results[index].results
        self.min_x, self.max_x, self.x_c, self.radius = self.__get_coord_info__()
        self.min_u = min(self.results)
        self.max_u = max(self.results)
        self.is_light = False
        self.is_legend = True
        self.is_fe_border = False
        self.is_object_border = False
        self.diffuse = 0.8
        self.ambient = 0.8
        self.specular = 0.6
        self.num_color = 16
        self.color_table = []
        self.gl = QGLWidget(self)
        self.gl.initializeGL()
        self.gl.resizeGL = self.__resize__
        vb_layout = QVBoxLayout(self)
        vb_layout.addWidget(self.gl)
        self.__init_color_table__()

    def __get_coord_info__(self):
        min_x = [0, 0, 0]
        max_x = [0, 0, 0]
        for k in range(0, len(self.x[0])):
            min_x[k] = min(self.x[i][k] for i in range(0, len(self.x)))
            max_x[k] = max(self.x[i][k] for i in range(0, len(self.x)))
        x_c = [(max_x[0] + min_x[0])/2, (max_x[1] + min_x[1])/2, (max_x[2] + min_x[2])/2]
        radius = ((max_x[0] - min_x[0])**2 + (max_x[1] - min_x[1])**2 + (max_x[2] - min_x[2])**2)**0.5
        return min_x, max_x, x_c, radius

    def __init_color_table__(self):
        step = self.num_color/6
        h = 1.0/step
        if self.min_u == self.max_u:
            self.max_u += 1
        green = 0
        blue = red = 1
        u = self.min_u
        h_u = (self.max_u - self.min_u)/float(self.num_color)
        for i in range(0, self.num_color):
            if i < step:
                # фиолетовый-синий
                self.color_table.append([red, 0, 1, u])
                red -= h
                if red < 0:
                    red = 0
            elif step <= i < 2*step:
                # синий-голубой
                self.color_table.append([0, green, 1, u])
                green += h
                if green > 1:
                    green = 1
            elif 2*step <= i < 3*step:
                # голубой-зеленый
                self.color_table.append([0, 1, blue, u])
                blue -= h
                if blue < 0:
                    blue = 0
            elif 3*step <= i < 4*step:
                # зеленый-желтый
                self.color_table.append([red, 1, 0, u])
                red += h
                if red > 1:
                    red = 1
            elif i > 4*step:
                # желтый-оранжевый-красный
                self.color_table.append([1, green, 0, u])
                green -= 0.5*h
                if green < 0:
                    green = 0
            u += h_u

    def set_color(self, r, g, b, a):
        if self.is_light:
            self.__make_material__(r, g, b, a)
        else:
            glColor4f(r, g, b, a)

    def __resize__(self, w, h):
        aspect = w/h
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(60.0, aspect, 0.01*self.radius, 10*self.radius)
        gluLookAt(0, 0, self.radius, 0, 0, 0, 0, 1, 0)
        glViewport(0, 0, w, h)
        glMatrixMode(GL_MODELVIEW)

    def __make_material__(self, r, g, b, a):
        glLightfv(GL_LIGHT0, GL_AMBIENT, [0.4, 0.4, 0.4, 1.0])
        glLightfv(GL_LIGHT0, GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
        glLightfv(GL_LIGHT0, GL_SPECULAR, [0.7, 0.7, 0.7, 1.0])
        glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, [r*self.diffuse, g*self.diffuse, b*self.diffuse, a])
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [r*self.ambient, g*self.ambient, b*self.ambient, a])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [self.specular, self.specular, self.specular, a])
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, [0, 0, 0, a])

    def __get_color_index__(self, u):
        ret = int(floor((u - self.min_u)/((self.max_u - self.min_u)/self.num_color))) - 1
        return ret - 1 if ret > 0 else 0

    def show_legend(self):
        font = QFont('Times', 12, QFont.Normal)
        fm = QFontMetrics(font)
        font_w1 = fm.width("12")
        font_w2 = fm.width("1234567890")
        font_h = fm.height()
        start = self.max_u
        stop = self.min_u
        cy = self.rect().top() + 20
        v = start
        n = 8 if self.min_u != self.max_u else 1
        step = (self.max_u - self.min_u)/n
        for k in range(0, n):
            if k == n - 1:
                v = stop
            i = self.__get_color_index__(v)
            glColor3f(self.color_table[i][0], self.color_table[i][1], self.color_table[i][2])
            font.setStyleStrategy(QFont.OpenGLCompatible)
            self.gl.renderText(self.rect().width() - font_w1 - font_w2 - 50, cy, '█', font)
            glColor3f(1, 1, 1)
            self.gl.renderText(self.rect().width() - font_w2 - 50, cy, '{:+0.5E}'.format(v), font)
            cy += font_h
            v -= step

    def __paint_triangle__(self, tri):
        tri = sorted(tri, key=lambda item: item[3])
        ind = []
        for i in range(0, 3):
            ind.append(self.__get_color_index__(tri[i][3]))
        # Треугольник одного цвета
        if ind[0] == ind[1] == ind[2]:
            self.set_color(self.color_table[ind[0]][0], self.color_table[ind[0]][1], self.color_table[ind[0]][2], 1)
            glBegin(GL_TRIANGLES)
            for i in range(0, 3):
                glVertex3f(tri[i][0] - self.x_c[0], tri[i][1] - self.x_c[1], tri[i][2] - self.x_c[2])
            glEnd()
        else:
            # Изолинии проходят по треугольнику
            step = ind[2] - ind[0] + 1
            p02 = []
            x = [tri[0][0], tri[0][1], tri[0][2], ind[0]]
            h = [(tri[2][0] - tri[0][0])/step, (tri[2][1] - tri[0][1])/step, (tri[2][2] - tri[0][2])/step,
                 (ind[2] - ind[0])/step]
            for i in range(0, step):
                p02.append([x[0] + i*h[0], x[1] + i*h[1], x[2] + i*h[2], ind[0] + i*h[3]])
            p02.append([tri[2][0], tri[2][1], tri[2][2], ind[2]])

            step = ind[1] - ind[0] + 1
            p012 = []
            x = [tri[0][0], tri[0][1], tri[0][2], ind[0]]
            h = [(tri[1][0] - tri[0][0])/step, (tri[1][1] - tri[0][1])/step, (tri[1][2] - tri[0][2])/step,
                 (ind[1] - ind[0])/step]
            for i in range(1, step):
                p012.append([x[0] + i*h[0], x[1] + i*h[1], x[2] + i*h[2], ind[0] + i*h[3]])
            p012.append([tri[1][0], tri[1][1], tri[1][2], ind[1]])

            step = ind[2] - ind[1] + 1
            x = [tri[1][0], tri[1][1], tri[1][2], ind[1]]
            h = [(tri[2][0] - tri[1][0])/step, (tri[2][1] - tri[1][1])/step, (tri[2][2] - tri[1][2])/step,
                 (ind[2] - ind[1])/step]
            for i in range(1, step):
                p012.append([x[0] + i*h[0], x[1] + i*h[1], x[2] + i*h[2], ind[1] + i*h[3]])

            for i in range(0, len(p02) - 1):
                if i < len(p012):
                    clr = round((p02[i][3] + p02[i + 1][3] + p012[i][3])/3)
                    self.set_color(self.color_table[clr][0], self.color_table[clr][1], self.color_table[clr][2], 1)
                    glBegin(GL_TRIANGLES)
                    glVertex3f(p02[i][0] - self.x_c[0], p02[i][1] - self.x_c[1], p02[i][2] - self.x_c[2])
                    glVertex3f(p02[i + 1][0] - self.x_c[0], p02[i + 1][1] - self.x_c[1], p02[i + 1][2] - self.x_c[2])
                    glVertex3f(p012[i][0] - self.x_c[0], p012[i][1] - self.x_c[1], p012[i][2] - self.x_c[2])
                    glEnd()
                    if i + 1 < len(p012):
                        clr = round((p02[i + 1][3] + p012[i][3] + p012[i + 1][3])/3)
                        self.set_color(self.color_table[clr][0], self.color_table[clr][1], self.color_table[clr][2], 1)
                        glBegin(GL_TRIANGLES)
                        glVertex3f(p02[i + 1][0] - self.x_c[0], p02[i + 1][1] - self.x_c[1], p02[i + 1][2] -
                                   self.x_c[2])
                        glVertex3f(p012[i][0] - self.x_c[0], p012[i][1] - self.x_c[1], p012[i][2] - self.x_c[2])
                        glVertex3f(p012[i + 1][0] - self.x_c[0], p012[i + 1][1] - self.x_c[1], p012[i + 1][2] -
                                   self.x_c[2])
                        glEnd()


# Визуализация одномерной задачи
class TPlot1d(TOpenGLWidget):
    def __init__(self, x, fe, be, results, index):
        super(TPlot1d, self).__init__(x, fe, be, results, index)
        self.gl.paintGL = self.__paint__

    def __paint__(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        glEnable(GL_POLYGON_OFFSET_FILL)
        glPolygonOffset(1.0, 1.0)
        if self.is_light:
            glDisable(GL_COLOR_MATERIAL)
            glEnable(GL_LIGHTING)
        else:
            glDisable(GL_LIGHTING)
            glEnable(GL_COLOR_MATERIAL)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glNormal3d(1, 0, 0)
        # Изображение КЭ
        for i in range(0, len(self.fe)):
            rod = []
            for j in range(0, len(self.fe[0])):
                rod.append([self.x[self.fe[i][j]][0], self.results[self.fe[i][j]]])
            rod = sorted(rod, key=lambda item: item[1])
            ind = []
            for j in range(0, 2):
                ind.append(self.__get_color_index__(rod[j][1]))
            if ind[0] == ind[1]:
                self.set_color(self.color_table[ind[0]][0], self.color_table[ind[0]][1], self.color_table[ind[0]][2], 1)
                glBegin(GL_LINES)
                glVertex2f(rod[0][0] - self.x_c[0])
                glVertex2f(rod[1][0])
                glEnd()
            else:
                step = abs(ind[1] - ind[0]) + 1
                x = rod[0][0]
                h = (rod[1][0] - rod[0][0])/step
                clr = ind[0]
                for j in range(0, step):
                    self.set_color(self.color_table[clr][0], self.color_table[clr][1], self.color_table[clr][2], 1)
                    glBegin(GL_LINES)
                    glVertex2f(x - self.x_c[0], 0)
                    glVertex2f(x + h - self.x_c[0], 0)
                    glEnd()
                    x += h
                    clr += 1
        # Изображение цветовой шкалы
        if self.is_legend:
            self.show_legend()


# Визуализация плоской задачи
class TPlot2d(TOpenGLWidget):
    def __init__(self, x, fe, be, results, index):
        super(TPlot2d, self).__init__(x, fe, be, results, index)
        self.gl.paintGL = self.__paint__

    def __paint__(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        glEnable(GL_POLYGON_OFFSET_FILL)
        glPolygonOffset(1.0, 1.0)
        if self.is_light:
            glDisable(GL_COLOR_MATERIAL)
            glEnable(GL_LIGHTING)
        else:
            glDisable(GL_LIGHTING)
            glEnable(GL_COLOR_MATERIAL)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glNormal3d(0, 0, 1)
        # Изображение КЭ
        for i in range(0, len(self.fe)):
            tri = []
            for j in range(0, len(self.fe[0])):
                tri.append([self.x[self.fe[i][j]][0], self.x[self.fe[i][j]][1], 0, self.results[self.fe[i][j]]])
            if len(tri) == 3:
                self.__paint_triangle__(tri)
            else:
                self.__paint_triangle__([tri[0], tri[1], tri[2]])
                self.__paint_triangle__([tri[0], tri[3], tri[2]])
            if self.is_fe_border:
                self.__paint_fe_border__(tri)
        # Изображение границы области
        if self.is_object_border:
            for i in range(0, len(self.be)):
                self.set_color(1, 1, 1, 1)
                glBegin(GL_LINE_LOOP)
                glVertex2f(self.x[self.be[i][0]][0] - self.x_c[0], self.x[self.be[i][0]][1] - self.x_c[1])
                glVertex2f(self.x[self.be[i][1]][0] - self.x_c[0], self.x[self.be[i][1]][1] - self.x_c[1])
                glEnd()
        # Изображение цветовой шкалы
        if self.is_legend:
            self.show_legend()

    def __paint_fe_border__(self, tri):
        self.set_color(1, 1, 1, 1)
        glBegin(GL_LINE_LOOP)
        for i in range(0, len(tri)):
            glVertex2f(tri[i][0] - self.x_c[0], tri[i][1] - self.x_c[1])
        glEnd()


# Класс, реализующий визуализацию расчета
class TPlot3d:
    def __init__(self):
        self.__name__ = ''          # Название задачи
        self.__task_type__ = ''     # Тип задачи
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
            print_error('Unable open file ' + file_name)
            return False

        with open(file_name, 'r') as file:
            data = json.loads(json.load(file))
        header = data['header']
        self.__name__ = header['name']
        self.__task_type__ = header['type']
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
            print_error('Error: \'%s\' is not a recognized function name or incorrect time' % fun_name)
            return
        # Задание заголовка
        if self.__task_type__ == 'dynamic':
            fun_name += ' (t = %5.2f)' % t
        # Визуализация результата
        app = QApplication(sys.argv)
        if self.__fe_type__ == 'fe_1d_2':
            window = TPlot1d(self.__x__, self.__fe__, self.__be__, self.__results__, index)
        elif self.__fe_type__ == 'fe_2d_3' or self.__fe_type__ == 'fe_2d_4':
            window = TPlot2d(self.__x__, self.__fe__, self.__be__, self.__results__, index)
        elif self.__fe_type__ == 'fe_3d_4':
            window = TPlot3d4(self.__x__, self.__fe__, self.__be__, self.__results__, index)
        else:
            window = TPlot3d8(self.__x__, self.__fe__, self.__be__, self.__results__, index)
        window.resize(500, 500)
        window.setWindowTitle(fun_name)
        window.show()
        sys.exit(app.exec_())
