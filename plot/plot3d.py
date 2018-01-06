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
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QAction)
from PyQt5.QtGui import (QFont, QFontMetrics, QIcon)
from PyQt5.QtOpenGL import QGLWidget
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *


class TApp(QMainWindow):
    def __init__(self, name, fe_type, x, fe, be, results):
        super().__init__()
        self.title = name
        self.resize(640, 480)
        self.__init_window__()
        self.setCentralWidget(TGLWidget(fe_type, x, fe, be, results, 0))

    def __init_window__(self):
        self.setWindowTitle(self.title)
#        self.setGeometry(self.left, self.top, self.width, self.height)

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        editMenu = mainMenu.addMenu('Edit')
        viewMenu = mainMenu.addMenu('View')
        searchMenu = mainMenu.addMenu('Search')
        toolsMenu = mainMenu.addMenu('Tools')
        helpMenu = mainMenu.addMenu('Help')

        exitButton = QAction(QIcon('exit24.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)
        self.show()

# Базовый класс, реализующий основной функционал OpenGL
class TGLWidget(QWidget):
    def __init__(self, fe_type, x, fe, be, results, index):
        super(TGLWidget, self).__init__()
        self.fe_type = fe_type
        self.x = x
        self.fe = fe
        self.be = be
        self.results = results[index].results
        self.min_x, self.max_x, self.x_c, self.radius = self.__get_coord_info__()
        self.min_u = min(self.results)
        self.max_u = max(self.results)
        self.is_light = False
        self.is_legend = True
        self.is_fe_border = True
        self.angle_x = 0
        self.angle_y = 0
        self.angle_z = 0
        self.alpha = 1.0
        self.diffuse = 0.8
        self.ambient = 0.8
        self.specular = 0.6
        self.num_color = 16
        self.__color_table__ = []
        self.__gl__ = QGLWidget(self)
        self.__gl__.initializeGL()
        self.__gl__.resizeGL = self.__resize__
        self.__gl__.paintGL = self.__paint__
        self.__init_color_table__()
        QVBoxLayout(self).addWidget(self.__gl__)
        self.mousePressEvent = self.__mouse_press_event
        self.__gl__.mouseMoveEvent = self.__mouse_move__
        # Создание нормалей
        if self.fe_type == 'fe_3d_4' or self.fe_type == 'fe_3d_8':
            self.__create_normal__()

    def __mouse_press_event(self, event):
        self.__last_pos__ = event.pos()

    def __mouse_move__(self, event):
        dx = event.x() - self.__last_pos__.x()
        dy = event.x() - self.__last_pos__.y()

        self.angle_x += dx/2
        self.angle_y += dy/2
        self.__gl__.repaint()

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
                self.__color_table__.append([red, 0, 1, u])
                red -= h
                if red < 0:
                    red = 0
            elif step <= i < 2*step:
                # синий-голубой
                self.__color_table__.append([0, green, 1, u])
                green += h
                if green > 1:
                    green = 1
            elif 2*step <= i < 3*step:
                # голубой-зеленый
                self.__color_table__.append([0, 1, blue, u])
                blue -= h
                if blue < 0:
                    blue = 0
            elif 3*step <= i < 4*step:
                # зеленый-желтый
                self.__color_table__.append([red, 1, 0, u])
                red += h
                if red > 1:
                    red = 1
            elif i > 4*step:
                # желтый-оранжевый-красный
                self.__color_table__.append([1, green, 0, u])
                green -= 0.5*h
                if green < 0:
                    green = 0
            u += h_u

    def color(self, index):
        r = self.__color_table__[index][0] if index >= 0 else 1
        g = self.__color_table__[index][1] if index >= 0 else 1
        b = self.__color_table__[index][2] if index >= 0 else 1
        if self.is_light:
            self.__make_material__(r, g, b, self.alpha)
        else:
            glColor4f(r, g, b, self.alpha)

    def __resize__(self, w, h):
        aspect = w/h
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(60.0, aspect, 0.01*self.radius, 10*self.radius)
        gluLookAt(0, 0, self.radius, 0, 0, 0, 0, 1, 0)
        glViewport(0, 0, w, h)
        glMatrixMode(GL_MODELVIEW)
        if self.is_light:
            self.__setup_light__()

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
            glColor3f(self.__color_table__[i][0], self.__color_table__[i][1], self.__color_table__[i][2])
            font.setStyleStrategy(QFont.OpenGLCompatible)
            self.__gl__.renderText(self.rect().width() - font_w1 - font_w2 - 50, cy, '█', font)
            glColor3f(1, 1, 1)
            self.__gl__.renderText(self.rect().width() - font_w2 - 50, cy, '{:+0.5E}'.format(v), font)
            cy += font_h
            v -= step

    def draw_triangle_3d(self, tri):
        tri = sorted(tri, key=lambda item: item[3])
        ind = []
        for i in range(0, 3):
            ind.append(self.__get_color_index__(tri[i][3]))
        # Треугольник одного цвета
        if ind[0] == ind[1] == ind[2]:
            self.color(ind[0])
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
                    self.color(clr)
                    glBegin(GL_TRIANGLES)
                    glVertex3f(p02[i][0] - self.x_c[0], p02[i][1] - self.x_c[1], p02[i][2] - self.x_c[2])
                    glVertex3f(p02[i + 1][0] - self.x_c[0], p02[i + 1][1] - self.x_c[1], p02[i + 1][2] - self.x_c[2])
                    glVertex3f(p012[i][0] - self.x_c[0], p012[i][1] - self.x_c[1], p012[i][2] - self.x_c[2])
                    glEnd()
                    if i + 1 < len(p012):
                        clr = round((p02[i + 1][3] + p012[i][3] + p012[i + 1][3])/3)
                        self.color(clr)
                        glBegin(GL_TRIANGLES)
                        glVertex3f(p02[i + 1][0] - self.x_c[0], p02[i + 1][1] - self.x_c[1], p02[i + 1][2] -
                                   self.x_c[2])
                        glVertex3f(p012[i][0] - self.x_c[0], p012[i][1] - self.x_c[1], p012[i][2] - self.x_c[2])
                        glVertex3f(p012[i + 1][0] - self.x_c[0], p012[i + 1][1] - self.x_c[1], p012[i + 1][2] -
                                   self.x_c[2])
                        glEnd()

    def __setup_light__(self):
        glLightfv(GL_LIGHT0, GL_AMBIENT, [self.ambient, self.ambient, self.ambient])
        glLightfv(GL_LIGHT0, GL_DIFFUSE, [self.diffuse, self.diffuse, self.diffuse])
        glLightfv(GL_LIGHT0, GL_SPECULAR, [self.specular, self.specular, self.specular])
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50)
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [1, 1, 1, 1])
        glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE)
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)

    def draw_fe_border(self, tri):
        glDisable(GL_LIGHTING)
        glEnable(GL_COLOR_MATERIAL)
        glColor4f(0, 0, 0, self.alpha)
        glBegin(GL_LINE_LOOP)
        for i in range(0, len(tri)):
            glVertex3f(tri[i][0] - self.x_c[0], tri[i][1] - self.x_c[1], tri[i][2] - self.x_c[2])
        glEnd()
        if self.is_light:
            glDisable(GL_COLOR_MATERIAL)
            glEnable(GL_LIGHTING)

    def __paint__(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glEnable(GL_DEPTH_TEST)
        glLoadIdentity()
        glEnable(GL_POLYGON_OFFSET_FILL)
        glPolygonOffset(1.0, 1.0)
        if self.is_light:
            glDisable(GL_COLOR_MATERIAL)
            glEnable(GL_LIGHTING)
        else:
            glDisable(GL_LIGHTING)
            glEnable(GL_COLOR_MATERIAL)

        glPushMatrix()
        glRotatef(self.angle_x, 1, 0, 0)
        glRotatef(self.angle_y, 0, 1, 0)
        glRotatef(self.angle_z, 0, 0, 1)

        self.__display_object__()
        glPopMatrix()
        # Изображение цветовой шкалы
        if self.is_legend:
            self.show_legend()

    # Визуализация результата
    def __display_object__(self):
        if self.fe_type == 'fe_1d_2':
            self.__paint_1d__()
        elif self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_4':
            self.__paint_2d__()
        else:
            self.__paint_3d__()

    # Визуализация одномерной задачи
    def __paint_1d__(self):
        # Изображение КЭ
        for i in range(0, len(self.fe)):
            rod = []
            for j in range(0, len(self.fe[0])):
                rod.append([self.x[self.fe[i][j]][0], self.results[self.fe[i][j]]])
            rod = sorted(rod, key=lambda item: item[1])
            ind = []
            for j in range(0, 2):
                ind.append(self.__get_color_index__(rod[j][1]))
            clr = ind[0]
            if ind[0] == ind[1]:
                self.color(clr)
                glBegin(GL_LINES)
                glVertex2f(rod[0][0] - self.x_c[0])
                glVertex2f(rod[1][0])
                glEnd()
            else:
                step = abs(ind[1] - ind[0]) + 1
                x = rod[0][0]
                h = (rod[1][0] - rod[0][0])/step
                for j in range(0, step):
                    self.color(clr)
                    glBegin(GL_LINES)
                    glVertex2f(x - self.x_c[0], 0)
                    glVertex2f(x + h - self.x_c[0], 0)
                    glEnd()
                    x += h
                    clr += 1

    # Визуализация плоской задачи
    def __paint_2d__(self):
        # Изображение КЭ
        for i in range(0, len(self.fe)):
            tri = []
            for j in range(0, len(self.fe[0])):
                tri.append([self.x[self.fe[i][j]][0], self.x[self.fe[i][j]][1], 0, self.results[self.fe[i][j]]])
            if len(tri) == 3:
                self.draw_triangle_3d(tri)
            else:
                self.draw_triangle_3d([tri[0], tri[1], tri[2]])
                self.draw_triangle_3d([tri[0], tri[3], tri[2]])
            if self.is_fe_border:
                self.draw_fe_border(tri)

    # Визуализация пространственной задачи
    def __paint_3d__(self):
        # Изображение поверхности
        for i in range(0, len(self.be)):
            if self.is_light:
                glNormal3d(self.__normal__[i][0], self.__normal__[i][1], self.__normal__[i][2])
            tri = []
            for j in range(0, len(self.be[0])):
                tri.append([self.x[self.be[i][j]][0], self.x[self.be[i][j]][1], self.x[self.be[i][j]][2],
                            self.results[self.be[i][j]]])
            if len(tri) == 3:
                self.draw_triangle_3d(tri)
            else:
                self.draw_triangle_3d([tri[0], tri[1], tri[2]])
                self.draw_triangle_3d([tri[0], tri[3], tri[2]])
            # Изображение границы ГЭ
            if self.is_fe_border:
                self.draw_fe_border(tri)

    def __create_normal__(self):
        normal = []
        for i in range(0, len(self.be)):
            x = []
            for j in range(0, 3):
                x.append([self.x[self.be[i][j]][0], self.x[self.be[i][j]][1], self.x[self.be[i][j]][2]])
            a = (x[1][1] - x[0][1])*(x[2][2] - x[0][2]) - (x[2][1] - x[0][1])*(x[1][2] - x[0][2])
            b = (x[2][0] - x[0][0])*(x[1][2] - x[0][2]) - (x[1][0] - x[0][0])*(x[2][2] - x[0][2])
            c = (x[1][0] - x[0][0])*(x[2][1] - x[0][1]) - (x[2][0] - x[0][0])*(x[1][1] - x[0][1])
            d = (a**2 + b**2 + c**2)**0.5
            if d == 0:
                d = 1
            normal.append([a/d, b/d, c/d])
        return normal


# Класс, реализующий визуализацию расчета
class TPlot:
    def __init__(self, file_name):
        self.__name__ = ''          # Название задачи
        self.__task_type__ = ''     # Тип задачи
        self.__date_time__ = ''     # Время и дата рачета задачи
        self.__fe_type__ = ''       # Тип КЭ
        self.__x__ = []             # Координаты узлов
        self.__fe__ = []            # Связи КЭ
        self.__be__ = []            # ... ГЭ
        self.__results__ = []       # Результаты расчета

        # Считывание файла
        if self.__load_file__(file_name):
            self.__create_window__()

    # Загрузка данных из файла
    def __load_file__(self, file_name):
        # Подготовка имени файла
        if len(file_name) < 5 or file_name[len(file_name) - 5:] != '.json':
            file_name += '.json'
        # Проверка наличия файла
        if not os.path.exists(file_name):
            print_error('Unable open file ' + file_name)
            return False
        # Чтение файла
        try:
            with open(file_name, 'r') as file:
                data = json.loads(json.load(file))
        except IOError:
            print_error('Unable read file ' + file_name)
            return False
        # Обработка данных
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

    # Создание главного окна
    def __create_window__(self):
        # Создание окна
        app = QApplication(sys.argv)
        window = TApp(self.__name__, self.__fe_type__, self.__x__, self.__fe__, self.__be__, self.__results__)
        window.show()
        sys.exit(app.exec_())


