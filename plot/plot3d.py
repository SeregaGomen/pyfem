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
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QMessageBox, QVBoxLayout, QAction, QActionGroup,
                             QMenu, QFileDialog)
from PyQt5.QtGui import (QFont, QFontMetrics)
from PyQt5.QtCore import (Qt, QObject, QPoint)
from PyQt5.QtOpenGL import QGLWidget
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *


class TMainWindow(QMainWindow):
    def __init__(self, file_name):
        super().__init__()
        self.file_name = file_name  # Имя файла
        self.fe_type = ''           # Тип КЭ
        self.x = []                 # Координаты узлов
        self.fe = []                # Связи КЭ
        self.be = []                # ... ГЭ
        self.results = []           # Результаты расчета

        self.__gl_widget__ = TGLWidget()
        self.setCentralWidget(self.__gl_widget__)
        self.resize(640, 480)
        # Загрузка данных
        if self.__load_file__():
            self.__gl_widget__.set_data(self.fe_type, self.x, self.fe, self.be, self.results[0].results)
        # Настройка окна
        self.__init_main_menu__()

    # Загрузка данных из файла
    def __load_file__(self):
        # Подготовка имени файла
        if len(self.file_name) < 5 or self.file_name[len(self.file_name) - 5:] != '.json':
            self.file_name += '.json'
        # Проверка наличия файла
        if not os.path.exists(self.file_name):
            print_error('Unable open file ' + self.file_name)
            self.statusBar().showMessage('Error read data from file ' + self.file_name)
            self.setWindowTitle('PyFEM results viewer')
            return False
        self.setWindowTitle('PyFEM results viewer - ' + self.file_name)
        self.statusBar().showMessage('Success load file ' + self.file_name)
        # Чтение файла
        try:
            with open(self.file_name, 'r') as file:
                data = json.loads(json.load(file))
        except IOError:
            print_error('Unable read file ' + self.file_name)
            return False
        # Обработка данных
        mesh = data['mesh']
        self.fe_type = mesh['fe_type']
        self.x = mesh['vertex']
        self.fe = mesh['fe']
        self.be = mesh['be']
        results = data['results']
        for i in range(0, len(results)):
            res = TResult()
            res.name = results[i]['function']
            res.t = results[i]['t']
            res.results = results[i]['results']
            self.results.append(res)
        return True

    def __init_main_menu__(self):
        # self.setWindowTitle(self.file_name)
        # self.statusBar().showMessage('')

        main_menu = self.menuBar()
        file_menu = main_menu.addMenu('&File')
        function_menu = main_menu.addMenu('F&unction')
        options_menu = main_menu.addMenu('&Option')
        help_menu = main_menu.addMenu('&?')

        # Настройка File
        open_action = QAction('&Open...', self)
        open_action.setStatusTip('Open a data file')
        open_action.triggered.connect(self.__open__)
        file_menu.addAction(open_action)
        close_action = QAction('&Close', self)
        close_action.setStatusTip('Close current file')
        close_action.triggered.connect(self.__close__)
        file_menu.addAction(close_action)
        file_menu.addSeparator()
        exit_action = QAction('E&xit', self)
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Настройка Function
        qa = QActionGroup(self)
        for i in range(0, len(self.results)):
            if i < 9:
                fun = '&' + str(i + 1) + ' ' + self.results[i].name
            else:
                fun = '&' + chr(ord('A') + i - 9) + ' ' + self.results[i].name
            function_action = QAction(fun, self)
            function_action.setStatusTip('Visualization function ' + self.results[i].name)
            function_action.setActionGroup(qa)
            function_action.setCheckable(True)
            function_action.setChecked(True if i == 0 else False)
            function_action.triggered.connect(self.__fun_action__)
            function_menu.addAction(function_action)

        # Настройка Options
        light_action = QAction('&Light', self)
        light_action.setStatusTip('Enable light')
        light_action.triggered.connect(self.__light_action__)
        light_action.setCheckable(True)
        light_action.setChecked(True)
        options_menu.addAction(light_action)
        fe_border_action = QAction('&FE border', self)
        fe_border_action.setStatusTip('Enable drawing FE border')
        fe_border_action.triggered.connect(self.__fe_border_action__)
        fe_border_action.setCheckable(True)
        fe_border_action.setChecked(False)
        options_menu.addAction(fe_border_action)
        legend_action = QAction('&Legend', self)
        legend_action.setStatusTip('Enable drawing legend')
        legend_action.triggered.connect(self.__legend_action__)
        legend_action.setCheckable(True)
        legend_action.setChecked(True)
        options_menu.addAction(legend_action)
        options_menu.addSeparator()
        qa = QActionGroup(self)
        color_16_action = QAction('&16', self)
        color_16_action.setStatusTip('Set 16 grading colors')
        color_16_action.setActionGroup(qa)
        color_16_action.setCheckable(True)
        color_16_action.setChecked(True)
        color_16_action.triggered.connect(self.__colors_action__)
        color_32_action = QAction('3&2', self)
        color_32_action.setStatusTip('Set 32 grading colors')
        color_32_action.setActionGroup(qa)
        color_32_action.setCheckable(True)
        color_32_action.setChecked(False)
        color_32_action.triggered.connect(self.__colors_action__)
        color_64_action = QAction('6&4', self)
        color_64_action.setStatusTip('Set 64 grading colors')
        color_64_action.setActionGroup(qa)
        color_64_action.setCheckable(True)
        color_64_action.setChecked(False)
        color_64_action.triggered.connect(self.__colors_action__)
        color_128_action = QAction('1&28', self)
        color_128_action.setStatusTip('Set 128 grading colors')
        color_128_action.setActionGroup(qa)
        color_128_action.setCheckable(True)
        color_128_action.setChecked(False)
        color_128_action.triggered.connect(self.__colors_action__)
        color_256_action = QAction('25&6', self)
        color_256_action.setStatusTip('Set 256 grading colors')
        color_256_action.setActionGroup(qa)
        color_256_action.setCheckable(True)
        color_256_action.setChecked(False)
        color_256_action.triggered.connect(self.__colors_action__)

        color_menu = QMenu('&Colors', self)
        color_menu.setStatusTip('Set grading colors')
        color_menu.addAction(color_16_action)
        color_menu.addAction(color_32_action)
        color_menu.addAction(color_64_action)
        color_menu.addAction(color_128_action)
        color_menu.addAction(color_256_action)
        options_menu.addMenu(color_menu)

        self.show()

    def __open__(self):
        dlg = QFileDialog(self, 'Open data file', '', 'JSON data files (*.json)')
        if dlg.exec_():
            self.__set_file__(dlg.selectedFiles()[0])

    def __close__(self):
        self.setWindowTitle('PyFEM results viewer')
        self.file_name = ''
        self.fe_type = ''
        self.x.clear()
        self.fe.clear()
        self.be.clear()
        self.results.clear()
        self.__gl_widget__.clear()

    def __colors_action__(self):
        num = int(QObject.sender(self).text().replace('&', ''))
        self.__gl_widget__.set_colors(num)

    def __fun_action__(self):
        index = self.__get_fun_index__(QObject.sender(self).text()[3:])
        self.__gl_widget__.set_results(self.results[index].results)

    def __light_action__(self):
        self.__gl_widget__.trigger_light()
        self.repaint()

    def __fe_border_action__(self):
        self.__gl_widget__.trigger_fe_border()
        self.repaint()

    def __legend_action__(self):
        self.__gl_widget__.trigger_legend()
        self.repaint()

    # Поиск индекса функции по ее имени
    def __get_fun_index__(self, fun_name, t=0):
        # Поиск индекса функции в списке результатов
        index = -1
        for i in range(0, len(self.results)):
            if self.results[i].name == fun_name and self.results[i].t == t:
                index = i
                break
        return index

    # Задание нового файла
    def __set_file__(self, file_name):
        self.file_name = file_name
        self.fe_type = ''
        self.x.clear()
        self.fe.clear()
        self.be.clear()
        self.results.clear()
        # Загрузка данных
        if self.__load_file__() is False:
            QMessageBox.critical(self, 'Error', 'Error read data from file ' + file_name, QMessageBox.Ok)
            return
        self.__gl_widget__.set_data(self.fe_type, self.x, self.fe, self.be, self.results[0].results)


# Базовый класс, реализующий основной функционал OpenGL
class TGLWidget(QWidget):
    def __init__(self):
        super(TGLWidget, self).__init__()
        self.fe_type = ''
        self.x = []
        self.fe = []
        self.be = []
        self.results = []
        self.min_x = []
        self.max_x = []
        self.x_c = []
        self.radius = 0
        self.min_u = 0
        self.max_u = 0
        self.is_light = True
        self.is_legend = True
        self.is_fe_border = False
        self.angle_x = 0
        self.angle_y = 0
        self.angle_z = 0
        self.num_color = 16
        self.__is_idle__ = True
        self.__last_pos__ = QPoint()
        self.__color_table__ = []
        self.__gl__ = QGLWidget(self)
        self.__gl__.initializeGL()
        self.__gl__.resizeGL = self.__resize__
        self.__gl__.paintGL = self.__paint__
        self.__init_color_table__()
        QVBoxLayout(self).addWidget(self.__gl__)
        self.mousePressEvent = self.__mouse_press_event
        self.mouseReleaseEvent = self.__mouse_release_event
        self.__gl__.mouseMoveEvent = self.__mouse_move__
        self.wheelEvent = self.__wheel_event__
        self.__xlist_object__ = 0
        self.__xlist_sceleton__ = 0

    """def __del__(self):
        print('Delete')
        if self.__xlist_object__ != 0:
            glDeleteLists(self.__xlist_object__, 1)
        if self.__xlist_sceleton__ != 0:
            glDeleteLists(self.__xlist_sceleton__, 1)"""

    def clear(self):
        self.fe_type = ''
        self.x.clear()
        self.fe.clear()
        self.be.clear()
        self.results.clear()
        self.min_x.clear()
        self.max_x.clear()
        self.x_c.clear()
        self.__color_table__.clear()
        if self.__xlist_object__ != 0:
            glDeleteLists(self.__xlist_object__, 1)
            self.__xlist_object__ = 0
        if self.__xlist_sceleton__ != 0:
            glDeleteLists(self.__xlist_sceleton__, 1)
            self.__xlist_sceleton__ = 0
        self.__gl__.updateGL()

    def set_data(self, fe_type, x, fe, be, results):
        self.__gl__.glInit()
        self.fe_type = fe_type
        self.x = x
        self.fe = fe
        self.be = be
        self.results = results
        self.min_x, self.max_x, self.x_c, self.radius = self.__get_coord_info__()
        self.min_u = min(self.results)
        self.max_u = max(self.results)
        self.is_light = True
        self.is_legend = True
        self.is_fe_border = False
        self.angle_x = 0
        self.angle_y = 0
        self.angle_z = 0
        self.scale = 1
        self.num_color = 16
        self.__is_idle__ = True
        self.__last_pos__ = QPoint()
        self.__color_table__ = []
        self.__init_color_table__()
        if self.__xlist_object__ != 0:
            glDeleteLists(self.__xlist_object__, 1)
            self.__xlist_object__ = 0
        if self.__xlist_sceleton__ != 0:
            glDeleteLists(self.__xlist_sceleton__, 1)
            self.__xlist_sceleton__ = 0
        self.__resize__(self.width(), self.height())
        self.__gl__.updateGL()

    def redraw(self):
        if self.__xlist_object__ != 0:
            glDeleteLists(self.__xlist_object__, 1)
            self.__xlist_object__ = 0
        if self.__xlist_sceleton__ != 0:
            glDeleteLists(self.__xlist_sceleton__, 1)
            self.__xlist_sceleton__ = 0
        self.__gl__.updateGL()

    def trigger_light(self):
        self.is_light = not self.is_light
        self.redraw()
        self.__gl__.update()

    def trigger_fe_border(self):
        self.is_fe_border = not self.is_fe_border
        self.redraw()
        self.__gl__.update()

    def trigger_legend(self):
        self.is_legend = not self.is_legend
        self.redraw()
        self.__gl__.update()

    def set_colors(self, colors):
        self.num_color = colors
        self.__init_color_table__()
        self.redraw()
        self.__gl__.update()

    def set_results(self, results):
        self.results = results
        self.min_u = min(self.results)
        self.max_u = max(self.results)
        self.__init_color_table__()
        self.redraw()
        self.__gl__.update()

    def __mouse_press_event(self, event):
        super(TGLWidget, self).mousePressEvent(event)
        if event.buttons() & Qt.LeftButton:
            self.__is_idle__ = False

    def __mouse_release_event(self, event):
        super(TGLWidget, self).mouseReleaseEvent(event)
        if self.__is_idle__ is False:
            self.__is_idle__ = True
            self.__gl__.repaint()

    def __mouse_move__(self, event):
        dx = event.x() - self.__last_pos__.x()
        dy = event.y() - self.__last_pos__.y()
        self.__last_pos__ = event.pos()
        if event.buttons() & Qt.LeftButton:
            if event.modifiers() & Qt.ShiftModifier:
                self.angle_z += (dy/abs(dy) if dy != 0 else 0) + (dx/abs(dx) if dx != 0 else 0)
            else:
                self.angle_x += dy/abs(dy) if dy != 0 else 0
                self.angle_y += dx/abs(dx) if dx != 0 else 0
            self.__gl__.repaint()

    def __wheel_event__(self, event):
        if event.angleDelta().y() > 0:
            self.scale *= 1.05
        else:
            self.scale /= 1.05
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
        self.__color_table__.clear()
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

    def color(self, i):
        if self.is_light:
            self.__make_material__(self.__color_table__[i][0], self.__color_table__[i][1], self.__color_table__[i][2])
        else:
            glColor3f(self.__color_table__[i][0], self.__color_table__[i][1], self.__color_table__[i][2])

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

    @staticmethod
    def __make_material__(r, g, b, a=1.0):
        diffuse = 0.5
        ambient = 0.4
        specular = 0.7
        shininess = 50
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, [r*diffuse, g*diffuse, b*diffuse, a])
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [r*ambient, g*ambient, b*ambient, a])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [specular, specular, specular, a])
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess)

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
        if self.is_light:
            # Задание нормали
            a = (tri[1][1] - tri[0][1])*(tri[2][2] - tri[0][2]) - (tri[2][1] - tri[0][1])*(tri[1][2] - tri[0][2])
            b = (tri[2][0] - tri[0][0])*(tri[1][2] - tri[0][2]) - (tri[1][0] - tri[0][0])*(tri[2][2] - tri[0][2])
            c = (tri[1][0] - tri[0][0])*(tri[2][1] - tri[0][1]) - (tri[2][0] - tri[0][0])*(tri[1][1] - tri[0][1])
            glNormal3f(a, b, c)
        if ind[0] == ind[1] == ind[2]:
            # Треугольник одного цвета
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
                    clr = round(min(p02[i][3], p02[i + 1][3], p012[i][3]))
                    self.color(clr)
                    glBegin(GL_TRIANGLES)
                    glVertex3f(p02[i][0] - self.x_c[0], p02[i][1] - self.x_c[1], p02[i][2] - self.x_c[2])
                    glVertex3f(p012[i][0] - self.x_c[0], p012[i][1] - self.x_c[1], p012[i][2] - self.x_c[2])
                    glVertex3f(p02[i + 1][0] - self.x_c[0], p02[i + 1][1] - self.x_c[1], p02[i + 1][2] - self.x_c[2])
                    glEnd()
                    if i + 1 < len(p012):
                        clr = round(min(p02[i + 1][3], p012[i][3], p012[i + 1][3]))
                        self.color(clr)
                        glBegin(GL_TRIANGLES)
                        glVertex3f(p02[i + 1][0] - self.x_c[0], p02[i + 1][1] - self.x_c[1], p02[i + 1][2] -
                                   self.x_c[2])
                        glVertex3f(p012[i][0] - self.x_c[0], p012[i][1] - self.x_c[1], p012[i][2] - self.x_c[2])
                        glVertex3f(p012[i + 1][0] - self.x_c[0], p012[i + 1][1] - self.x_c[1], p012[i + 1][2] -
                                   self.x_c[2])
                        glEnd()

    @staticmethod
    def __setup_light__():
        diffuse = 0.8
        ambient = 0.8
        specular = 0.6

        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
        glLightfv(GL_LIGHT0, GL_AMBIENT, [ambient, ambient, ambient])
        glLightfv(GL_LIGHT0, GL_DIFFUSE, [diffuse, diffuse, diffuse])
        glLightfv(GL_LIGHT0, GL_SPECULAR, [specular, specular, specular])

        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_NORMALIZE)

    def draw_fe_border(self, tri):
        glDisable(GL_LIGHTING)
        glEnable(GL_COLOR_MATERIAL)
        glColor3f(0, 0, 0)
        glBegin(GL_LINE_LOOP)
        for i in range(0, len(tri)):
            glVertex3f(tri[i][0] - self.x_c[0], tri[i][1] - self.x_c[1], tri[i][2] - self.x_c[2])
        glEnd()
        if self.is_light:
            glDisable(GL_COLOR_MATERIAL)
            glEnable(GL_LIGHTING)

    def __paint__(self):
        glClearColor(0.39, 0.39, 0.6, 0.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        if self.fe_type == '':
            return
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

        glScalef(self.scale, self.scale, self.scale)

        if self.__is_idle__:
            self.__display_object__()
        else:
            self.__display_sceleton__()

        glPopMatrix()
        # Изображение цветовой шкалы
        if self.is_legend and self.__is_idle__:
            self.show_legend()

    # Визуализация результата
    def __display_object__(self):
        if self.__xlist_object__ == 0:
            self.__xlist_object__ = glGenLists(1)
            glNewList(self.__xlist_object__, GL_COMPILE)
            if self.fe_type == 'fe_1d_2':
                self.__paint_1d__()
            elif self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_4':
                self.__paint_2d__()
            else:
                self.__paint_3d__()
            glEndList()
        else:
            glCallList(self.__xlist_object__)

    # Изображение каркаса объекта
    def __display_sceleton__(self):
        if self.__xlist_sceleton__ == 0:
            self.__xlist_sceleton__ = glGenLists(2)
            glNewList(self.__xlist_sceleton__, GL_COMPILE)
            if self.fe_type == 'fe_1d_2':
                glBegin(GL_LINES)
                glVertex2f(self.min_x[0] - self.x_c[0], 0)
                glVertex2f(self.max_x[0] - self.x_c[0], 0)
                glEnd()
            elif self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_4':
                for i in range(0, len(self.be)):
                    glBegin(GL_LINES)
                    for j in range(0, len(self.be[0])):
                        glVertex2f(self.x[self.be[i][j]][0] - self.x_c[0], self.x[self.be[i][j]][1] - self.x_c[1])
                    glEnd()
            else:
                for i in range(0, len(self.be)):
                    glBegin(GL_LINE_LOOP)
                    for j in range(0, len(self.be[0])):
                        glVertex3f(self.x[self.be[i][j]][0] - self.x_c[0], self.x[self.be[i][j]][1] - self.x_c[1],
                                   self.x[self.be[i][j]][2] - self.x_c[2])
                    glEnd()
            glEndList()
        else:
            glCallList(self.__xlist_sceleton__)

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
                glVertex2f(rod[0][0] - self.x_c[0], 0)
                glVertex2f(rod[1][0], 0)
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
#        for i in range(212, 214):
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

# Класс, реализующий визуализацию расчета
class TPlot:
    def __init__(self, file_name):
        self.file_name = file_name
        # Создание главного окна приложения
        app = QApplication(sys.argv)
        window = TMainWindow(file_name)
        window.show()
        sys.exit(app.exec_())


""" 1. При открытии нового объекта сбрасывать меню параметров изображения в начальное состоянии """
""" 2. Добавить кнопку Ctrl при вращении объекта (чтобы можно было крутить по каждой оси отдельно) """
""" 3. Сделать translate """
""" 4. Запись в файл """