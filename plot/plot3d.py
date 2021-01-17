#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#     Визуализация результатов расчета с использованием OpenGL
###################################################################

import os.path
import simplejson as json
from math import floor
from PyQt5.QtCore import Qt, QObject, QPoint, pyqtSlot
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QMessageBox, QVBoxLayout, QAction, QActionGroup, \
    QMenu, QFileDialog, QDialog, QDialogButtonBox, QListWidget, QLabel
from PyQt5.QtGui import QFont, QFontMetrics
from PyQt5.QtOpenGL import QGLWidget
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from core.fem_result import TResult
from core.fem_object import print_error


class TFunctionListDialog(QDialog):
    def __init__(self, results, index):
        super().__init__()
        self.fun_list = QListWidget()   # Список функций для визуализации
        self.current_index = index      # Текущий номер функции в списке

        self.setWindowTitle('Choose function')
        self.resize(250, 200)
        self.setMaximumSize(250, 200)
        self.setMinimumSize(250, 200)
        # Настройка кнопок
        dlg_btn_box= QDialogButtonBox()
        dlg_btn_box.setStandardButtons(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        dlg_btn_box.accepted.connect(self.ok_click)
        dlg_btn_box.rejected.connect(self.close)

        self.fun_list.doubleClicked.connect(self.ok_click)

        # Создание списка функций
        task_type = 'static' if results[len(results) - 1].t == 0 else 'dynamic'
        for i in range(len(results)):
            fun_name = results[i].name
            if task_type == 'dynamic':
                fun_name += ('(' + str(results[i].t) + ')')
            self.fun_list.addItem(fun_name)
        self.fun_list.item(self.current_index).setSelected(True)

        main_layout = QVBoxLayout()
        main_layout.addWidget(self.fun_list)
        main_layout.addWidget(dlg_btn_box)
        self.setLayout(main_layout)

    @pyqtSlot()
    def ok_click(self):
        self.current_index = self.fun_list.currentRow()
        self.accept()


class TMainWindow(QMainWindow):
    def __init__(self, file_name):
        super().__init__()
        self.file_name = file_name  # Имя файла
        self.fe_type = ''           # Тип КЭ
        self.x = []                 # Координаты узлов
        self.fe = []                # Связи КЭ
        self.be = []                # ... ГЭ
        self.results = []           # Результаты расчета
        self.current_index = 0      # Номер текущей визуализируемой функции

        self.__gl_widget = TGLWidget()
        self.setCentralWidget(self.__gl_widget)
        self.resize(640, 480)
        # Загрузка данных
        if self.__load_file():
            self.__gl_widget.set_data(self.fe_type, self.x, self.fe, self.be, self.results, 0)
        # Настройка окна
        self.__init_main_menu()

    # Загрузка данных из файла
    def __load_file(self):
        ret = False
        ext = os.path.splitext(self.file_name)[1]
        if ext == '':
            # По-уморлчанию работаем с json-файлом
            self.file_name += '.json'
            ext = '.json'
        if ext == '.json':
            ret = self.__load_file_json()
        elif ext == '.qres':
            ret = self.__load_file_qres()
        elif ext == '.res':
            ret = self.__load_file_res()
        return ret

    # Загрузка данных из json-файла
    def __load_file_json(self):
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

    # Загрузка данных из qres-файла
    def __load_file_qres(self):
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
            file = open(self.file_name)
            lines = file.readlines()
            file.close()
        except IOError:
            print_error('Unable read file ' + self.file_name)
            return False
        # Обработка данных
        if lines[0].split('\n')[0] != 'QFEM results file':
            print_error('Wrong file format: ' + self.file_name)
            return False
        # Тип КЭ
        type = int(lines[1])
        if type == 3:
            self.fe_type = 'fe_2d_3'
        elif type == 4:
            self.fe_type = 'fe_3d_4'
        elif type == 6:
            self.fe_type = 'fe_2d_6'
        elif type == 8:
            self.fe_type = 'fe_3d_8'
        elif type == 10:
            self.fe_type = 'fe_3d_10'
        elif type == 24:
            self.fe_type = 'fe_2d_4'
        elif type == 34:
            self.fe_type = 'fe_1d_2'
        elif type == 123:
            self.fe_type = 'fe_2d_3_p'
        elif type == 124:
            self.fe_type = 'fe_2d_4_p'
        elif type == 125:
            self.fe_type = 'fe_2d_6_p'
        elif type == 223:
            self.fe_type = 'fe_2d_3_s'
        elif type == 224:
            self.fe_type = 'fe_2d_4_s'
        elif type == 225:
            self.fe_type = 'fe_2d_6_s'
        else:
            print_error('Wrong file format: ' + self.file_name)
            return False
        # Кол-во узлов
        n = int(lines[2])
        # Считываем узлы
        index = 3
        for i in range(n):
            row = list()
            coord = lines[3 + i].split()
            for j in range(len(coord)):
                row.append(float(coord[j]))
            self.x.append(row)
            index += 1
        # Кол-во КЭ
        n = int(lines[index])
        index += 1
        # Считываем КЭ
        for i in range(n):
            row = []
            fe = lines[index].split()
            for j in range(len(fe)):
                row.append(int(fe[j]))
            self.fe.append(row)
            index += 1
        # Кол-во ГЭ
        n = int(lines[index])
        index += 1
        # Считываем ГЭ
        for i in range(n):
            row = []
            be = lines[index].split()
            for j in range(len(be)):
                row.append(int(be[j]))
            self.be.append(row)
            index += 1
        if self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_2d_3_s' or \
                self.fe_type == 'fe_2d_4_p' or self.fe_type == 'fe_2d_6_s' or self.fe_type == 'fe_2d_4_s':
            self.be = self.fe
        # Пропускаекм время расчета
        index += 1
        # Считываем кол-во функций
        n = int(lines[index])
        index += 1
        for i in range(n):
            res = TResult()
            res.name = lines[index].split('\n')[0]
            res.t = float(lines[index + 1])
            m = int(lines[index + 2])
            index += 3
            for j in range(m):
                res.results.append(float(lines[index]))
                index += 1
            self.results.append(res)
        return True

    # Загрузка данных из res-файла
    def __load_file_res(self):
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
            file = open(self.file_name)
            lines = file.readlines()
            file.close()
        except IOError:
            print_error('Unable read file ' + self.file_name)
            return False
        # Обработка данных
        if lines[0].split('\n')[0] != 'FEM Solver Results File':
            print_error('Wrong file format: ' + self.file_name)
            return False
        if lines[1].split('\n')[0] != 'Mesh':
            print_error('Wrong file format: ' + self.file_name)
            return False
        # Тип КЭ
        type = lines[2].split('\n')[0]
        if type == 'fe1d2':
            self.fe_type = 'fe_1d_2'
        elif type == 'fe2d3':
            self.fe_type = 'fe_2d_3'
        elif type == 'fe2d4':
            self.fe_type = 'fe_2d_4'
        elif type == 'fe2d6':
            self.fe_type = 'fe_2d_6'
        elif type == 'fe2d3p':
            self.fe_type = 'fe_2d_3_p'
        elif type == 'fe2d4p':
            self.fe_type = 'fe_2d_4_p'
        elif type == 'fe2d6p':
            self.fe_type = 'fe_2d_6_p'
        elif type == 'fe3d4':
            self.fe_type = 'fe_3d_4'
        elif type == 'fe3d8':
            self.fe_type = 'fe_3d_8'
        elif type == 'fe3d10':
            self.fe_type = 'fe_3d_10'
        elif type == 'fe3d3s':
            self.fe_type = 'fe_2d_3_s'
        elif type == 'fe3d4s':
            self.fe_type = 'fe_2d_4_s'
        elif type == 'fe3d6s':
            self.fe_type = 'fe_2d_6_s'
        else:
            print_error('Wrong file format: ' + self.file_name)
            return False
        # Кол-во узлов
        n = int(lines[3])
        # Считываем узлы
        index = 4
        for i in range(n):
            row = list()
            coord = lines[4 + i].split()
            for j in range(len(coord)):
                row.append(float(coord[j]))
            self.x.append(row)
            index += 1
        # Кол-во КЭ
        n = int(lines[index])
        index += 1
        # Считываем КЭ
        for i in range(n):
            row = []
            fe = lines[index].split()
            for j in range(len(fe)):
                row.append(int(fe[j]))
            self.fe.append(row)
            index += 1
        # Кол-во ГЭ
        n = int(lines[index])
        index += 1
        # Считываем ГЭ
        for i in range(n):
            row = []
            be = lines[index].split()
            for j in range(len(be)):
                row.append(int(be[j]))
            self.be.append(row)
            index += 1
        if self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_2d_3_s' or \
                self.fe_type == 'fe_2d_4_p' or self.fe_type == 'fe_2d_6_s' or self.fe_type == 'fe_2d_4_s':
            self.be = self.fe
        # Пропускаем время расчета и заголовок
        index += 2
        # Считываем кол-во функций
        n = int(lines[index])
        index += 1
        for i in range(n):
            res = TResult()
            res.name = lines[index].split('\n')[0]
            res.t = float(lines[index + 1])
            m = int(lines[index + 2])
            index += 3
            for j in range(m):
                res.results.append(float(lines[index]))
                index += 1
            self.results.append(res)
        return True

    def __init_main_menu(self):
        # Пункты главного меню
        file_menu = self.menuBar().addMenu('&File')
        analyse_menu = self.menuBar().addMenu('&Analyse')
        options_menu = self.menuBar().addMenu('&Option')
        help_menu = self.menuBar().addMenu('&?')
        # Настройка File
        open_action = QAction('&Open...', self)
        open_action.setStatusTip('Open a data file')
        open_action.triggered.connect(self.__open_action)
        file_menu.addAction(open_action)
        close_action = QAction('&Close', self)
        close_action.setStatusTip('Close current file')
        close_action.triggered.connect(self.__close_action)
        file_menu.addAction(close_action)
        file_menu.addSeparator()
        exit_action = QAction('E&xit', self)
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Настройка Analyse
        function_action = QAction('&Function...', self)
        function_action.setStatusTip('Choose function to analyse')
        function_action.triggered.connect(self.__function_action)
        analyse_menu.addAction(function_action)

        # Настройка Options
        light_action = QAction('&Light', self)
        light_action.setStatusTip('Enable light')
        light_action.triggered.connect(self.__light_action)
        light_action.setCheckable(True)
        light_action.setChecked(True)
        options_menu.addAction(light_action)
        fe_border_action = QAction('&FE border', self)
        fe_border_action.setStatusTip('Enable drawing FE border')
        fe_border_action.triggered.connect(self.__fe_border_action)
        fe_border_action.setCheckable(True)
        fe_border_action.setChecked(False)
        options_menu.addAction(fe_border_action)
        legend_action = QAction('Le&gend', self)
        legend_action.setStatusTip('Enable drawing legend')
        legend_action.triggered.connect(self.__legend_action)
        legend_action.setCheckable(True)
        legend_action.setChecked(True)
        options_menu.addAction(legend_action)
        options_menu.addSeparator()
        invert_normal_action = QAction('&Invert normal', self)
        invert_normal_action.setStatusTip('Invert polygon normnal vector')
        invert_normal_action.triggered.connect(self.__invert_normal_action)
        invert_normal_action.setCheckable(True)
        invert_normal_action.setChecked(False)
        options_menu.addAction(invert_normal_action)

        light_two_side_action = QAction('&Two-sided lighting', self)
        light_two_side_action.setStatusTip('Invert polygon normnal vector')
        light_two_side_action.triggered.connect(self.__light_two_side_action)
        light_two_side_action.setCheckable(True)
        light_two_side_action.setChecked(True)
        options_menu.addAction(light_two_side_action)

        options_menu.addSeparator()
        qa = QActionGroup(self)
        color_16_action = QAction('&16', self)
        color_16_action.setStatusTip('Set 16 grading colors')
        color_16_action.setActionGroup(qa)
        color_16_action.setCheckable(True)
        color_16_action.setChecked(True)
        color_16_action.triggered.connect(self.__colors_action)
        color_32_action = QAction('3&2', self)
        color_32_action.setStatusTip('Set 32 grading colors')
        color_32_action.setActionGroup(qa)
        color_32_action.setCheckable(True)
        color_32_action.setChecked(False)
        color_32_action.triggered.connect(self.__colors_action)
        color_64_action = QAction('6&4', self)
        color_64_action.setStatusTip('Set 64 grading colors')
        color_64_action.setActionGroup(qa)
        color_64_action.setCheckable(True)
        color_64_action.setChecked(False)
        color_64_action.triggered.connect(self.__colors_action)
        color_128_action = QAction('1&28', self)
        color_128_action.setStatusTip('Set 128 grading colors')
        color_128_action.setActionGroup(qa)
        color_128_action.setCheckable(True)
        color_128_action.setChecked(False)
        color_128_action.triggered.connect(self.__colors_action)
        color_256_action = QAction('25&6', self)
        color_256_action.setStatusTip('Set 256 grading colors')
        color_256_action.setActionGroup(qa)
        color_256_action.setCheckable(True)
        color_256_action.setChecked(False)
        color_256_action.triggered.connect(self.__colors_action)

        color_menu = QMenu('&Colors', self)
        color_menu.setStatusTip('Set grading colors')
        color_menu.addAction(color_16_action)
        color_menu.addAction(color_32_action)
        color_menu.addAction(color_64_action)
        color_menu.addAction(color_128_action)
        color_menu.addAction(color_256_action)
        options_menu.addMenu(color_menu)

        # Настройка Help
        about_action = QAction('&About...', self)
        about_action.setStatusTip('Show information about Plot3d')
        about_action.triggered.connect(self.__about_action)
        help_menu.addAction(about_action)

        self.show()

    def __function_action(self):
        dialog = TFunctionListDialog(self.results, self.current_index)
        if dialog.exec() == QDialog.Accepted:
            self.current_index = dialog.current_index
            self.__gl_widget.set_fun_index(self.current_index)

    def __open_action(self):
        dlg = QFileDialog(self, 'Open data file', '', 'JSON data files (*.json);;QFEM results file (*.qres);;FEM Solver results file (*.res)')
        if dlg.exec_():
            if len(self.x):
                self.__close_action()
            self.__set_file(dlg.selectedFiles()[0])

    def __close_action(self):
        self.setWindowTitle('PyFEM results viewer')
        self.file_name = ''
        self.fe_type = ''
        self.x.clear()
        self.fe.clear()
        self.be.clear()
        self.results.clear()
        self.__gl_widget.clear()
        # Сброс меню
        self.menuBar().actions()[0].menu().actions()[1].setEnabled(False)   # Close
        # self.menuBar().actions()[1].menu().clear()                          # Очистка меню function
        self.menuBar().actions()[1].setEnabled(False)                       # Function
        self.menuBar().actions()[2].setEnabled(False)                       # Options
        self.menuBar().actions()[2].menu().actions()[0].setChecked(True)
        self.menuBar().actions()[2].menu().actions()[1].setChecked(False)
        self.menuBar().actions()[2].menu().actions()[2].setChecked(True)
        self.menuBar().actions()[2].menu().actions()[4].setChecked(False)
        self.menuBar().actions()[2].menu().actions()[5].setChecked(False)
        self.menuBar().actions()[2].menu().actions()[7].menu().actions()[0].setChecked(True)

    def __colors_action(self):
        num = int(QObject.sender(self).text().replace('&', ''))
        self.__gl_widget.set_colors(num)

    def __fun_action(self):
        index = self.__get_fun_index(QObject.sender(self).text()[3:])
        self.__gl_widget.set_fun_index(index)

    def __light_action(self):
        self.menuBar().actions()[2].menu().actions()[4].setEnabled(not self.menuBar().actions()[2].
                                                                   menu().actions()[4].isEnabled())
        self.menuBar().actions()[2].menu().actions()[5].setEnabled(not self.menuBar().actions()[2].
                                                                   menu().actions()[5].isEnabled())
        self.__gl_widget.trigger_light()
        self.repaint()

    def __fe_border_action(self):
        self.__gl_widget.trigger_fe_border()
        self.repaint()

    def __legend_action(self):
        self.__gl_widget.trigger_legend()
        self.repaint()

    def __invert_normal_action(self):
        self.__gl_widget.trigger_invert_normal()
        self.repaint()

    def __light_two_side_action(self):
        self.__gl_widget.trigger_light_two_side()
        self.repaint()

    # Задание нового файла
    def __set_file(self, file_name):
        self.file_name = file_name
        # Загрузка данных
        if self.__load_file() is False:
            QMessageBox.critical(self, 'Error', 'Error read data from file ' + file_name, QMessageBox.Ok)
            return
        # self.__init_function_menu()
        self.__gl_widget.set_data(self.fe_type, self.x, self.fe, self.be, self.results, 0)
        # Настройка меню
        self.menuBar().actions()[0].menu().actions()[1].setEnabled(True)
        self.menuBar().actions()[1].setEnabled(True)                       # Function
        self.menuBar().actions()[2].setEnabled(True)                       # Options

    @staticmethod
    def __about_action(self):
        dialog = QDialog()
        dialog.setWindowTitle('About Plot3d')
        dialog.resize(250, 100)
        dialog.setMaximumSize(250, 100)
        dialog.setMinimumSize(250, 100)
        # Настройка кнопок
        dlg_btn_box= QDialogButtonBox()
        dlg_btn_box.setStandardButtons(QDialogButtonBox.Ok)
        dlg_btn_box.accepted.connect(dialog.close)

        main_layout = QVBoxLayout()
        main_layout.addWidget(QLabel('Plot3d - simple visualizer for PyFEM'))
        main_layout.addWidget(dlg_btn_box)
        dialog.setLayout(main_layout)
        dialog.exec()


# Базовый класс, реализующий основной функционал OpenGL
class TGLWidget(QWidget):
    def __init__(self):
        super(TGLWidget, self).__init__()
        self.fe_type = ''
        self.x = []
        self.dx = []
        self.fe = []
        self.be = []
        self.results = []
        self.min_x = []
        self.max_x = []
        self.x_c = []
        self.fun_index = 0
        self.radius = 0
        self.min_u = 0
        self.max_u = 0
        self.is_light = True
        self.is_legend = True
        self.is_fe_border = False
        self.is_invert_normal = False
        self.is_light_two_side = True
        self.transform_coeff = 0
        self.angle_x = 0
        self.angle_y = 0
        self.angle_z = 0
        self.scale = 1
        self.num_color = 16
        self.__is_idle = True
        self.__last_pos = QPoint()
        self.__color_table = []
        self.__gl = QGLWidget(self)
        self.__gl.initializeGL()
        self.__gl.resizeGL = self.__resize
        self.__gl.paintGL = self.__paint
        self.__init_color_table()
        QVBoxLayout(self).addWidget(self.__gl)
        self.mousePressEvent = self.__mouse_press_event
        self.mouseReleaseEvent = self.__mouse_release_event
        self.__gl.mouseMoveEvent = self.__mouse_move
        self.wheelEvent = self.__wheel_event
        self.__xlist_object = 0
        self.__xlist_skeleton = 0

    def clear(self):
        self.is_light = True
        self.is_legend = True
        self.is_fe_border = False
        self.is_invert_normal = False
        self.is_light_two_side = False
        self.transform_coeff = 0
        self.fe_type = ''
        self.x.clear()
        self.fe.clear()
        self.be.clear()
        self.results.clear()
        self.min_x.clear()
        self.max_x.clear()
        self.x_c.clear()
        self.__color_table.clear()
        if self.__xlist_object != 0:
            glDeleteLists(self.__xlist_object, 1)
            self.__xlist_object = 0
        if self.__xlist_skeleton != 0:
            glDeleteLists(self.__xlist_skeleton, 1)
            self.__xlist_skeleton = 0
        self.__gl.updateGL()

    def set_data(self, fe_type, px, fe, be, results, fun_index):
        self.__gl.glInit()
        self.fe_type = fe_type
        self.x = px
        self.fe = fe
        self.be = be
        self.results = results
        self.fun_index = fun_index
        self.min_x, self.max_x, self.x_c, self.radius = self.__get_coord_info()
        self.min_u = min(self.results[self.fun_index].results)
        self.max_u = max(self.results[self.fun_index].results)
        self.is_light = True
        self.is_legend = True
        self.is_invert_normal = False
        self.is_light_two_side = True
        self.is_fe_border = False
        self.transform_coeff = 0
        self.angle_x = 0
        self.angle_y = 0
        self.angle_z = 0
        self.scale = 1
        self.num_color = 32
        self.__is_idle = True
        self.__last_pos = QPoint()
        self.__color_table = []
        self.__init_color_table()
        if self.__xlist_object != 0:
            glDeleteLists(self.__xlist_object, 1)
            self.__xlist_object = 0
        if self.__xlist_skeleton != 0:
            glDeleteLists(self.__xlist_skeleton, 1)
            self.__xlist_skeleton = 0
        self.__resize(self.width(), self.height())
        self.__gl.updateGL()

    def redraw(self):
        if self.__xlist_object != 0:
            glDeleteLists(self.__xlist_object, 1)
            self.__xlist_object = 0
        if self.__xlist_skeleton != 0:
            glDeleteLists(self.__xlist_skeleton, 1)
            self.__xlist_skeleton = 0
        self.__gl.updateGL()

    def trigger_light(self):
        self.is_light = not self.is_light
        self.redraw()
        self.__gl.update()

    def trigger_fe_border(self):
        self.is_fe_border = not self.is_fe_border
        self.redraw()
        self.__gl.update()

    def trigger_legend(self):
        self.is_legend = not self.is_legend
        self.redraw()
        self.__gl.update()

    def trigger_invert_normal(self):
        self.is_invert_normal = not self.is_invert_normal
        self.redraw()
        self.__gl.update()

    def trigger_light_two_side(self):
        self.is_light_two_side = not self.is_light_two_side
        self.__setup_light()
        self.redraw()
        self.__gl.update()

    def set_colors(self, colors):
        self.num_color = 2 * colors
        self.__init_color_table()
        self.redraw()
        self.__gl.update()

    def set_fun_index(self, fun_index):
        self.fun_index = fun_index
        self.min_u = min(self.results[self.fun_index].results)
        self.max_u = max(self.results[self.fun_index].results)
        self.__init_color_table()
        self.redraw()
        self.__gl.update()

    def __mouse_press_event(self, event):
        super(TGLWidget, self).mousePressEvent(event)
        if event.buttons() & Qt.LeftButton:
            self.__is_idle = False

    def __mouse_release_event(self, event):
        super(TGLWidget, self).mouseReleaseEvent(event)
        if self.__is_idle is False:
            self.__is_idle = True
            self.__gl.repaint()

    def __mouse_move(self, event):
        dx = event.x() - self.__last_pos.x()
        dy = event.y() - self.__last_pos.y()
        self.__last_pos = event.pos()
        if event.buttons() & Qt.LeftButton:
            if event.modifiers() & Qt.ShiftModifier:
                self.angle_z += (dy / abs(dy) if dy != 0 else 0) + (dx / abs(dx) if dx != 0 else 0)
            else:
                self.angle_x += dy / abs(dy) if dy != 0 else 0
                self.angle_y += dx / abs(dx) if dx != 0 else 0
            self.__gl.repaint()

    def __wheel_event(self, event):
        if event.angleDelta().y() > 0:
            self.scale *= 1.05
        else:
            self.scale /= 1.05
        self.__gl.repaint()

    def __get_coord_info(self):
        min_x = [0, 0, 0]
        max_x = [0, 0, 0]
        for k in range(0, len(self.x[0])):
            min_x[k] = min(self.x[i][k] for i in range(0, len(self.x)))
            max_x[k] = max(self.x[i][k] for i in range(0, len(self.x)))
        x_c = [(max_x[0] + min_x[0]) / 2, (max_x[1] + min_x[1]) / 2, (max_x[2] + min_x[2]) / 2]
        radius = ((max_x[0] - min_x[0]) ** 2 + (max_x[1] - min_x[1]) ** 2 + (max_x[2] - min_x[2]) ** 2) ** 0.5
        return min_x, max_x, x_c, radius

    def __init_color_table(self):
        self.__color_table.clear()
        step = self.num_color / 6
        h = 1.0 / step
        green = 0
        blue = 1
        red = 0.24  # Темно-фиолетовый
        u = self.min_u
        h_u = (self.max_u - self.min_u) / float(self.num_color)

        for i in range(0, self.num_color):
            if i < step:
                # фиолетовый-синий
                self.__color_table.append([red, 0, 1, u])
                red -= h
                if red < 0:
                    red = 0
            elif step <= i < 2 * step:
                # синий-голубой
                self.__color_table.append([0, green, 1, u])
                green += h
                if green > 1:
                    green = 1
            elif 2 * step <= i < 3 * step:
                # голубой-зеленый
                self.__color_table.append([0, 1, blue, u])
                blue -= h
                if blue < 0:
                    blue = 0
            elif 3 * step <= i < 4 * step:
                # зеленый-желтый
                self.__color_table.append([red, 1, 0, u])
                red += h
                if red > 1:
                    red = 1
            elif i > 4 * step:
                # желтый-оранжевый-красный
                self.__color_table.append([1, green, 0, u])
                green -= 0.5 * h
                if green < 0:
                    green = 0
            u += h_u

    def color(self, i):
        if self.is_light:
            self.__make_material(self.__color_table[i][0], self.__color_table[i][1], self.__color_table[i][2])
        else:
            glColor3f(self.__color_table[i][0], self.__color_table[i][1], self.__color_table[i][2])

    def __resize(self, w, h):
        aspect = w / h
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(60.0, aspect, 0.01 * self.radius, 10 * self.radius)
        gluLookAt(0, 0, self.radius, 0, 0, 0, 0, 1, 0)
        glViewport(0, 0, w, h)
        glMatrixMode(GL_MODELVIEW)
        if self.is_light:
            self.__setup_light()

    @staticmethod
    def __make_material(r, g, b, a=1.0):
        diffuse = 0.5
        ambient = 0.4
        specular = 0.7
        shininess = 50
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, [r * diffuse, g * diffuse, b * diffuse, a])
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [r * ambient, g * ambient, b * ambient, a])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [specular, specular, specular, a])
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess)

    def __get_color_index(self, u):
        ret = 1
        if self.min_u != self.max_u:
            ret = int(floor((u - self.min_u) / ((self.max_u - self.min_u) / self.num_color))) - 1
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
        step = (self.max_u - self.min_u) / n
        for k in range(0, n):
            if k == n - 1:
                v = stop
            i = self.__get_color_index(v)
            glColor3f(self.__color_table[i][0], self.__color_table[i][1], self.__color_table[i][2])
            font.setStyleStrategy(QFont.OpenGLCompatible)
            self.__gl.renderText(self.rect().width() - font_w1 - font_w2 - 50, cy, '█', font)
            glColor3f(1, 1, 1)
            self.__gl.renderText(self.rect().width() - font_w2 - 50, cy, '{:+0.5E}'.format(v), font)
            cy += font_h
            v -= step
        self.__show_header()

    def __show_header(self):
        fun_name = self.results[self.fun_index].name
        if self.results[len(self.results) - 1].t != 0:
            fun_name += ('(' + str(self.results[self.fun_index].t) + ')')
        font = QFont('Times', 14, QFont.Normal)
        fm = QFontMetrics(font)
        w = fm.width(fun_name)
        h = fm.height()
        self.__gl.renderText(self.rect().width() // 2 - w // 2, h, fun_name, font)

    @staticmethod
    def __sort(tri):
        inverse = 1.0
        min_i = 0
        min_u = tri[0][3]
        for i in range(1, 3):
            if tri[i][3] < min_u:
                min_i = i
                min_u = tri[i][3]
        for i in range(0, min_i):
            tri.append(tri.pop(0))
        if tri[1][3] > tri[2][3]:
            inverse = -1.0
            tri.append(tri.pop(1))
        return inverse

    def draw_triangle_3d(self, tri):
        # tri = sorted(tri, key=lambda item: item[3])
        inv = self.__sort(tri)
        glFrontFace(GL_CCW if inv == 1 else GL_CW)
        if self.is_light:
            # Задание нормали
            a = (tri[1][1] - tri[0][1]) * (tri[2][2] - tri[0][2]) - (tri[2][1] - tri[0][1]) * (tri[1][2] - tri[0][2])
            b = (tri[2][0] - tri[0][0]) * (tri[1][2] - tri[0][2]) - (tri[1][0] - tri[0][0]) * (tri[2][2] - tri[0][2])
            c = (tri[1][0] - tri[0][0]) * (tri[2][1] - tri[0][1]) - (tri[2][0] - tri[0][0]) * (tri[1][1] - tri[0][1])
            if self.is_invert_normal:
                inv *= -1
            glNormal3f(inv * a, inv * b, inv * c)
        color_index = []
        for i in range(0, 3):
            color_index.append(self.__get_color_index(tri[i][3]))
        if color_index[0] == color_index[1] == color_index[2]:
            # Треугольник одного цвета
            self.color(color_index[0])
            glBegin(GL_TRIANGLES)
            for i in range(0, 3):
                glVertex3f(tri[i][0] - self.x_c[0], tri[i][1] - self.x_c[1], tri[i][2] - self.x_c[2])
            glEnd()
        else:
            # Изолинии проходят по треугольнику
            step = color_index[2] - color_index[0] + 1
            p02 = []
            x = [tri[0][0], tri[0][1], tri[0][2], color_index[0]]
            h = [(tri[2][0] - tri[0][0]) / step, (tri[2][1] - tri[0][1]) / step, (tri[2][2] - tri[0][2]) / step,
                 (color_index[2] - color_index[0]) / step]
            for i in range(0, step):
                p02.append([x[0] + i * h[0], x[1] + i * h[1], x[2] + i * h[2], color_index[0] + i * h[3]])
            p02.append([tri[2][0], tri[2][1], tri[2][2], color_index[2]])

            step = color_index[1] - color_index[0] + 1
            p012 = []
            x = [tri[0][0], tri[0][1], tri[0][2], color_index[0]]
            h = [(tri[1][0] - tri[0][0]) / step, (tri[1][1] - tri[0][1]) / step, (tri[1][2] - tri[0][2]) / step,
                 (color_index[1] - color_index[0]) / step]
            for i in range(1, step):
                p012.append([x[0] + i * h[0], x[1] + i * h[1], x[2] + i * h[2], color_index[0] + i * h[3]])
            p012.append([tri[1][0], tri[1][1], tri[1][2], color_index[1]])

            step = color_index[2] - color_index[1] + 1
            x = [tri[1][0], tri[1][1], tri[1][2], color_index[1]]
            h = [(tri[2][0] - tri[1][0]) / step, (tri[2][1] - tri[1][1]) / step, (tri[2][2] - tri[1][2]) / step,
                 (color_index[2] - color_index[1]) / step]
            for i in range(1, step):
                p012.append([x[0] + i * h[0], x[1] + i * h[1], x[2] + i * h[2], color_index[1] + i * h[3]])

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

    def __setup_light(self):
        diffuse = 0.8
        ambient = 0.8
        specular = 0.6

        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE if self.is_light_two_side else GL_FALSE)
#        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
        glLightfv(GL_LIGHT0, GL_AMBIENT, [ambient, ambient, ambient])
        glLightfv(GL_LIGHT0, GL_DIFFUSE, [diffuse, diffuse, diffuse])
        glLightfv(GL_LIGHT0, GL_SPECULAR, [specular, specular, specular])

        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_NORMALIZE)

    def draw_fe_border(self, tri):
        if self.is_light:
            self.__make_material(0, 0, 0)
        else:
            glColor3f(0, 0, 0)
        glBegin(GL_LINE_LOOP)
        n = len(tri) if len(tri) <= 4 else 3
        for i in range(0, n):
            glVertex3f(tri[i][0] - self.x_c[0], tri[i][1] - self.x_c[1], tri[i][2] - self.x_c[2])
        glEnd()

    def __paint(self):
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

        if self.__is_idle:
            self.__display_object()
        else:
            self.__display_skeleton()

        glPopMatrix()
        # Изображение цветовой шкалы
        if self.is_legend and self.__is_idle:
            self.show_legend()

    # Визуализация результата
    def __display_object(self):
        if self.__xlist_object == 0:
            self.__xlist_object = glGenLists(1)
            glNewList(self.__xlist_object, GL_COMPILE)
            if self.fe_type == 'fe_1d_2':
                self.__paint_1d()
            elif self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_4' or self.fe_type == 'fe_2d_3_p' or \
                    self.fe_type == 'fe_2d_4_p' or self.fe_type == 'fe_2d_6' or self.fe_type == 'fe_2d_6_p':
                self.__paint_2d()
            else:
                self.__paint_3d()
            glEndList()
        else:
            glCallList(self.__xlist_object)

    # Изображение каркаса объекта
    def __display_skeleton(self):
        if self.__xlist_skeleton == 0:
            self.__xlist_skeleton = glGenLists(2)
            glNewList(self.__xlist_skeleton, GL_COMPILE)
            if self.fe_type == 'fe_1d_2':
                glBegin(GL_LINES)
                glVertex2f(self.min_x[0] - self.x_c[0], 0)
                glVertex2f(self.max_x[0] - self.x_c[0], 0)
                glEnd()
            elif self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_6' or self.fe_type == 'fe_2d_4':
                for i in range(0, len(self.be)):
                    glBegin(GL_LINES)
                    for j in range(0, len(self.be[0])):
                        glVertex2f(self.x[self.be[i][j]][0] - self.x_c[0], self.x[self.be[i][j]][1] - self.x_c[1])
                    glEnd()
            elif self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_4_p' or self.fe_type == 'fe_2d_6_p':
                for i in range(len(self.fe)):
                    size = len(self.fe[i]) if self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_4_p' else 3
                    glBegin(GL_LINE_LOOP)
                    for j in range(size):
                        glVertex2f(self.x[self.fe[i][j]][0] - self.x_c[0], self.x[self.fe[i][j]][1] - self.x_c[1])
                    glEnd()
            else:
                for i in range(0, len(self.be)):
                    glBegin(GL_LINE_LOOP)
                    size = 4 if self.fe_type == 'fe_3d_10' else 3 if self.fe_type == 'fe_2d_6_s' else len(self.be[0])
                    for j in range(size):
                        glVertex3f(self.x[self.be[i][j]][0] - self.x_c[0], self.x[self.be[i][j]][1] - self.x_c[1],
                                   self.x[self.be[i][j]][2] - self.x_c[2])
                    glEnd()
            glEndList()
        else:
            glCallList(self.__xlist_skeleton)

    # Визуализация одномерной задачи
    def __paint_1d(self):
        # Изображение КЭ
        for i in range(0, len(self.fe)):
            data = []
            for j in range(0, len(self.fe[0])):
                data.append([self.x[self.fe[i][j]][0] + self.transform_coeff * self.results[0].results[self.fe[i][j]],
                            self.results[self.fun_index].results[self.fe[i][j]]])
            color1 = self.__get_color_index(data[0][1])
            color2 = self.__get_color_index(data[1][1])
            delta_color = 1 if color1 < color2 else -1
            step = abs(color1 - color2)
            h = (data[1][0] - data[0][0]) / step if step > 0 else 0
            glBegin(GL_LINE_STRIP)
            color = color1
            self.color(color)
            x = data[0][0]
            glVertex2f(x - self.x_c[0], 0)
            for j in range(0, step - 1):
                x += h
                color += delta_color
                self.color(color)
                glVertex2f(x - self.x_c[0], 0)
            self.color(color2)
            glVertex2f(data[1][0] - self.x_c[0], 0)
            glEnd()

    # Визуализация плоской задачи
    def __paint_2d(self):
        # self.transform_coeff = 1.0E+2
        # Изображение КЭ
        for i in range(0, len(self.fe)):
            # for i in range(212, 214):
            tri = []
            for j in range(0, len(self.fe[0])):
                dx = self.transform_coeff * self.results[0].results[self.fe[i][j]] if self.fe_type != 'fe_2d_3_p' and \
                                                                                      self.fe_type != 'fe_2d_6_p' and \
                                                                                      self.fe_type != 'fe_2d_4_p' else 0
                dy = self.transform_coeff * self.results[1].results[self.fe[i][j]] if self.fe_type != 'fe_2d_3_p' and \
                                                                                      self.fe_type != 'fe_2d_6_p' and \
                                                                                      self.fe_type != 'fe_2d_4_p' else 0
                dz = 0 if self.fe_type != 'fe_2d_3_p' and self.fe_type != 'fe_2d_6_p' and self.fe_type != 'fe_2d_4_p' \
                    else self.transform_coeff * self.results[2].results[self.fe[i][j]]
                tri.append([self.x[self.fe[i][j]][0] + dx, self.x[self.fe[i][j]][1] + dy, dz,
                            self.results[self.fun_index].results[self.fe[i][j]]])
            if len(tri) == 3:
                # Линейный треугольник
                self.draw_triangle_3d(tri)
            elif len(tri) == 4:
                # Четырехугольник
                self.draw_triangle_3d([tri[0], tri[1], tri[2]])
                self.draw_triangle_3d([tri[0], tri[2], tri[3]])
            else:
                # Квадратичный треугольник
                self.draw_triangle_3d([tri[0], tri[3], tri[5]])
                self.draw_triangle_3d([tri[3], tri[1], tri[4]])
                self.draw_triangle_3d([tri[3], tri[4], tri[5]])
                self.draw_triangle_3d([tri[5], tri[4], tri[2]])
            if self.is_fe_border:
                self.draw_fe_border(tri)

    # Визуализация пространственной задачи
    def __paint_3d(self):
        # Изображение поверхности
        # self.transform_coeff = 1.0E+5
        for i in range(0, len(self.be)):
            tri = []
            for j in range(len(self.be[0])):
                tri.append([self.x[self.be[i][j]][0] + self.transform_coeff * self.results[0].results[self.be[i][j]],
                            self.x[self.be[i][j]][1] + self.transform_coeff * self.results[1].results[self.be[i][j]],
                            self.x[self.be[i][j]][2] + self.transform_coeff * self.results[2].results[self.be[i][j]],
                            self.results[self.fun_index].results[self.be[i][j]]])
            if len(tri) == 3:
                self.draw_triangle_3d(tri)
            elif len(tri) == 4:
                self.draw_triangle_3d([tri[0], tri[1], tri[2]])
                self.draw_triangle_3d([tri[0], tri[2], tri[3]])
            elif len(tri) == 6:
                self.draw_triangle_3d([tri[0], tri[3], tri[5]])
                self.draw_triangle_3d([tri[3], tri[4], tri[5]])
                self.draw_triangle_3d([tri[3], tri[1], tri[4]])
                self.draw_triangle_3d([tri[5], tri[4], tri[2]])

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
