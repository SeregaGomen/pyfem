#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
from core.fem_defs import DIR_1, DIR_2, DIR_3, INIT_U, INIT_V, INIT_W, INIT_U_T, INIT_V_T, INIT_W_T, INIT_U_T_T, \
    INIT_V_T_T, INIT_W_T_T
from core.fem_object import TObject
from plot.plot3d import TPlot


def body1d(res_name):
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/body1d.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'x=0', DIR_1)
    #    obj.add_volume_load('-1.0E+5', '', DIR_X)
        obj.add_concentrated_load('-1.0E+5', 'x=1', DIR_1)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def cube(res_name):
    obj = TObject()
    e = [203200]
    m = [0.27]
    if obj.set_mesh('mesh/cube.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'z=0', DIR_1 | DIR_2 | DIR_3)
    #    obj.add_volume_load('-1000', '', DIR_Z)
        obj.add_surface_load('-1000', 'z=1', DIR_3)
    #    obj.add_concentrated_load('-1000', 'z=1', DIR_Z)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def cube_test(res_name):
    obj = TObject()
    e = [203200]
    m = [0.27]
    if obj.set_mesh('mesh/cube_test.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'y=0', DIR_1 | DIR_2)
        obj.add_surface_load('-1000', 'y=1', DIR_2)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def beam(res_name):
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/beam.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('iterative')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'y=0', DIR_1 | DIR_2 | DIR_3)
        obj.add_volume_load('-1.0E+5', '', DIR_2)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def beam_dynamic(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/beam.trpa'):
        obj.set_problem_type('dynamic')
        obj.set_solve_method('iterative')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity([6.5E+10], [0.3])
        obj.set_density(1.0E+3)
        obj.set_time(0, 1.0, 0.25)
        obj.add_boundary_condition('0', 'y=0', DIR_1 | DIR_2 | DIR_3)
        obj.add_volume_load('-1.0E+5*cos(t)', '', DIR_2)
        obj.add_initial_condition('0', INIT_U)
        obj.add_initial_condition('0', INIT_V)
        obj.add_initial_condition('0', INIT_W)
        obj.add_initial_condition('0', INIT_U_T)
        obj.add_initial_condition('0', INIT_V_T)
        obj.add_initial_condition('0', INIT_W_T)
        obj.add_initial_condition('0', INIT_U_T_T)
        obj.add_initial_condition('0', INIT_V_T_T)
        obj.add_initial_condition('0', INIT_W_T_T)
        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def console(res_name):
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/console.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
    #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'x=0', DIR_1 | DIR_2)
        obj.add_concentrated_load('-1.0E+6', 'x=10', DIR_2)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def console4(res_name):
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/console4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
    #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'x=0', DIR_1 | DIR_2)
    #    obj.add_concentrated_load('-1.0E+5', 'x=10', DIR_Y)
        obj.add_volume_load('-1.0E+5', '', DIR_2)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def quad(res_name):
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/quad.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'y=0', DIR_1 | DIR_2)
        obj.add_concentrated_load('-1.0E+5', 'y=1', DIR_2)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def cylinder(res_name):
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/cyl.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
    #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_variable('eps', 1.0E-6)
        obj.add_boundary_condition('0', 'x=0', DIR_1 | DIR_2 | DIR_3)
        obj.add_boundary_condition('0', 'x=2', DIR_1 | DIR_2 | DIR_3)
        obj.add_concentrated_load('-1.0e+4*cos(atan2(y,z))', 'abs(y^2 + z^2 - 0.5^2) <= eps', DIR_3)
        obj.add_concentrated_load('-1.0e+4*sin(atan2(y,z))', 'abs(y^2 + z^2 - 0.5^2) <= eps', DIR_2)
        obj.add_concentrated_load('1.0e+4*cos(atan2(y,z))', 'abs(y^2 + z^2 - 0.25^2) <= eps', DIR_3)
        obj.add_concentrated_load('1.0e+4*sin(atan2(y,z))', 'abs(y^2 + z^2 - 0.25^2) <= eps', DIR_2)
        if obj.calc():
            # obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def tank3(res_name):
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/tank3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_elasticity(e, m)
        obj.add_variable('eps', 1.0E-6)
        obj.add_variable('min', 0.0015)
        obj.set_width(10)
        obj.set_precision(5)
        obj.add_boundary_condition('0', 'y=-0.598 and abs(x^2+z^2-1.6635^2)<=eps', DIR_1 | DIR_2 | DIR_3)
        obj.add_boundary_condition('0', 'x=0', DIR_1)
        obj.add_boundary_condition('0', 'z=0', DIR_3)
        obj.add_surface_load('1.0e+4*cos(atan2(z,x))',
                             '(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037-min)^2) <= eps)', DIR_1)
        obj.add_surface_load('1.0e+4*sin(atan2(z,x))',
                             '(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037-min)^2) <= eps)', DIR_3)
        obj.add_surface_load('1.0e+4*cos(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037-min)^2) <= eps)', DIR_1)
        obj.add_surface_load('1.0e+4*cos(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037-min)^2) <= eps)', DIR_2)
        obj.add_surface_load('1.0e+4*sin(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037-min)^2) <= eps)', DIR_3)
        obj.add_surface_load('1.0e+4*cos(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,y))',
                             '(y > 0) and (abs(x^2 + y^2 + z^2 - (1.037-min)^2) <= eps)', DIR_1)
        obj.add_surface_load('1.0e+4*cos(atan2((x^2+z^2)^0.5,y))',
                             '(y > 0) and (abs(x^2 + y^2 + z^2 - (1.037-min)^2) <= eps)', DIR_2)
        obj.add_surface_load('1.0e+4*sin(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,y))',
                             '(y > 0) and (abs(x^2 + y^2 + z^2 - (1.037-min)^2) <= eps)', DIR_3)
        obj.add_surface_load('-5.0e+3*cos(atan2(z,x))', '(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - \
        (1.037)^2) <= eps)', DIR_1)
        obj.add_surface_load('-5.0e+3*sin(atan2(z,x))', '(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - \
        (1.037)^2) <= eps)', DIR_3)
        obj.add_surface_load('-5.0e+3*cos(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037)^2) <= eps)', DIR_1)
        obj.add_surface_load('-5.0e+3*cos(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037)^2) <= eps)', DIR_2)
        obj.add_surface_load('-5.0e+3*sin(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037)^2) <= eps)', DIR_3)
        obj.add_surface_load('-5.e+3', '(y=-1.7235) and (x^2+z^2 - 0.34205^2 <= eps)', DIR_2)
        obj.add_surface_load('-5.e+3', '(y=-1.944) and (x^2+z^2 - 0.657857^2 <= eps and x^2+z^2 - 0.562143^2 >= eps)',
                             DIR_2)
        obj.add_surface_load('5.0e+3*cos(atan2(z,x))', 'abs(y+0.6431) <= eps and abs(x^2 + z^2 - 1.6389^2) <= eps',
                             DIR_1)
        obj.add_surface_load('5.0e+3*sin(atan2(z,x))', 'abs(y+0.6431) <= eps and abs(x^2 + z^2 - 1.6389^2) <= eps',
                             DIR_3)
        obj.add_surface_load('5.0e+3*x*(1.0644108554^2)/(((x*(1.0644108554^2))^2+(y+1.1013629509)^2+\
        (z*(1.0644108554^2))^2)^0.5)',
                             '(y>-0.6431 and y <-0.0234) and abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps',
                             DIR_1)
        obj.add_surface_load('5.0e+3*(y+1.1013629509)/(((x*(1.0644108554^2))^2+(y+1.1013629509)^2+\
        (z*(1.0644108554^2))^2)^0.5)', '(y>-0.6431 and y <-0.0234) and abs(y-((x^2+z^2)^0.5)*\
        (-1.0644108554)-1.1013629509)<=eps', DIR_2)
        obj.add_surface_load('5.0e+3*z*(1.0644108554^2)/(((x*(1.0644108554^2))^2+(y+1.1013629509)^2+\
        (z*(1.0644108554^2))^2)^0.5)', '(y>-0.6431 and y <-0.0234) and abs(y-((x^2+z^2)^0.5)*\
        (-1.0644108554)-1.1013629509)<=eps', DIR_3)
        obj.add_surface_load('-5.0e+3*x*(1.0018498686^2)/(((x*(1.0018498686^2))^2+(z*(1.0018498686^2))^2+\
        (y-1.3808172524)^2)^0.5)', '(y>-1.944 and y <-1.7235) and abs(y - ((x^2+z^2)^0.5)*\
        (-1.0018498686)+1.3808172524)<=eps', DIR_1)
        obj.add_surface_load('5.0e+3*(y-1.3808172524)/(((x*(1.0018498686^2))^2+(z*(1.0018498686^2))^2+\
        (y-1.3808172524)^2)^0.5)', '(y>-1.944 and y <-1.7235) and abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+\
        1.3808172524)<=eps', DIR_2)
        obj.add_surface_load('-5.0e+3*z*(1.0018498686^2)/(((x*(1.0018498686^2))^2+(z*(1.0018498686^2))^2+\
        (y-1.3808172524)^2)^0.5)', '(y>-1.944 and y <-1.7235) and abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+\
        1.3808172524)<=eps', DIR_3)
        obj.add_surface_load('5.0e+3*x*(1.3260378897^2)/(((3*x*(1.3260378897^2))^2+(y-2.8163434974)^2+\
        (3*z*(1.3260378897^2))^2)^0.5)', '(y>-1.944 and y < -0.6431) and abs(y-((x^2+z^2)^0.5)*(1.3260378897)+\
        2.8163434974)<=eps', DIR_1)
        obj.add_surface_load('5.0e+3*(y-2.8163434974)/(((3*x*(1.3260378897^2))^2+(y-2.8163434974)^2+\
        (3*z*(1.3260378897^2))^2)^0.5)', '(y>-1.944 and y < -0.6431) and abs(y-((x^2+z^2)^0.5)*(1.3260378897)+\
        2.8163434974)<=eps', DIR_2)
        obj.add_surface_load('5.0e+3*z*(1.3260378897^2)/(((3*x*(1.3260378897^2))^2+(y-2.8163434974)^2+\
        (3*z*(1.3260378897^2))^2)^0.5)', '(y>-1.944 and y < -0.6431) and abs(y-((x^2+z^2)^0.5)*(1.3260378897)+\
        2.8163434974)<=eps', DIR_3)
        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def head3d(res_name):
    obj = TObject()
    e = [1000]
    m = [0.3]
    if obj.set_mesh('mesh/head3d.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
    #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_variable('eps', 1.0E-6)
        obj.add_boundary_condition('0', 'y=0', DIR_1 | DIR_3)
        obj.add_boundary_condition('0', 'y=991.3', DIR_1 | DIR_2 | DIR_3)
        obj.add_surface_load('-1*cos(atan2(z,x))', 'abs(x^2 + z^2 - 210^2) <=0.001', DIR_1)
        obj.add_surface_load('-1*sin(atan2(z,x))', 'abs(x^2 + z^2 - 210^2) <= 0.001', DIR_3)
        if obj.calc():
            obj.save_result('head3d')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def console_dynamic(res_name):
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/console.trpa'):
        obj.set_problem_type('dynamic')
        obj.set_solve_method('direct')
        obj.set_density(1.0E+3)
        obj.set_damping([3.383, 0.00206])
        obj.set_time(0, 1.0, 0.25)
        obj.set_width(10)
        obj.set_precision(5)
    #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'x=0', DIR_1 | DIR_2)
        obj.add_concentrated_load('-1.0E+5*cos(t)', 'x=10', DIR_1)
        obj.add_initial_condition('0', INIT_U)
        obj.add_initial_condition('0', INIT_V)
        obj.add_initial_condition('0', INIT_U_T)
        obj.add_initial_condition('0', INIT_V_T)
        obj.add_initial_condition('0', INIT_U_T_T)
        obj.add_initial_condition('0', INIT_V_T_T)
        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def plate4(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/plate4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_h(0.01)
        obj.set_precision(5)
        obj.set_elasticity([2E+6], [0.3])
#        obj.set_names(['W', 'Tx', 'Ty', 'Exx', 'Eyy', 'Exy', 'Sxx', 'Syy', 'Sxy'])
        obj.add_boundary_condition('0', 'x = -0.1 or x = 0.1 or y = -0.1 or y = 0.1', DIR_1 | DIR_2 | DIR_3)
#        obj.add_boundary_condition('0', 'y = -0.1 or y = 0.1', DIR_1 | DIR_2)
#        obj.add_concentrated_load('-1.0E+5', 'x = 0 and y = 0', DIR_1)
#        obj.add_volume_load('-1.0E+5', '', DIR_1)
        obj.add_surface_load('-2000', '', DIR_1)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def plate3(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/plate3_1.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_h(0.01)
        obj.set_precision(5)
        obj.set_elasticity([2E+6], [0.3])
        obj.add_boundary_condition('0', 'x = -0.1 or x = 0.1 or y = -0.1 or y = 0.1', DIR_1 | DIR_2 | DIR_3)
#        obj.add_boundary_condition('0', 'x = -0.1 or x = 0.1 or y = -0.1 or y = 0.1', DIR_1)
#        obj.add_concentrated_load('-1.0E+5', 'x = 0 and y = 0', DIR_1)
        obj.add_surface_load('-2000', '', DIR_1)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell3(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/shell3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_h(0.01)
        obj.set_precision(5)
        obj.set_elasticity([2E+6], [0.3])
#        obj.add_boundary_condition('0', 'y = 0', DIR_1 | DIR_2 | DIR_3)
#        obj.add_volume_load('-1.0E+2', '', DIR_2)
        obj.add_boundary_condition('0', 'y = 0', DIR_2)
        obj.add_boundary_condition('0', 'z = 0', DIR_1 | DIR_2 | DIR_3)
        obj.add_concentrated_load('-1.0E+2', 'z = 5', DIR_2)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell4(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/shell4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_h(0.01)
        obj.set_precision(5)
        obj.set_elasticity([2E+6], [0.3])
#        obj.set_names(['x1', 'x2',  'x3', 't', 'u1', 'u2', 'u3', 'e11', 'e22', 'e33', 'e12', 'e13', 'e23', 's11',
#                       's22', 's33', 's12', 's13', 's23', 'u1t', 'u2t', 'u3t', 'u1tt', 'u2tt', 'u3tt'])
#        obj.add_boundary_condition('0', 'x1 = -0.1 or x1 = 0.1 or x2 = -0.1 or x2 = 0.1', DIR_1 | DIR_2 | DIR_3)
        obj.add_boundary_condition('0', 'x = -0.1 or x = 0.1 or y = -0.1 or y = 0.1', DIR_1 | DIR_2 | DIR_3)
        obj.add_surface_load('-2000', '', DIR_1)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell_plate3(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/shell_plate3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_h(0.01)
        obj.set_precision(5)
        obj.set_elasticity([2E+6], [0.3])
        obj.add_boundary_condition('0', 'x = -0.1 or x = 0.1 or y = -0.1 or y = 0.1', DIR_1 | DIR_2 | DIR_3)
#        obj.add_concentrated_load('-1.0E+5', 'x = 0 and y = 0', DIR_1)
#        obj.add_surface_load('-2000', '', DIR_1)
        obj.add_concentrated_load('-2000', 'z = 0', DIR_1)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell3_test(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/shell-tube.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_h(0.0369)
        obj.set_elasticity([203200], [0.27])
        obj.add_boundary_condition('0', 'z = 0 or z = 4.014', DIR_1 | DIR_2 | DIR_3)
        obj.add_surface_load('0.05*cos(atan2(y,x))', '(abs(x^2 + y^2 - 1.99^2) <= 1.0E-3)', DIR_1)
        obj.add_surface_load('0.05*sin(atan2(y,x))', '(abs(x^2 + y^2 - 1.99^2) <= 1.0E-3)', DIR_2)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def plate3_test(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/plate3_1_0.trpa'):
#    if obj.set_mesh('mesh/plate3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_h(0.01)
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity([203200], [0.27])
        obj.add_boundary_condition('0', 'x = -0.5 or x = 0.5 or y = -0.5 or y = 0.5', DIR_1 | DIR_2 | DIR_3)
        obj.add_surface_load('0.05', '', DIR_1)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def plate4_test(res_name):
    obj = TObject()
#    if obj.set_mesh('mesh/plate4-1.0.trpa'):
    if obj.set_mesh('mesh/plate4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_h(0.01)
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity([203200], [0.27])
        obj.add_boundary_condition('0', 'x = -0.5 or x = 0.5 or y = -0.5 or y = 0.5', DIR_1 | DIR_2 | DIR_3)
        obj.add_surface_load('0.05', '', DIR_1)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell4_test(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/shell4_1_0.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_h(0.0369)
        obj.set_elasticity([203200], [0.27])
        obj.add_boundary_condition('0', 'z = 0 or z = 4.014', DIR_1 | DIR_2 | DIR_3)
        obj.add_surface_load('0.05*cos(atan2(y,x))', '(abs(x^2 + y^2 - 1.99^2) <= 1.0E-3)', DIR_1)
        obj.add_surface_load('0.05*sin(atan2(y,x))', '(abs(x^2 + y^2 - 1.99^2) <= 1.0E-3)', DIR_2)
#        obj.add_concentrated_load('50000*cos(atan2(y,x))', '(abs(x^2 + y^2 - 1.99^2) <= 1.0E-5)', DIR_1)
#        obj.add_concentrated_load('50000*sin(atan2(y,x))', '(abs(x^2 + y^2 - 1.99^2) <= 1.0E-5)', DIR_2)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def tube_test(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/tube-solid-test.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity([203200000000], [0.27])
        obj.add_boundary_condition('0', 'z = 0 or z = 4.014', DIR_1 | DIR_2 | DIR_3)
        obj.add_surface_load('50000*cos(atan2(y,x))', '(abs(x^2 + y^2 - 1.99^2) <= 1.0E-3)', DIR_1)
        obj.add_surface_load('50000*sin(atan2(y,x))', '(abs(x^2 + y^2 - 1.99^2) <= 1.0E-3)', DIR_2)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def tank3s(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/tank3s.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
#        obj.set_solve_method('iterative')
        obj.set_elasticity([6.5E+10], [0.3])
        obj.set_h(0.0015)
        obj.add_variable('eps', 0.01)
        obj.add_variable('L', 0.269)
        obj.add_variable('P', 10000)
        obj.add_variable('p', 5000)
        obj.add_variable('R', 1.037)
        obj.set_width(10)
        obj.set_precision(5)
        obj.add_boundary_condition('0', 'y = -0.643 and abs(x^2 + z^2 -1.641^2)<=eps', DIR_1 | DIR_2 | DIR_3)
        obj.add_boundary_condition('0', 'abs(x) <= eps', DIR_1)
        obj.add_boundary_condition('0', 'abs(z) <= eps', DIR_3)

        obj.add_surface_load('P*cos(atan2(z,x))', '(y <= 0 and y>=-L) and (abs(x^2 + z^2 - R^2) <= eps)', DIR_1)
        obj.add_surface_load('P*sin(atan2(z,x))', '(y <= 0 and y>=-L) and (abs(x^2 + z^2 - R^2) <= eps)', DIR_3)

        obj.add_surface_load('P*cos(atan2(z,x))*sin(atan2((x^2 + z^2)^0.5,(y + L)))',
                             '(y < -L) and (abs(x^2 + z^2 + (y + L)^2 - R^2) <= eps)', DIR_1)
        obj.add_surface_load('P*cos(atan2((x^2 + z^2)^0.5,(y + L)))',
                             '(y < -L) and (abs(x^2 + z^2 + (y + L)^2 - R^2) <= eps)', DIR_2)
        obj.add_surface_load('P*sin(atan2(z,x))*sin(atan2((x^2 + z^2)^0.5,(y + L)))',
                             '(y < -L) and (abs(x^2 + z^2 + (y + L)^2 - R^2) <= eps)', DIR_3)

        obj.add_surface_load('P*cos(atan2(z,x))*sin(atan2((x^2 + z^2)^0.5,y))',
                             '(y > 0) and (abs(x^2 + y^2 + z^2 - R^2) <= eps)', DIR_1)
        obj.add_surface_load('P*cos(atan2((x^2 + z^2)^0.5,y))',
                             '(y > 0) and (abs(x^2 + y^2 + z^2 - R^2) <= eps)', DIR_2)
        obj.add_surface_load('P*sin(atan2(z,x))*sin(atan2((x^2 + z^2)^0.5,y))',
                             '(y > 0) and (abs(x^2 + y^2 + z^2 - R^2) <= eps)', DIR_3)

        obj.add_surface_load('-p*cos(atan2(z,x))', '(y <= 0 and y>=-L) and (abs(x^2 + z^2 - R^2) <= eps)', DIR_1)
        obj.add_surface_load('-p*sin(atan2(z,x))', '(y <= 0 and y>=-L) and (abs(x^2 + z^2 - R^2) <= eps)', DIR_2)

        obj.add_surface_load('-p*cos(atan2(z,x))*sin(atan2((x^2 + z^2)^0.5,(y + L)))',
                             '(y < -L) and (abs(x^2 + z^2 + (y + L)^2 - R^2) <= eps)', DIR_1)
        obj.add_surface_load('-p*cos(atan2((x^2 + z^2)^0.5,(y + L)))',
                             '(y < -L) and (abs(x^2 + z^2 + (y + L)^2 - R^2) <= eps)', DIR_2)
        obj.add_surface_load('-p*sin(atan2(z,x))*sin(atan2((x^2 + z^2)^0.5,(y + L)))',
                             '(y < -L) and (abs(x^2 + z^2 + (y + L)^2 - R^2) <= eps)', DIR_3)

        obj.add_surface_load('-p', '(y = -1.724) and (x^2+z^2 - 0.342^2 <= eps)', DIR_2)
        obj.add_surface_load('-p', '(y = -1.944) and (x^2+z^2 - 0.660^2 <= eps)', DIR_2)

        obj.add_surface_load('p*cos(atan2(z,x))', 'abs(y + 0.641) <= eps and abs(x^2 + z^2 - 1.643^2) <= eps', DIR_1)
        obj.add_surface_load('p*sin(atan2(z,x))', 'abs(y + 0.641) <= eps and abs(x^2 + z^2 - 1.643^2) <= eps', DIR_2)

        obj.add_surface_load('p*x*(1.0644108554^2)/(((x*(1.0644108554^2))^2+(y + 1.1013629509)^2 + '
                             '(z*(1.0644108554^2))^2)^0.5)', '(y > -0.641 and y <-0.0234) and '
                             'abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps', DIR_1)
        obj.add_surface_load('p*(y+1.1013629509)/(((x*(1.0644108554^2))^2 + (y + 1.1013629509)^2 + '
                             '(z*(1.0644108554^2))^2)^0.5)', '(y > -0.6431 and y < -0.0234) and '
                             'abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps', DIR_2)
        obj.add_surface_load('p*z*(1.0644108554^2)/(((x*(1.0644108554^2))^2 + (y + 1.1013629509)^2 + '
                             '(z*(1.0644108554^2))^2)^0.5)', '(y>-0.6431 and y <-0.0234) and '
                             'abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps', DIR_3)

        obj.add_surface_load('-p*x*(1.0018498686^2)/(((x*(1.0018498686^2))^2 + (z*(1.0018498686^2))^2 + '
                             '(y-1.3808172524)^2)^0.5)', '(y>-1.944 and y <-1.7235) and '
                             'abs(y - ((x^2 + z^2)^0.5)*(-1.0018498686)+1.3808172524)<=eps', DIR_1)
        obj.add_surface_load('p*(y-1.3808172524)/(((x*(1.0018498686^2))^2 + (z*(1.0018498686^2))^2 + '
                             '(y - 1.3808172524)^2)^0.5)', '(y>-1.944 and y <-1.7235) and '
                             'abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+1.3808172524)<=eps', DIR_2)
        obj.add_surface_load('-p*z*(1.0018498686^2)/(((x*(1.0018498686^2))^2 + (z*(1.0018498686^2))^2 + '
                             '(y - 1.3808172524)^2)^0.5)', '(y>-1.944 and y <-1.7235) and '
                             'abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+1.3808172524)<=eps', DIR_3)

        obj.add_surface_load('p*x*(1.3260378897^2)/(((3*x*(1.3260378897^2))^2 + (y - 2.8163434974)^2 + '
                             '(3*z*(1.3260378897^2))^2)^0.5)', '(y>-1.944 and y < -0.6431) and '
                             'abs(y-((x^2+z^2)^0.5)*(1.3260378897)+2.8163434974)<=eps', DIR_1)
        obj.add_surface_load('p*(y-2.8163434974)/(((3*x*(1.3260378897^2))^2 + (y - 2.8163434974)^2 + '
                             '(3*z*(1.3260378897^2))^2)^0.5)', '(y>-1.944 and y < -0.6431) and '
                             'abs(y-((x^2+z^2)^0.5)*(1.3260378897)+2.8163434974)<=eps', DIR_2)
        obj.add_surface_load('p*z*(1.3260378897^2)/(((3*x*(1.3260378897^2))^2+(y - 2.8163434974)^2 + '
                             '(3*z*(1.3260378897^2))^2)^0.5)', '(y>-1.944 and y < -0.6431) and '
                             'abs(y - ((x^2 + z^2)^0.5)*(1.3260378897) + 2.8163434974)<=eps', DIR_3)

        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def plate3d(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/plate3d.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity([203200000000], [0.27])
        obj.add_boundary_condition('0', 'x = -0.5 or x = 0.5 or y = -0.5 or y = 0.5', DIR_1 | DIR_2 | DIR_3)
        obj.add_surface_load('-50000', 'z = 0', DIR_3)
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def create_plate_mesh_4():
    x_min = [-0.5, -0.5]
    x_max = [0.5, 0.5]
    n = 200
    h = [(x_max[0] - x_min[0])/n, (x_max[1] - x_min[1])/n]
    index = []
    x = []
    counter = 0
    for i in range(0, n + 1):
        c_index = []
        for j in range(0, n + 1):
            c_index.append(counter)
            counter += 1
            x.append([x_min[0] + i*h[0], x_min[1] + j*h[1]])
        index.append(c_index)
    with open('mesh/plate4-1.0.trpa', 'w') as file:
        file.write('124\n')
        file.write(str(counter) + '\n')
        for i in range(0, len(x)):
            file.write(str(x[i][0]) + ' ' + str(x[i][1]) + '\n')
        file.write(str(n**2) + '\n')
        for i in range(0, len(index) - 1):
            for j in range(0, len(index) - 1):
                file.write(str(index[i][j]) + ' ' +str(index[i][j + 1]) + ' ' + str(index[i + 1][j + 1]) + ' ' +
                           str(index[i + 1][j]) + '\n')
        file.write('0\n')
    return


def create_shell_mesh_4():
    r = 3.98/2
    height = 4.014
    n_xy = 100
    n_z = 50
    d_fi = 2*math.pi/n_xy
    d_h = height/n_z
    index = []
    x = []
    counter = 0
    for i in range(0, n_z + 1):
        c_index = []
        for j in range(0, n_xy):
            c_index.append(counter)
            counter += 1
            x.append([r*math.cos(j*d_fi), r*math.sin(j*d_fi), i*d_h])
        index.append(c_index)
    with open('mesh/shell4-1.0.trpa', 'w') as file:
        file.write('224\n')
        file.write(str(counter) + '\n')
        for i in range(0, len(x)):
            file.write(str(x[i][0]) + ' ' + str(x[i][1]) + ' ' + str(x[i][2]) + '\n')
        file.write(str(n_xy*n_z) + '\n')
        for i in range(0, len(index) - 1):
            for j in range(0, len(index[0]) - 1):
                file.write(str(index[i][j]) + ' ' +str(index[i][j + 1]) + ' ' + str(index[i + 1][j + 1]) + ' ' +
                           str(index[i + 1][j]) + '\n')
            file.write(str(index[i][j + 1]) + ' ' + str(index[i][0]) + ' ' + str(index[i + 1][0]) + ' ' +
                       str(index[i + 1][j + 1]) + '\n')
        file.write('0\n')
    return


if __name__ == '__main__':
    # create_shell_mesh_4()
    # create_plate_mesh_4()
    # beam('beam')
    # head3d('head3d')
    # cube('cube')
    # tank3('tank3')
    # cylinder('cylinder')
    # quad('quad')
    # cube_test('cube_test')
    # console_dynamic('console_dynamic')
    # console('console')
    # console4('console4')
    # body1d('body1d')
    # plate4('plate4')
    # body1d('body1d')
    # beam_dynamic('beam_dynamic')
    # shell3('shell3')
    # shell4('shell4')
    # shell_plate3('shell_plate3')
    # plate3('plate3')
    # plate3_test('plate3_test')
    # plate4_test('plate4_test')
    # shell4_test('shell4_test')
    # shell3_test('shell3_test')
    # tube_test('tube_test')
    # plate3d('plate3d')
    tank3s('tank3s')



'''
2. Правильно отображать динамическую задачу в plot3d
3. Рисовать оси координат
4. Правильно рисовать маленькие величины (стержень) 
'''
