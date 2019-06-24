#!/usr/bin/env python
# -*- coding: utf-8 -*-

from core.fem_defs import DIR_1, DIR_2, DIR_3, INIT_U, INIT_V, INIT_W, INIT_U_T, INIT_V_T, INIT_W_T, INIT_U_T_T, \
    INIT_V_T_T, INIT_W_T_T
from core.fem_object import TObject
from plot.plot3d import TPlot


def body1d(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/body1d.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0')
        obj.add_boundary_condition(DIR_1, '0', 'x == 0')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        obj.add_volume_load(DIR_1, '-0.5')
        # obj.add_surface_load(DIR_1, '0.5', 'x == 1')
        # obj.add_concentrated_load(DIR_1, '0.5', 'x == 1')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def cube(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/cube.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')

        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z == 0')
        obj.add_volume_load(DIR_3, '-0.5')
        # obj.add_surface_load(DIR_3, '-0.5', 'z == 1')
        # obj.add_concentrated_load(DIR_3, '-1000', 'z == 1 and x == 0 and y == 0')
        # obj.add_concentrated_load(DIR_3, '-1000', 'z == 1 and x == 1 and y == 0')
        # obj.add_concentrated_load(DIR_3, '-1000', 'z == 1 and x == 0 and y == 1')
        # obj.add_concentrated_load(DIR_3, '-1000', 'z == 1 and x == 1 and y == 1')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def cube4(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/cube-4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z == 0')
        # obj.add_volume_load(DIR_3, '-0.5')
        obj.add_surface_load(DIR_3, '-0.5', 'z == 1')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def cube10(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/cube-10.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z == 0')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        obj.add_volume_load(DIR_3, '-0.5')
        # obj.add_surface_load(DIR_3, '-0.5', 'z == 1')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def beam(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/beam.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('iterative')
        obj.add_young_modulus('6.5E+10')
        obj.add_poisson_ratio('0.3')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'y == 0')
        obj.add_volume_load(DIR_2, '-1.0E+5')
        # obj.add_surface_load(DIR_2, '-1.0E+5', 'y=4')
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
        obj.add_young_modulus('6.5E+10')
        obj.add_poisson_ratio('0.3')
        obj.add_density('1.0E+3')
        obj.set_time(0, 1.0, 0.25)
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'y == 0')
        obj.add_volume_load(DIR_2, '-1.0E+5*cos(t)')
        obj.add_initial_condition(INIT_U, '0')
        obj.add_initial_condition(INIT_V, '0')
        obj.add_initial_condition(INIT_W, '0')
        obj.add_initial_condition(INIT_U_T, '0')
        obj.add_initial_condition(INIT_V_T, '0')
        obj.add_initial_condition(INIT_W_T, '0')
        obj.add_initial_condition(INIT_U_T_T, '0')
        obj.add_initial_condition(INIT_V_T_T, '0')
        obj.add_initial_condition(INIT_W_T_T, '0')
        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def console(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/console.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
    #    obj.set_solve_method('iterative')
        obj.add_young_modulus('6.5E+10')
        obj.add_poisson_ratio('0.3')
        obj.add_boundary_condition(DIR_1 | DIR_2, '0', 'x == 0')
        obj.add_concentrated_load(DIR_2, '-1.0E+6', 'x == 10 and y == -0.25')
        # obj.add_surface_load(DIR_2, '-1.0E+6', 'y == 0.25')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def console4(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/console4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
    #    obj.set_solve_method('iterative')
        obj.add_young_modulus('6.5E+10')
        obj.add_poisson_ratio('0.3')
        obj.add_boundary_condition(DIR_1 | DIR_2, '0', 'x == 0')
    #    obj.add_concentrated_load(DIR_2, '-1.0E+5', 'x=10')
        obj.add_volume_load(DIR_2, '-1.0E+5')
        # obj.add_surface_load(DIR_2, '-1.0E+6', 'y == 1')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def rod4(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/rod-4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z = 0')
        obj.add_volume_load(DIR_3, '-0.0765')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def rod10(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/rod-10.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z = 0')
        obj.add_volume_load(DIR_3, '-0.0765')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def quad(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/quad.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('6.5E+10')
        obj.add_poisson_ratio('0.3')
        obj.add_boundary_condition(DIR_1 | DIR_2, '0', 'y == 0')
        obj.add_concentrated_load(DIR_2, '-1.0E+5', 'y == 1')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def cylinder(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/cyl.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
    #    obj.set_solve_method('iterative')
        obj.add_young_modulus('6.5E+10')
        obj.add_poisson_ratio('0.3')
        obj.add_variable('eps', 1.0E-6)
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == 0')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == 2')
        obj.add_concentrated_load(DIR_3, '-1.0e+4*cos(atan2(y,z))', 'abs(y**2 + z**2 - 0.5**2) <= eps')
        obj.add_concentrated_load(DIR_2, '-1.0e+4*sin(atan2(y,z))', 'abs(y**2 + z**2 - 0.5**2) <= eps')
        obj.add_concentrated_load(DIR_3, '1.0e+4*cos(atan2(y,z))', 'abs(y**2 + z**2 - 0.25**2) <= eps')
        obj.add_concentrated_load(DIR_2, '1.0e+4*sin(atan2(y,z))', 'abs(y**2 + z**2 - 0.25**2) <= eps')
        if obj.calc():
            # obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def tank3(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/tank3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('6.5E+10')
        obj.add_poisson_ratio('0.3')
        obj.add_variable('eps', 1.0E-6)
        obj.add_variable('min', 0.0015)
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'y == -0.598 and abs(x**2 + z**2 - 1.6635**2) <= eps')
        obj.add_boundary_condition(DIR_1, '0', 'x == 0')
        obj.add_boundary_condition(DIR_3, '0', 'z == 0')
        obj.add_surface_load(DIR_1, '1.0e+4*cos(atan2(z,x))',
                             '(y <= 0 and y >= -0.2690) and (abs(x**2 + z**2 - (1.03 - min)**2) <= eps)')
        obj.add_surface_load(DIR_3, '1.0e+4*sin(atan2(z,x))',
                             '(y <= 0 and y >= -0.2690) and (abs(x**2 + z**2 - (1.037 - min)**2) <= eps)')
        obj.add_surface_load(DIR_1, '1.0e+4*cos(atan2(z,x))*sin(atan2((x**2 + z**2)**0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x**2 + z**2 + (y + 0.2690)**2 - (1.037-min)**2) <= eps)')
        obj.add_surface_load(DIR_2, '1.0e+4*cos(atan2((x**2+z**2)**0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x**2 + z**2 + (y + 0.2690)**2 - (1.037-min)**2) <= eps)')
        obj.add_surface_load(DIR_3, '1.0e+4*sin(atan2(z,x))*sin(atan2((x**2+z**2)**0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x**2 + z**2 + (y + 0.2690)**2 - (1.037-min)**2) <= eps)')
        obj.add_surface_load(DIR_1, '1.0e+4*cos(atan2(z,x))*sin(atan2((x**2+z**2)**0.5,y))',
                             '(y > 0) and (abs(x**2 + y**2 + z**2 - (1.037-min)**2) <= eps)')
        obj.add_surface_load(DIR_2, '1.0e+4*cos(atan2((x**2+z**2)**0.5,y))',
                             '(y > 0) and (abs(x**2 + y**2 + z**2 - (1.037-min)**2) <= eps)')
        obj.add_surface_load(DIR_3, '1.0e+4*sin(atan2(z,x))*sin(atan2((x**2+z**2)**0.5,y))',
                             '(y > 0) and (abs(x**2 + y**2 + z**2 - (1.037-min)**2) <= eps)')
        obj.add_surface_load(DIR_1, '-5.0e+3*cos(atan2(z,x))', '(y <= 0 and y>=-0.2690) and (abs(x**2 + z**2 - \
        (1.037)**2) <= eps)')
        obj.add_surface_load(DIR_3, '-5.0e+3*sin(atan2(z,x))', '(y <= 0 and y>=-0.2690) and (abs(x**2 + z**2 - \
        (1.037)**2) <= eps)')
        obj.add_surface_load(DIR_1, '-5.0e+3*cos(atan2(z,x))*sin(atan2((x**2+z**2)**0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x**2 + z**2 + (y + 0.2690)**2 - (1.037)**2) <= eps)')
        obj.add_surface_load(DIR_2, '-5.0e+3*cos(atan2((x**2+z**2)**0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x**2 + z**2 + (y + 0.2690)**2 - (1.037)**2) <= eps)')
        obj.add_surface_load(DIR_3, '-5.0e+3*sin(atan2(z,x))*sin(atan2((x**2+z**2)**0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x**2 + z**2 + (y + 0.2690)**2 - (1.037)**2) <= eps)')
        obj.add_surface_load(DIR_2, '-5.e+3', '(y == -1.7235) and (x**2+z**2 - 0.34205**2 <= eps)')
        obj.add_surface_load(DIR_2, '-5.e+3', '(y == -1.944) and (x**2 + z**2 - 0.657857**2 <= eps and x**2 + '
                                              'z**2 - 0.562143**2 >= eps)')
        obj.add_surface_load(DIR_1, '5.0e+3*cos(atan2(z,x))',
                             'abs(y+0.6431) <= eps and abs(x**2 + z**2 - 1.6389**2) <= eps')
        obj.add_surface_load(DIR_3, '5.0e+3*sin(atan2(z,x))',
                             'abs(y+0.6431) <= eps and abs(x**2 + z**2 - 1.6389**2) <= eps')
        obj.add_surface_load(DIR_1, '5.0e+3*x*(1.0644108554**2)/(((x*(1.0644108554**2))**2+(y+1.1013629509)**2 + '
                             '(z*(1.0644108554**2))**2)**0.5)',
                             '(y>-0.6431 and y <-0.0234) and abs(y-((x**2+z**2)**0.5)*(-1.0644108554) - '
                             '1.1013629509) <= eps')
        obj.add_surface_load(DIR_2, '5.0e+3*(y+1.1013629509)/(((x*(1.0644108554**2))**2+(y+1.1013629509)**2+\
        (z*(1.0644108554**2))**2)**0.5)', '(y>-0.6431 and y <-0.0234) and abs(y-((x**2+z**2)**0.5)*\
        (-1.0644108554)-1.1013629509) <= eps')
        obj.add_surface_load(DIR_3, '5.0e+3*z*(1.0644108554**2)/(((x*(1.0644108554**2))**2+(y+1.1013629509)**2+\
        (z*(1.0644108554**2))**2)**0.5)', '(y>-0.6431 and y <-0.0234) and abs(y-((x**2+z**2)**0.5)*\
        (-1.0644108554)-1.1013629509) <= eps')
        obj.add_surface_load(DIR_1, '-5.0e+3*x*(1.0018498686**2)/(((x*(1.0018498686**2))**2+(z*(1.0018498686**2))**2+\
        (y-1.3808172524)**2)**0.5)', '(y>-1.944 and y <-1.7235) and abs(y - ((x**2+z**2)**0.5)*\
        (-1.0018498686)+1.3808172524) <= eps')
        obj.add_surface_load(DIR_2, '5.0e+3*(y-1.3808172524)/(((x*(1.0018498686**2))**2+(z*(1.0018498686**2))**2+\
        (y-1.3808172524)**2)**0.5)', '(y>-1.944 and y <-1.7235) and abs(y - ((x**2+z**2)**0.5)*(-1.0018498686)+\
        1.3808172524) <= eps')
        obj.add_surface_load(DIR_3, '-5.0e+3*z*(1.0018498686**2)/(((x*(1.0018498686**2))**2+(z*(1.0018498686**2))**2+\
        (y-1.3808172524)**2)**0.5)', '(y>-1.944 and y <-1.7235) and abs(y - ((x**2+z**2)**0.5)*(-1.0018498686)+\
        1.3808172524) <= eps')
        obj.add_surface_load(DIR_1, '5.0e+3*x*(1.3260378897**2)/(((3*x*(1.3260378897**2))**2+(y-2.8163434974)**2+\
        (3*z*(1.3260378897**2))**2)**0.5)', '(y>-1.944 and y < -0.6431) and abs(y-((x**2+z**2)**0.5)*(1.3260378897)+\
        2.8163434974) <= eps')
        obj.add_surface_load(DIR_2, '5.0e+3*(y-2.8163434974)/(((3*x*(1.3260378897**2))**2+(y-2.8163434974)**2+\
        (3*z*(1.3260378897**2))**2)**0.5)', '(y>-1.944 and y < -0.6431) and abs(y-((x**2+z**2)**0.5)*(1.3260378897)+\
        2.8163434974) <= eps')
        obj.add_surface_load(DIR_3, '5.0e+3*z*(1.3260378897**2)/(((3*x*(1.3260378897**2))**2+(y-2.8163434974)**2+\
        (3*z*(1.3260378897**2))**2)**0.5)', '(y>-1.944 and y < -0.6431) and abs(y-((x**2+z**2)**0.5)*(1.3260378897)+\
        2.8163434974) <= eps')
        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def head3d(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/head3d.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
    #    obj.set_solve_method('iterative')
        obj.add_young_modulus('1000')
        obj.add_poisson_ratio('0.3')
        obj.add_variable('eps', 1.0E-6)
        obj.add_boundary_condition(DIR_1 | DIR_3, '0', 'y == 0')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'y == 991.3')
        obj.add_surface_load(DIR_1, '-1*cos(atan2(z,x))', 'abs(x**2 + z**2 - 210**2) <=0.001')
        obj.add_surface_load(DIR_3, '-1*sin(atan2(z,x))', 'abs(x**2 + z**2 - 210**2) <= 0.001')
        if obj.calc():
            obj.save_result('head3d')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def console_dynamic(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/console.trpa'):
        obj.set_problem_type('dynamic')
        obj.set_solve_method('direct')
        #    obj.set_solve_method('iterative')
        obj.add_density('1.0E+3')
        obj.add_damping('3.383')
        obj.set_time(0, 1.0, 0.25)
        obj.add_young_modulus('6.5E+10')
        obj.add_poisson_ratio('0.3')
        obj.add_boundary_condition(DIR_1 | DIR_2, '0', 'x == 0')
        obj.add_concentrated_load(DIR_1, '-1.0E+5*cos(t)', 'x == 10')
        obj.add_initial_condition(INIT_U, '0')
        obj.add_initial_condition(INIT_V, '0')
        obj.add_initial_condition(INIT_U_T, '0')
        obj.add_initial_condition(INIT_V_T, '0')
        obj.add_initial_condition(INIT_U_T_T, '0')
        obj.add_initial_condition(INIT_V_T_T, '0')
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
        obj.add_young_modulus('2E+6')
        obj.add_poisson_ratio('0.3')
#        obj.set_names(['W', 'Tx', 'Ty', 'Exx', 'Eyy', 'Exy', 'Sxx', 'Syy', 'Sxy'])
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == -0.1 or x == 0.1 or y == -0.1 or y == 0.1')
        obj.add_surface_load(DIR_1, '-2000')
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
        obj.add_young_modulus('2E+6')
        obj.add_poisson_ratio('0.3')
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == -0.1 or x == 0.1 or y == -0.1 or y == 0.1')
#        obj.add_boundary_condition(DIR_1, '0', 'x = -0.1 or x = 0.1 or y = -0.1 or y = 0.1')
#        obj.add_concentrated_load(DIR_1, '-1.0E+5', 'x = 0 and y = 0')
        obj.add_surface_load(DIR_1, '-2000')
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
        obj.add_young_modulus('2E+6')
        obj.add_poisson_ratio('0.3')
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_2, '0', 'y = 0')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z = 0')
        obj.add_concentrated_load(DIR_2, '-1.0E+2', 'z = 5')
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
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x = -0.1 or x = 0.1 or y = -0.1 or y = 0.1')
        obj.add_surface_load(DIR_1, '-2000')
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
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x = -0.1 or x = 0.1 or y = -0.1 or y = 0.1')
#        obj.add_concentrated_load(DIR_1, '-1.0E+5', 'x = 0 and y = 0')
#        obj.add_surface_load(DIR_1, '-2000')
        obj.add_concentrated_load(DIR_1, '-2000', 'z = 0')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell3_test(res_name):
    obj = TObject()
    # if obj.set_mesh('mesh/shell-tube.trpa'):
    if obj.set_mesh('mesh/shell-tube-3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.0369')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z == 0 or z == 4.014')
        obj.add_surface_load(DIR_1, '0.05*cos(atan2(y,x))', '(abs(x**2 + y**2 - 1.99**2) <= 1.0E-3)')
        obj.add_surface_load(DIR_2, '0.05*sin(atan2(y,x))', '(abs(x**2 + y**2 - 1.99**2) <= 1.0E-3)')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell6_test(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/shell-tube6.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.0369')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z == 0 or z == 4.014')
        obj.add_surface_load(DIR_1, '0.05*cos(atan2(y,x))', '(abs(x**2 + y**2 - 1.99**2) <= 1.0E-3)')
        obj.add_surface_load(DIR_2, '0.05*sin(atan2(y,x))', '(abs(x**2 + y**2 - 1.99**2) <= 1.0E-3)')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def plate3_test(res_name):
    obj = TObject()
    # if obj.set_mesh('mesh/plate3_1_0.trpa'):
    if obj.set_mesh('mesh/plate3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_output(20, 10)
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        # obj.add_alpha('1.25E-5')
        # obj.add_temperature('100')
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == -0.5 or x == 0.5 or y == -0.5 or y == 0.5')
        obj.add_surface_load(DIR_1, '0.05')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def plate6_test(res_name):
    obj = TObject()
    # if obj.set_mesh('mesh/plate3_1_0.trpa'):
    if obj.set_mesh('mesh/plate6.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == -0.5 or x == 0.5 or y == -0.5 or y == 0.5')
        obj.add_surface_load(DIR_1, '0.05')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell_plate3_test(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/shell-plate3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == -0.5 or x == 0.5 or y == -0.5 or y == 0.5')
        obj.add_surface_load(DIR_3, '0.05')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell_plate6_test(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/shell-plate6.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.01')
        # obj.add_alpha('1.25E-5')
        # obj.add_temperature('100')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == -0.5 or x == 0.5 or y == -0.5 or y == 0.5')
        obj.add_surface_load(DIR_3, '0.05')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell_plate4_test(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/shell-plate4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.01')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == -0.5 or x == 0.5 or y == -0.5 or y == 0.5')
        obj.add_surface_load(DIR_3, '0.05')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def plate4_test(res_name):
    obj = TObject()
    # if obj.set_mesh('mesh/plate4-1.0.trpa'):
    if obj.set_mesh('mesh/plate4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == -0.5 or x == 0.5 or y == -0.5 or y == 0.5')
        obj.add_surface_load(DIR_1, '0.05')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def shell4_test(res_name):
    obj = TObject()
    # if obj.set_mesh('mesh/shell4-2.0.trpa'):
    if obj.set_mesh('mesh/shell4_1_0.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.0369')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z == 0 or z == 4.014')
        obj.add_surface_load(DIR_1, '0.05*cos(atan2(y,x))', '(abs(x**2 + y**2 - 1.99**2) <= 1.0E-3)')
        obj.add_surface_load(DIR_2, '0.05*sin(atan2(y,x))', '(abs(x**2 + y**2 - 1.99**2) <= 1.0E-3)')
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
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'z == 0 or z == 4.014')
        obj.add_surface_load(DIR_1, '0.05*cos(atan2(y,x))', '(abs(x**2 + y**2 - 1.99**2) <= 1.0E-3)')
        obj.add_surface_load(DIR_2, '0.05*sin(atan2(y,x))', '(abs(x**2 + y**2 - 1.99**2) <= 1.0E-3)')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def tank3s(res_name):
    obj = TObject()
    # if obj.set_mesh('mesh/tank6s.trpa'):
    if obj.set_mesh('mesh/tank3s.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
#        obj.set_solve_method('iterative')
        obj.add_variable('eps', 0.01)
        obj.add_variable('L', 0.269)
        obj.add_variable('P', 10000)
        obj.add_variable('p', 5000)
        obj.add_variable('R', 1.037)

        obj.add_young_modulus('6.5E+10')
        obj.add_poisson_ratio('0.3')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')

        # obj.add_thickness('0.025', 'abs(y + 1.724) <= eps and (x**2+z**2 - 0.342**2 <= eps)')
        obj.add_thickness('0.0015')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'y == -0.643 and abs(x**2 + z**2 -1.641**2) <= eps')
        # obj.add_boundary_condition(DIR_1, '0', 'abs(x) <= eps')
        # obj.add_boundary_condition(DIR_3, '0', 'abs(z) <= eps')

        obj.add_surface_load(DIR_1, 'P*cos(atan2(z,x))', '(y <= 0 and y>=-L) and (abs(x**2 + z**2 - R**2) <= eps)')
        obj.add_surface_load(DIR_3, 'P*sin(atan2(z,x))', '(y <= 0 and y>=-L) and (abs(x**2 + z**2 - R**2) <= eps)')

        obj.add_surface_load(DIR_1, 'P*cos(atan2(z,x))*sin(atan2((x**2 + z**2)**0.5,(y + L)))',
                             '(y < -L) and (abs(x**2 + z**2 + (y + L)**2 - R**2) <= eps)')
        obj.add_surface_load(DIR_2, 'P*cos(atan2((x**2 + z**2)**0.5,(y + L)))',
                             '(y < -L) and (abs(x**2 + z**2 + (y + L)**2 - R**2) <= eps)')
        obj.add_surface_load(DIR_3, 'P*sin(atan2(z,x))*sin(atan2((x**2 + z**2)**0.5,(y + L)))',
                             '(y < -L) and (abs(x**2 + z**2 + (y + L)**2 - R**2) <= eps)')

        obj.add_surface_load(DIR_1, 'P*cos(atan2(z,x))*sin(atan2((x**2 + z**2)**0.5,y))',
                             '(y > 0) and (abs(x**2 + y**2 + z**2 - R**2) <= eps)')
        obj.add_surface_load(DIR_2, 'P*cos(atan2((x**2 + z**2)**0.5,y))',
                             '(y > 0) and (abs(x**2 + y**2 + z**2 - R**2) <= eps)')
        obj.add_surface_load(DIR_3, 'P*sin(atan2(z,x))*sin(atan2((x**2 + z**2)**0.5,y))',
                             '(y > 0) and (abs(x**2 + y**2 + z**2 - R**2) <= eps)')

        obj.add_surface_load(DIR_1, '-p*cos(atan2(z,x))', '(y <= 0 and y>=-L) and (abs(x**2 + z**2 - R**2) <= eps)')
        obj.add_surface_load(DIR_3, '-p*sin(atan2(z,x))', '(y <= 0 and y>=-L) and (abs(x**2 + z**2 - R**2) <= eps)')

        obj.add_surface_load(DIR_1, '-p*cos(atan2(z,x))*sin(atan2((x**2 + z**2)**0.5,(y + L)))',
                             '(y < -L) and (abs(x**2 + z**2 + (y + L)**2 - R**2) <= eps)')
        obj.add_surface_load(DIR_2, '-p*cos(atan2((x**2 + z**2)**0.5,(y + L)))',
                             '(y < -L) and (abs(x**2 + z**2 + (y + L)**2 - R**2) <= eps)')
        obj.add_surface_load(DIR_3, '-p*sin(atan2(z,x))*sin(atan2((x**2 + z**2)**0.5,(y + L)))',
                             '(y < -L) and (abs(x**2 + z**2 + (y + L)**2 - R**2) <= eps)')

        obj.add_surface_load(DIR_2, '-p', '(y == -1.724) and (x**2+z**2 - 0.342**2 <= eps)')
        obj.add_surface_load(DIR_2, '-p', '(y == -1.944) and (x**2+z**2 - 0.660**2 <= eps)')

        obj.add_surface_load(DIR_1, 'p*cos(atan2(z,x))', 'abs(y + 0.641) <= eps and abs(x**2 + z**2 - 1.643**2) <= eps')
        obj.add_surface_load(DIR_3, 'p*sin(atan2(z,x))', 'abs(y + 0.641) <= eps and abs(x**2 + z**2 - 1.643**2) <= eps')

        obj.add_surface_load(DIR_1, 'p*x*(1.0644108554**2)/(((x*(1.0644108554**2))**2+(y + 1.1013629509)**2 + '
                             '(z*(1.0644108554**2))**2)**0.5)', '(y > -0.641 and y <-0.0234) and '
                             'abs(y-((x**2+z**2)**0.5)*(-1.0644108554)-1.1013629509)<=eps')
        obj.add_surface_load(DIR_2, 'p*(y+1.1013629509)/(((x*(1.0644108554**2))**2 + (y + 1.1013629509)**2 + '
                             '(z*(1.0644108554**2))**2)**0.5)', '(y > -0.6431 and y < -0.0234) and '
                             'abs(y-((x**2+z**2)**0.5)*(-1.0644108554)-1.1013629509)<=eps')
        obj.add_surface_load(DIR_3, 'p*z*(1.0644108554**2)/(((x*(1.0644108554**2))**2 + (y + 1.1013629509)**2 + '
                             '(z*(1.0644108554**2))**2)**0.5)', '(y>-0.6431 and y <-0.0234) and '
                             'abs(y-((x**2+z**2)**0.5)*(-1.0644108554)-1.1013629509)<=eps')

        obj.add_surface_load(DIR_1, '-p*x*(1.0018498686**2)/(((x*(1.0018498686**2))**2 + (z*(1.0018498686**2))**2 + '
                             '(y-1.3808172524)**2)**0.5)', '(y>-1.944 and y <-1.7235) and '
                             'abs(y - ((x**2 + z**2)**0.5)*(-1.0018498686)+1.3808172524)<=eps')
        obj.add_surface_load(DIR_2, 'p*(y-1.3808172524)/(((x*(1.0018498686**2))**2 + (z*(1.0018498686**2))**2 + '
                             '(y - 1.3808172524)**2)**0.5)', '(y>-1.944 and y <-1.7235) and '
                             'abs(y - ((x**2+z**2)**0.5)*(-1.0018498686)+1.3808172524)<=eps')
        obj.add_surface_load(DIR_3, '-p*z*(1.0018498686**2)/(((x*(1.0018498686**2))**2 + (z*(1.0018498686**2))**2 + '
                             '(y - 1.3808172524)**2)**0.5)', '(y>-1.944 and y <-1.7235) and '
                             'abs(y - ((x**2+z**2)**0.5)*(-1.0018498686)+1.3808172524)<=eps')

        obj.add_surface_load(DIR_1, 'p*x*(1.3260378897**2)/(((3*x*(1.3260378897**2))**2 + (y - 2.8163434974)**2 + '
                             '(3*z*(1.3260378897**2))**2)**0.5)', '(y>-1.944 and y < -0.6431) and '
                             'abs(y-((x**2+z**2)**0.5)*(1.3260378897)+2.8163434974)<=eps')
        obj.add_surface_load(DIR_2, 'p*(y-2.8163434974)/(((3*x*(1.3260378897**2))**2 + (y - 2.8163434974)**2 + '
                             '(3*z*(1.3260378897**2))**2)**0.5)', '(y>-1.944 and y < -0.6431) and '
                             'abs(y-((x**2+z**2)**0.5)*(1.3260378897)+2.8163434974)<=eps')
        obj.add_surface_load(DIR_3, 'p*z*(1.3260378897**2)/(((3*x*(1.3260378897**2))**2+(y - 2.8163434974)**2 + '
                             '(3*z*(1.3260378897**2))**2)**0.5)', '(y>-1.944 and y < -0.6431) and '
                             'abs(y - ((x**2 + z**2)**0.5)*(1.3260378897) + 2.8163434974)<=eps')

        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def plate3d(res_name):
    obj = TObject()
    # if obj.set_mesh('mesh/plate3d.trpa'):
    if obj.set_mesh('mesh/plate3d-1.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200000000')
        obj.add_poisson_ratio('0.27')
        obj.add_temperature('100')
        obj.add_alpha('1.25E-5')

        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'x == -0.5 or x == 0.5 or y == -0.5 or y == 0.5')
        # obj.add_surface_load(DIR_3, '-50000', 'z == 0')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def quad4(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/quad-4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_boundary_condition(DIR_1 | DIR_2, '0', 'y == -0.5')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        # obj.add_volume_load(DIR_2, '-0.05')
        obj.add_surface_load(DIR_2, '-0.05', 'y == 0.5')
        # obj.add_concentrated_load(DIR_2, '-0.05', 'y = 0.5 and (x = -0.5 or x = 0.5)')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def quad3(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/quad-3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_boundary_condition(DIR_1 | DIR_2, '0', 'y == -0.5')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        obj.add_volume_load(DIR_2, '-0.05')
        # obj.add_surface_load(DIR_2, '-0.05', 'y == 0.5')
        # obj.add_concentrated_load(DIR_2, '-0.05', 'y = 0.5 and (x = -0.5 or x = 0.5)')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def quad6(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/quad-6.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        # obj.add_temperature('100')
        # obj.add_alpha('1.25E-5')
        obj.add_boundary_condition(DIR_1 | DIR_2, '0', 'y == -0.5')
        obj.add_volume_load(DIR_2, '-0.05')
        # obj.add_surface_load(DIR_2, '-0.05', 'y == 0.5')
        # obj.add_concentrated_load(DIR_2, '-0.05', 'y = 0.5 and (x = -0.5 or x = 0.5)')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def beam2d3(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/beam2d-3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_boundary_condition(DIR_2, '0', '(x == -5 and y == -0.25) or (x == 5 and y == -0.25)')
        # obj.add_boundary_condition(DIR_1 | DIR_2, '0', '(x == -5) or (x == 5)')
        obj.add_thickness('0.01')
        obj.add_surface_load(DIR_2, '-1', 'y == 0.25')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def beam3d4(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/beam3d-10.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.add_young_modulus('203200')
        obj.add_poisson_ratio('0.27')
        obj.add_thickness('0.01')
        obj.add_boundary_condition(DIR_2, '0', '(x == -5 and y == -0.25) or (x == 5 and y == -0.25)')
        # obj.add_boundary_condition(DIR_1 | DIR_2, '0', '(x == -5) or (x == 5)')
        obj.add_surface_load(DIR_2, '-1', 'y == 0.25')
        if obj.calc():
            obj.print_result()
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


def tank3ds(res_name):
    obj = TObject()
    if obj.set_mesh('mesh/tank3ds.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
#        obj.set_solve_method('iterative')

        obj.add_young_modulus('6.5E+10', 'z < -7')
        obj.add_poisson_ratio('0.3', 'z < -7')
        obj.add_young_modulus('9.5E+10', 'z >= -7')
        obj.add_poisson_ratio('0.3', 'z >= -7')

        # obj.add_young_modulus('6.5E+10')
        # obj.add_poisson_ratio('0.3')

        obj.add_variable('p', 1000)         # Давление
        obj.add_variable('eps', 1.0E-2)
        obj.add_variable('l', 16.691)       # Высота обечайки
        obj.add_variable('h', 17.626)       # Высота бака
        obj.add_variable('r', 2.5)          # Радиус днищ
        obj.add_variable('d', 3.9)          # Диаметр бака
        obj.add_variable('c0', -1.565)      # Координата z центра верхнего днища
        obj.add_variable('c1', -15.126)     # ... нижнего днища

        obj.add_thickness('0.0028')
        obj.add_boundary_condition(DIR_1 | DIR_2 | DIR_3, '0', 'abs(z + h) <= eps')
        obj.add_boundary_condition(DIR_1, '0', 'abs(x) <= eps')
        obj.add_boundary_condition(DIR_2, '0', 'abs(y) <= eps')

        obj.add_surface_load(DIR_1, 'p*cos(atan2(y,x))', '(z <= 0 and z >= -l)', )
        obj.add_surface_load(DIR_2, 'p*sin(atan2(y,x))', '(z <= 0 and z >= -l)')

        obj.add_surface_load(DIR_1, 'p*cos(atan2(y,x))*sin(atan2((x**2+y**2)**0.5,(z -c0)))',
                             'abs(x**2 + y**2 + (z - c0)**2 - r**2) <= eps')
        obj.add_surface_load(DIR_2, 'p*sin(atan2(y,x))*sin(atan2((x**2+y**2)**0.5,(z - c0)))',
                             'abs(x**2 + y**2 + (z - c0)**2 - r**2) <= eps')
        obj.add_surface_load(DIR_3, 'p*cos(atan2((x**2+y**2)**0.5,(z - c0)))',
                             'abs(x**2 + y**2 + (z - c0)**2 - r**2) <= eps')

        obj.add_surface_load(DIR_1, 'p*cos(atan2(y,x))*sin(atan2((x**2+y**2)**0.5,(z -c1)))',
                             'abs(x**2 + y**2 + (z - c1)**2 - r**2) <= eps')
        obj.add_surface_load(DIR_2, 'p*sin(atan2(y,x))*sin(atan2((x**2+y**2)**0.5,(z - c1)))',
                             'abs(x**2 + y**2 + (z - c1)**2 - r**2) <= eps')
        obj.add_surface_load(DIR_3, 'p*cos(atan2((x**2+y**2)**0.5,(z - c1)))',
                             'abs(x**2 + y**2 + (z - c1)**2 - r**2) <= eps')

        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.save_result(res_name)
            TPlot(res_name)
            return True
        return False


if __name__ == '__main__':
    # -------------- 1d -----------------
    # body1d('body1d')

    # -------------- 2d -----------------
    # quad4('quad-4')
    # quad3('quad-3')
    # quad6('quad-6')
    # console('console')
    # console4('console4')
    # quad('quad')

    # -------------- 3d -----------------
    # cube('cube')
    # cube10('cube-10')
    # cube4('cube-4')

    # beam('beam')
    # head3d('head3d')
    # tank3('tank3')
    # cylinder('cylinder')
    # rod4('rod4')
    # rod10('rod10')
    # beam3d4('beam3d-4')
    # beam2d3('beam2d-3')

    # tube_test('tube_test')
    # plate3d('plate3d')

    # ------------- Plate ----------------
    # plate3('plate3')
    # plate3_test('plate3_test')
    # plate6_test('plate6_test')
    plate4_test('plate4_test')
    # plate4('plate4')

    # -------------- Shell ---------------
    # shell3('shell3')
    # shell4('shell4')
    # shell_plate3('shell_plate3')

    # shell_plate3_test('shell-plate3')
    # shell_plate6_test('shell-plate3')
    # shell_plate4_test('shell-plate4')

    # shell3_test('shell3_test')
    # shell4_test('shell4_test')
    # shell6_test('shell6_test')

    # tank3ds('tank3ds')
    # tank3s('tank3s')

    # -------------- Dynamic -----------------
    # console_dynamic('console_dynamic')
    # beam_dynamic('beam_dynamic')

    # TPlot('tank3ds')
    # TPlot('console_dynamic')
'''
2. Правильно отображать динамическую задачу в plot3d
3. Рисовать оси координат
4. Правильно рисовать маленькие величины (стержень) 
'''
