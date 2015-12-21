#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fem_defs import DIR_X, DIR_Y, DIR_Z
from fem_error import TFEMException
from fem_object import TObject


from scipy.sparse import csr_matrix, lil_matrix

def body1d():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/body1d.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
    obj.set_elasticity(e, m)
    obj.add_boundary_condition('0', 'x=0', DIR_X)
    obj.add_volume_load('-1.0E+5', '', DIR_X)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result()


def cube():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/cube.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
    obj.set_elasticity(e, m)
    obj.add_boundary_condition('0', 'z=0', DIR_X | DIR_Y | DIR_Z)
#    obj.add_volume_load('-1.0E+5', '', DIR_Z)
    obj.add_surface_load('-1.0E+5', '', DIR_Z)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result()


def beam():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/beam.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('iterative')
#    obj.set_solve_method('direct')
    obj.set_elasticity(e, m)
    obj.add_boundary_condition('0', 'y=0', DIR_X | DIR_Y | DIR_Z)
    obj.add_volume_load('-1.0E+5', '', DIR_Y)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result()


def console():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/console.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
#    obj.set_solve_method('iterative')
    obj.set_elasticity(e, m)
    obj.add_boundary_condition('0', 'x=0', DIR_X | DIR_Y)
    obj.add_volume_load('-1.0E+5', '', DIR_X)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result()


def cylinder():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/cyl.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
#    obj.set_solve_method('iterative')
    obj.set_elasticity(e, m)
    obj.add_variable('eps', 1.0E-6)
    obj.add_boundary_condition('0', 'x=0', DIR_X | DIR_Y | DIR_Z)
    obj.add_boundary_condition('0', 'x=2', DIR_X | DIR_Y | DIR_Z)
    obj.add_surface_load('-1.0e+4*cos(atan2(y,z))', 'abs(y^2 + z^2 - 0.5^2) <= eps', DIR_Z)
    obj.add_surface_load('-1.0e+4*sin(atan2(y,z))', 'abs(y^2 + z^2 - 0.5^2) <= eps', DIR_Y)
    obj.add_surface_load('1.0e+4*cos(atan2(y,z))', 'abs(y^2 + z^2 - 0.25^2) <= eps', DIR_Z)
    obj.add_surface_load('1.0e+4*sin(atan2(y,z))', 'abs(y^2 + z^2 - 0.25^2) <= eps',DIR_Y)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result()


def tank3():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/tank3.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
    obj.set_elasticity(e, m)
    obj.add_variable('eps', 1.0E-6)
    obj.add_variable('min', 0.0015)
    obj.add_boundary_condition('0', 'y=-0.598 and abs(x^2+z^2-1.6635^2)<=eps', DIR_X | DIR_Y | DIR_Z)
    obj.add_boundary_condition('0', 'x=0', DIR_X)
    obj.add_boundary_condition('0', 'z=0', DIR_Z)
    obj.add_surface_load('1.0e+4*cos(atan2(z,x))','(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037-min)^2) <= eps)',DIR_X)
    obj.add_surface_load('1.0e+4*sin(atan2(z,x))','(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037-min)^2) <= eps)',DIR_Z)
    obj.add_surface_load('1.0e+4*cos(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))','(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037-min)^2) <= eps)',DIR_X)
    obj.add_surface_load('1.0e+4*cos(atan2((x^2+z^2)^0.5,(y + 0.2690)))','(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037-min)^2) <= eps)',DIR_Y)
    obj.add_surface_load('1.0e+4*sin(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))','(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037-min)^2) <= eps)',DIR_Z)
    obj.add_surface_load('1.0e+4*cos(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,y))','(y > 0) and (abs(x^2 + y^2 + z^2 - (1.037-min)^2) <= eps)',DIR_X)
    obj.add_surface_load('1.0e+4*cos(atan2((x^2+z^2)^0.5,y))','(y > 0) and (abs(x^2 + y^2 + z^2 - (1.037-min)^2) <= eps)',DIR_Y)
    obj.add_surface_load('1.0e+4*sin(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,y))','(y > 0) and (abs(x^2 + y^2 + z^2 - (1.037-min)^2) <= eps)',DIR_Z)
    obj.add_surface_load('-5.0e+3*cos(atan2(z,x))','(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037)^2) <= eps)',DIR_X)
    obj.add_surface_load('-5.0e+3*sin(atan2(z,x))','(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037)^2) <= eps)',DIR_Z)
    obj.add_surface_load('-5.0e+3*cos(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))','(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037)^2) <= eps)',DIR_X)
    obj.add_surface_load('-5.0e+3*cos(atan2((x^2+z^2)^0.5,(y + 0.2690)))','(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037)^2) <= eps)',DIR_Y)
    obj.add_surface_load('-5.0e+3*sin(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))','(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037)^2) <= eps)',DIR_Z)
    obj.add_surface_load('-5.e+3','(y=-1.7235) and (x^2+z^2 - 0.34205^2 <= eps)',DIR_Y)
    obj.add_surface_load('-5.e+3','(y=-1.944) and (x^2+z^2 - 0.657857^2 <= eps and x^2+z^2 - 0.562143^2 >= eps)',DIR_Y)
    obj.add_surface_load('5.0e+3*cos(atan2(z,x))','abs(y+0.6431) <= eps and abs(x^2 + z^2 - 1.6389^2) <= eps',DIR_X)
    obj.add_surface_load('5.0e+3*sin(atan2(z,x))','abs(y+0.6431) <= eps and abs(x^2 + z^2 - 1.6389^2) <= eps',DIR_Z)
    obj.add_surface_load('5.0e+3*x*(1.0644108554^2)/(((x*(1.0644108554^2))^2+(y+1.1013629509)^2+(z*(1.0644108554^2))^2)^0.5)','(y>-0.6431 and y <-0.0234) and abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps',DIR_X)
    obj.add_surface_load('5.0e+3*(y+1.1013629509)/(((x*(1.0644108554^2))^2+(y+1.1013629509)^2+(z*(1.0644108554^2))^2)^0.5)','(y>-0.6431 and y <-0.0234) and abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps',DIR_Y)
    obj.add_surface_load('5.0e+3*z*(1.0644108554^2)/(((x*(1.0644108554^2))^2+(y+1.1013629509)^2+(z*(1.0644108554^2))^2)^0.5)','(y>-0.6431 and y <-0.0234) and abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps',DIR_Z)
    obj.add_surface_load('-5.0e+3*x*(1.0018498686^2)/(((x*(1.0018498686^2))^2+(z*(1.0018498686^2))^2+(y-1.3808172524)^2)^0.5)','(y>-1.944 and y <-1.7235) and abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+1.3808172524)<=eps',DIR_X)
    obj.add_surface_load('5.0e+3*(y-1.3808172524)/(((x*(1.0018498686^2))^2+(z*(1.0018498686^2))^2+(y-1.3808172524)^2)^0.5)','(y>-1.944 and y <-1.7235) and abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+1.3808172524)<=eps',DIR_Y)
    obj.add_surface_load('-5.0e+3*z*(1.0018498686^2)/(((x*(1.0018498686^2))^2+(z*(1.0018498686^2))^2+(y-1.3808172524)^2)^0.5)','(y>-1.944 and y <-1.7235) and abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+1.3808172524)<=eps',DIR_Z)
    obj.add_surface_load('5.0e+3*x*(1.3260378897^2)/(((3*x*(1.3260378897^2))^2+(y-2.8163434974)^2+(3*z*(1.3260378897^2))^2)^0.5)','(y>-1.944 and y < -0.6431) and abs(y-((x^2+z^2)^0.5)*(1.3260378897)+2.8163434974)<=eps',DIR_X)
    obj.add_surface_load('5.0e+3*(y-2.8163434974)/(((3*x*(1.3260378897^2))^2+(y-2.8163434974)^2+(3*z*(1.3260378897^2))^2)^0.5)','(y>-1.944 and y < -0.6431) and abs(y-((x^2+z^2)^0.5)*(1.3260378897)+2.8163434974)<=eps',DIR_Y)
    obj.add_surface_load('5.0e+3*z*(1.3260378897^2)/(((3*x*(1.3260378897^2))^2+(y-2.8163434974)^2+(3*z*(1.3260378897^2))^2)^0.5)','(y>-1.944 and y < -0.6431) and abs(y-((x^2+z^2)^0.5)*(1.3260378897)+2.8163434974)<=eps',DIR_Z)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result('mesh/' + obj.object_name() + '.res')

try:
    # body1d()
    # cube()
    # console()
    # beam()
    tank3()
    # cylinder()

    #A = lil_matrix([[1, 2, 0, 0, 0, 5.2],
    #                [0, 0, 3, 3, 3, 1],
    #                [4, 0, 5, 0, 1, 0],
    #                [4, 0, 5, 0, 1, 0],
    #                [4, 0, 5, 0, 1, 0],
    #                [4, 0, 5, 0, 1, 0]
    #                ])

    #print(A.nonzero())
    #for i in range(0, 6):
    #    for k in range(0, len(A[i].nonzero()[1])):
    #        row = i
    #        col = A[i].nonzero()[1][k]
    #        print('(', row, col, ')', A[row, col])

    #A = lil_matrix([[1, 2, 0, 0, 0, 5],
    #                [2, 0, 3, 3, 3, 1],
    #                [0, 3, 5, 0, 1, 0],
    #                [0, 3, 0, 4, 1, 0],
    #                [0, 3, 1, 1, 1, 1],
    #                [5, 1, 0, 0, 1, 9]
    #                ])

    #print(A.nonzero())
    #for i in range(0, 6):
    #    for k in A[i].nonzero()[1]:
    #        row = i
    #        col = k
    #        print('(', row, col, ')', A[row, col])

    #l = 2
    #for k in A[l].nonzero()[1]:
    #    row = l
    #    col = k
    #    if row != col:
    #        A[row, col] = A[col, row] = 0

    #for i in range(0, 6):
    #    print(A[i,0], A[i,1], A[i,2], A[i,3], A[i,4], A[i,5])



except TFEMException as err:
    err.print_error()

