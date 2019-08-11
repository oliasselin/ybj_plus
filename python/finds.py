#!/usr/bin/env python
import os
import re
import numpy as np

def find_resolution(run):
    #Read parameters.f90 from the original filee:                                                                                                                              
    with open (run+'/source/parameters.f90', 'rt') as in_file:  # Open file lorem.txt for reading of text data.                                                      
        for line in in_file: # Store each line in a string variable "line"                                                                                                           
            if line.find('integer, parameter :: n1') != -1:
                n1=int(re.findall(r'\d+',line)[1])
                n2=int(re.findall(r'\d+',line)[3])
                n3=int(re.findall(r'\d+',line)[5])
                break

    return n1,n2,n3


def find_scales(run):
    #Read parameters.f90 from the original filee:                                                                                                                              
    with open (run+'/source/parameters.f90', 'rt') as in_file:  # Open file lorem.txt for reading of text data.                                                      
        for line in in_file: # Store each line in a string variable "line"                                                                                                           
            if line.find('    double precision, parameter :: dom_x') != -1:
                Dx=float(re.findall(r'\d+',line)[0])
                L_scale=Dx/(2.*np.pi)
            if line.find('    double precision, parameter :: dom_z') != -1:
                Dz=float(re.findall(r'\d+',line)[0])
                H_scale=Dz/(2.*np.pi)
            if line.find('    double precision, parameter :: U_scale') != -1:
                U_scale=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
            if line.find('    double precision, parameter :: N2_scale') != -1:
                h_thermo=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])  #Not robust...
                
                
    return Dx,Dz,L_scale,H_scale,U_scale,h_thermo

def find_timestep(run):
    #Read runs_spec.dat from the original file:                                                                                                                            
    with open (run+'/output/run_specs.dat', 'rt') as in_file:  # Open file lorem.txt for reading of text data.                                                                          
        for line in in_file: # Store each line in a string variable "line"                                                                                                              
            if line.find(' dt  =') != -1:
                delt=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
            if line.find(' Period of total energy        output:') != -1:
                freq_etot=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])
            if line.find(' Period of flow z-profile      output:') != -1:
                freq_ez=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])
            if line.find(' Period of total wave energy   output:') != -1:
                freq_we=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])
            if line.find(' Period of wave z-profile      output:') != -1:
                freq_wz=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])



    return delt,freq_etot,freq_ez,freq_we,freq_wz


def find_nondim(run):
    #Read runs_spec.dat from the original file:                                                                                                                            
    with open (run+'/output/run_specs.dat', 'rt') as in_file:  # Open file lorem.txt for reading of text data.                                                                          
        for line in in_file: # Store each line in a string variable "line"                                                                                                              
            if line.find(' Ro =   ') != -1:
                Ro=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
                Fr=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])
                Bu = Fr*Fr/(Ro*Ro)


    with open (run+'/source/parameters.f90', 'rt') as in_file:  # Open file lorem.txt for reading of text data.                                                                                                           
        for line in in_file: # Store each line in a string variable "line"  
            if line.find('    double precision, parameter :: YBJ_criterion') != -1:
                YBJ_criterion=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])


    return Ro,Fr,Bu,YBJ_criterion

