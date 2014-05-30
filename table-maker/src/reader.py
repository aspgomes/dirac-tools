# Copyright (c) 2014, Andre Severo Pereira Gomes
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the {organization} nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import re
import math
import sys
from electronic_state import *

class fscc_results:
    """
    fscc_results fetched the relevant information from a dirac fock-space cc
    output
    """

    def __init__ (self) :
        self.states  = [ ]

    def process_output (self,outfile):
#
# what we want to read
#
        version_str       = "Release DIRAC([0-9]+)" 

        start_str         = "\s+Solving equations for sector ([0-9]+)"
        state_str         = "\s+Irrep\s+(-?[A-Za-z0-9]+)\s+State\s+(\d+)\s+([-.0-9]+)\s+([-.0-9]+)"

        det_s11_r_str     = "^\s+(-?\d+\.\d+)\s+(-?[0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\) ->\s+(-?[0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\)"
        det_s11_c_str     = "^\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?[0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\) ->\s+(-?[0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\)"

        det_s01_s10_r_str = "^\s+(-?\d+\.\d+)\s+\|\s+(-?[0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\)"
        det_s01_s10_c_str = "^\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\|\s+(-?[0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\)"

        det_s02_s20_r_str = "^\s+(-?\d+\.\d+)\s+\|\s+(-?[0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\)\, \s+(-?[-0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\)"
        det_s02_s20_c_str = "^\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\|\s+(-?[-0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\)\, \s+(-?[-0-9a-zA-Z]+) #\s+(\d+) \(\s*(-?\d*\.\d+)\)"

        version_re        = re.compile (r''+ version_str   +'', re.IGNORECASE)
        start_re          = re.compile (r''+ start_str   +'', re.IGNORECASE)
        state_re          = re.compile (r''+ state_str   +'', re.IGNORECASE)
        det_s11_r_re      = re.compile (r''+ det_s11_r_str +'', re.IGNORECASE)
        det_s11_c_re      = re.compile (r''+ det_s11_c_str +'', re.IGNORECASE)
        det_s01_s10_r_re  = re.compile (r''+ det_s01_s10_r_str +'', re.IGNORECASE)
        det_s01_s10_c_re  = re.compile (r''+ det_s01_s10_c_str +'', re.IGNORECASE)
        det_s02_s20_r_re  = re.compile (r''+ det_s02_s20_r_str +'', re.IGNORECASE)
        det_s02_s20_c_re  = re.compile (r''+ det_s02_s20_c_str +'', re.IGNORECASE)
#
# now we get to work
#
        program_version = 12

        print "reading fscc output from file: ",outfile,"\n"
        f = file(outfile,'r')
        lines = f.readlines()
        read_state = False
        for i, l in enumerate(lines) :
            if version_re.match (l) :
                program_version = version_re.match(l).group(1)
                print "Dirac version recognized as ",program_version
            if start_re.match (l) :
                current_sector = start_re.match(l).group(1)
                read_state = True
            if read_state :
                ca   = 0.0
                cb   = 0.0
                e_h  = 0
                e_p  = 0
                s_h  = ""
                s_p  = ""
                if state_re.match (l) :
                    label  = state_re.match(l).group(1)
                    index  = int(state_re.match(l).group(2))

                    if program_version <= 11 :
                       energy = float(state_re.match(l).group(4))
                       rel_en = float(state_re.match(l).group(3))
                    else :
                       energy = float(state_re.match(l).group(3))
                       rel_en = float(state_re.match(l).group(4))

                    state  = electronic_state()
                    state.set_sector(current_sector)
                    state.set_label(label)
                    state.set_index(index)
                    state.set_energy(energy)
                    state.set_relative_energy(rel_en)
                    self.states.append(state)
                elif det_s11_r_re.match (l) :
                    ca  = float(det_s11_r_re.match(l).group(1))
                    cb  = 0.0
                    i_h = float(det_s11_r_re.match(l).group(3))
                    i_p = float(det_s11_r_re.match(l).group(6))
                    e_h = float(det_s11_r_re.match(l).group(4))
                    e_p = float(det_s11_r_re.match(l).group(7))
                    s_h = det_s11_r_re.match(l).group(2)
                    s_p = det_s11_r_re.match(l).group(5)
                    state.add_determinant(ca,cb,i_h,e_h,i_p,e_p,s_h,s_p)
                elif det_s11_c_re.match (l) :
                    ca  = float(det_s11_c_re.match(l).group(1)) 
                    cb  = float(det_s11_c_re.match(l).group(2)) 
                    i_h = float(det_s11_c_re.match(l).group(4))
                    i_p = float(det_s11_c_re.match(l).group(7))
                    e_h = float(det_s11_c_re.match(l).group(5))
                    e_p = float(det_s11_c_re.match(l).group(8))
                    s_h = det_s11_c_re.match(l).group(3)
                    s_p = det_s11_c_re.match(l).group(6)
                    state.add_determinant(ca,cb,i_h,e_h,i_p,e_p,s_h,s_p)
                elif det_s02_s20_r_re.match (l) :
                    ca  = float(det_s02_s20_r_re.match(l).group(1))
                    cb  = 0.0
                    i_h = float(det_s02_s20_r_re.match(l).group(3))
                    i_p = float(det_s02_s20_r_re.match(l).group(6))
                    e_h = float(det_s02_s20_r_re.match(l).group(4))
                    e_p = float(det_s02_s20_r_re.match(l).group(7))
                    s_h = det_s02_s20_r_re.match(l).group(2)
                    s_p = det_s02_s20_r_re.match(l).group(5)
                    state.add_determinant(ca,cb,i_h,e_h,i_p,e_p,s_h,s_p)
                elif det_s02_s20_c_re.match (l) :
                    ca  = float(det_s02_s20_c_re.match(l).group(1)) 
                    cb  = float(det_s02_s20_c_re.match(l).group(2)) 
                    i_h = float(det_s02_s20_c_re.match(l).group(4))
                    i_p = float(det_s02_s20_c_re.match(l).group(7))
                    e_h = float(det_s02_s20_c_re.match(l).group(5))
                    e_p = float(det_s02_s20_c_re.match(l).group(8))
                    s_h = det_s02_s20_c_re.match(l).group(3)
                    s_p = det_s02_s20_c_re.match(l).group(6)
                    state.add_determinant(ca,cb,i_h,e_h,i_p,e_p,s_h,s_p)
                elif det_s01_s10_r_re.match (l) :
                    ca  = float(det_s01_s10_r_re.match(l).group(1))
                    cb  = 0.0
                    i   = float(det_s01_s10_r_re.match(l).group(3))
                    e   = float(det_s01_s10_r_re.match(l).group(4))
                    s   = det_s01_s10_r_re.match(l).group(2)
                    if current_sector == "01" :
                        i_p = i
                        e_p = e
                        s_p = s
                        i_h = 0
                        e_h = 0.0
                        s_h = ""
                    elif current_sector == "10" :
                        i_p = 0 
                        e_p = 0.0 
                        s_p = ""
                        i_h = i 
                        e_h = e 
                        s_h = s
                    state.add_determinant(ca,cb,i_h,e_h,i_p,e_p,s_h,s_p)
                elif det_s01_s10_c_re.match (l) :
                    ca  = float(det_s01_s10_c_re.match(l).group(1)) 
                    cb  = float(det_s01_s10_c_re.match(l).group(2)) 
                    i   = float(det_s01_s10_c_re.match(l).group(4))
                    e   = float(det_s01_s10_c_re.match(l).group(5))
                    s   = det_s01_s10_c_re.match(l).group(3)
                    if current_sector == "01" :
                        i_p = i
                        e_p = e
                        s_p = s
                        i_h = 0
                        e_h = 0.0
                        s_h = ""
                    elif current_sector == "10" :
                        i_p = 0 
                        e_p = 0.0 
                        s_p = ""
                        i_h = i 
                        e_h = e 
                        s_h = s
                    state.add_determinant(ca,cb,i_h,e_h,i_p,e_p,s_h,s_p)
        f.close()
# end

class tdrsp_results:
    """
    fscc_results fetched the relevant information from a dirac fock-space cc
    output
    """

    def __init__ (self) :
        self.states  = [ ]

    def process_output (self,outfile):
#
# what we want to read
#
        start_str         = "\s+Analysis of response solution vectors"
        tdrsp_symmetry    = "^.*solution vectors : PP EXCITATION([0-9a-zA-Z]+)\s+Irrep:\s*([0-9a-zA-Z]+)"
        tdrsp_energy      = "^ Freq.:\s+([0-9.-]+)\s+Norm:\s+([0-9.eEdD+-]+)\s+Residual norm:\s+([0-9.eEdD+-]+)"
        tdrsp_composition = "^\s+(\d+)\(i\:([0-9a-zA-Z]+)\)\s+--->\s+(\d+)\(v\:([0-9a-zA-Z]+)\)\s+([0-9.eEdD+-]+)"

        start_re          = re.compile (r''+ start_str   +'', re.IGNORECASE)
        tdrsp_sym_re      = re.compile (r''+ tdrsp_symmetry   +'', re.IGNORECASE)
        tdrsp_ener_re     = re.compile (r''+ tdrsp_energy   +'', re.IGNORECASE)
        tdrsp_compo_re    = re.compile (r''+ tdrsp_composition   +'', re.IGNORECASE)
#
# now we get to work
#
# for consistency, we define sector here, and since it is a single excitation we put it as 11=1h1p 
        current_sector = "11"

        print "reading time-dependend calculation (tddft,tdhf) output from file: ",outfile,"\n"
        f = file(outfile,'r')
        lines = f.readlines()
        read_state = False
        index      = 0
        label      = "a"
        label_prev = "a"
        for i, l in enumerate(lines) :
            if start_re.match (l) :
                read_state = True
            if read_state :
                ca   = 0.0
                cb   = 0.0
                e_h  = 0
                e_p  = 0
                s_h  = ""
                s_p  = ""
                if tdrsp_sym_re.match(l) :
                    label  = tdrsp_sym_re.match(l).group(2)
                elif tdrsp_ener_re.match(l) :
                    rel_en = float(tdrsp_ener_re.match(l).group(1))
                    norm   = float(tdrsp_ener_re.match(l).group(2))
                    if label != label_prev : 
                        index = 0
                        label_prev = label
                    index  = index + 1 

                    state  = electronic_state()
                    state.set_sector(current_sector)
                    state.set_label(label)
                    state.set_index(index)
                    state.set_energy(rel_en)
                    self.states.append(state)

                elif tdrsp_compo_re.match (l) :
#                   symm_h = tdrsp_compo_re.match(l).group(2)
#                   symm_p = tdrsp_compo_re.match(l).group(4)

                    weight = float(tdrsp_compo_re.match(l).group(5))
                    ca  = weight*math.sqrt(2)
                    cb  = 0.0
                    i_h = int(tdrsp_compo_re.match(l).group(1))
                    i_p = int(tdrsp_compo_re.match(l).group(3))
                    e_h = 0.0
                    e_p = 0.0
                    state.add_determinant(ca,cb,i_h,e_h,i_p,e_p,s_h,s_p)
        f.close()
# end
