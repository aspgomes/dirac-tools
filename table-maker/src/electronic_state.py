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


from determinant import *

class electronic_state:
    """
    electronic state class:

    in composition, we'll add determiants (=elements of the determinant class)
    """

    def __init__ (self) :
        self.label        = "E"
        self.energy       = 0.0  
        self.relative_energy = 0.0  
        self.index        = 0.0
        self.fs_sector    = "00"
        self.composition  = [ ]
        self.max_index_p  = {}
        self.min_index_p  = {}
        self.max_index_h  = {}
        self.min_index_h  = {}
        self.symmetries   = []

    def set_sector (self,sector) :
        self.fs_sector = sector
        self.max_index_h[sector] = 0
        self.min_index_h[sector] = 10000
        self.max_index_p[sector] = 0
        self.min_index_p[sector] = 10000

    def get_sector (self) :
        return self.fs_sector

    def set_index (self,index) :
        self.index = index

    def get_index (self) :
        return self.index

    def set_relative_energy (self,relative_energy):
        self.relative_energy = relative_energy

    def set_energy (self,energy):
        self.energy = energy

    def get_relative_energy (self) :
        return self.relative_energy

    def get_energy (self) :
        return self.energy

    def set_label (self,label) :
        self.label = label

    def get_label (self) :
        return self.label

    def add_symmetry(self,new) :
        if new not in self.symmetries :
           self.symmetries.append(new)

    def add_determinant (self,ca,cb,index_h,spinor_h,index_p,spinor_p,symmetry_h,symmetry_p) :
        det = determinant()
        det.set_coef_r(ca)
        det.set_coef_i(cb)
        det.set_from_index(index_h)
        det.set_to_index(index_p)
        det.set_from_ener(spinor_h)
        det.set_to_ener(spinor_p)
        det.set_from_symmetry(symmetry_h)
        det.set_to_symmetry(symmetry_p)
        det.set_weight()
        if self.max_index_h[self.fs_sector] < int(index_h) :
            self.max_index_h[self.fs_sector] = int(index_h)
        if self.min_index_h[self.fs_sector] > int(index_h) :
            self.min_index_h[self.fs_sector] = int(index_h)
        if self.max_index_p[self.fs_sector] < int(index_p) :
            self.max_index_p[self.fs_sector] = int(index_p)
        if self.min_index_p[self.fs_sector] > int(index_p) :
            self.min_index_p[self.fs_sector] = int(index_p)
        self.composition.append(det)

        self.add_symmetry(symmetry_h)
        self.add_symmetry(symmetry_p)

    def print_min_max_indexes(self,sector) :
        print "maximum and minimum indexes for sector ",sector
        print "h: max ",self.max_index_h[sector],"   min ",self.min_index_h[sector]
        print "p: max ",self.max_index_p[sector],"   min ",self.min_index_p[sector]


    def setup_template_dets(self,sector):
        # setup templates for unique determinants that span all possible
        # index combinations for the model (P=Pm+Pi) spaces

        template_dets = []

        for sym1 in self.symmetries :
            if sector == "10" :
                for i_h in range(self.min_index_h[sector],(self.max_index_h[sector])+1) :
                    new_det = [ i_h, 0, 0, 0, 0, sym1, "" ]
                    template_dets.append(new_det)
            if sector == "01" :
                for i_p in range(self.min_index_p[sector],(self.max_index_p[sector])+1) :
                    new_det = [ 0, i_p, 0, 0, 0, "", sym1 ]
                    template_dets.append(new_det)
            if sector == "11" or sector == "02" or sector == "20" :
                for i_h in range(self.min_index_h[sector],(self.max_index_h[sector])+1) :
                    for sym2 in self.symmetries :
                        for i_p in range(self.min_index_p[sector],(self.max_index_p[sector])+1) :
                            new_det = [ i_h, i_p, 0, 0, 0, sym1, sym2 ]
                            template_dets.append(new_det)

        return template_dets


    def get_non_unique_dets (self,sector): 
        non_unique_dets = []
        #
        # fill the templates with the content of the determiants read from the output
        # and where we've accumulted the weights whenever we have the same indexes.
        #
        # todo: verify the boson symmetry of each index, to be sure that we're actually finding the same determinants
        #

        for d in range(len(self.composition)) :
            e_h = self.composition[d].get_from_ener()
            e_p = self.composition[d].get_to_ener()
            i_h = self.composition[d].get_from_index()
            i_p = self.composition[d].get_to_index()
            s_h = self.composition[d].get_from_symmetry()
            s_p = self.composition[d].get_to_symmetry()
            w   = self.composition[d].get_weight()
                 
            new_det = [ i_h, i_p, w, e_h, e_p, s_h, s_p ]
            non_unique_dets.append(new_det)

        return non_unique_dets


    def get_unique_dets (self,sector,group_by_energy=False) :
        unique_dets = self.setup_template_dets(sector)
        #
        # fill the templates with the content of the determiants read from the output
        # and where we've accumulted the weights whenever we have the same indexes.
        #
        # todo: verify the boson symmetry of each index, to be sure that we're actually finding the same determinants
        #
         
        for d in range(len(self.composition)) :
            e_h = self.composition[d].get_from_ener()
            e_p = self.composition[d].get_to_ener()
            i_h = self.composition[d].get_from_index()
            i_p = self.composition[d].get_to_index()
            s_h = self.composition[d].get_from_symmetry()
            s_p = self.composition[d].get_to_symmetry()
            w   = self.composition[d].get_weight()
            for p in range(len(unique_dets)) :
                ud = unique_dets[p]

                if group_by_energy:
                    if ud[3] == e_h and ud[4] == e_p :
                        ud[2]+= w
                        ud[3] = e_h
                        ud[4] = e_p
                    
                else :
                    if i_h == ud[0] and i_p == ud[1] and ud[5] == s_h and ud[6] == s_p :
                        ud[2]+= w
                        ud[3] = e_h
                        ud[4] = e_p

        return unique_dets


    def print_determinant(self,det) :
        i_h = det[0]
        i_p = det[1]
        w   = det[2]
        e_h = det[3]
        e_p = det[4]
        s_h = det[5]
        s_p = det[6]

#  add here code to translate the symmetry and index in symmetry to a global identifier
        i_global_p = -1
        i_global_h = -1

        if (i_h == 0) :
           print "   % 5.1f    % 3d (%3d %3s, % 6.4f) " % (w*100,i_global_p,i_p,s_p,e_p)
        elif (i_p == 0) :
           print "   % 5.1f    % 3d (%3d %3s, % 6.4f) " % (w*100,i_global_h,i_h,s_h,e_h)
        else :
           print "   % 5.1f    % 3d (%3d %3s, % 6.4f);  % 3d (%3d %3s, % 6.4f) " % (w*100,i_global_h,i_h,s_h,e_h,i_global_p,i_p,s_p,e_p)


    def print_list (self,sector,threshold,max_states,unique=True) :
        if self.index > max_states:
           return

        print "\n electronic state #",self.index," in symmetry ",self.label," energy: ",self.energy

        total_w = 0.0

        if not unique:
            dets = self.get_non_unique_dets(sector)
        else:
            dets = self.get_unique_dets(sector)

        for p in range(len(dets)) :
            d = dets[p]
            if d[2] >= threshold :
                self.print_determinant(d)
                total_w += d[2]

        print "    ----\n   % 5.1f\n" % (total_w*100)


    def print_table (self,sector,threshold,max_states,range_h=[],range_p=[],unique=True) :
        if self.index > max_states:
           return 

        print "\n electronic state #",self.index," in symmetry ",self.label," energy: ",self.energy

        total_w = 0.0

        dets = self.get_unique_dets(sector)

        if range_h == []:
            range_h = range(1,(self.max_index_h[sector])+1)  
        if range_p == []:
            range_p = range(1,(self.max_index_p[sector])+1)  

        if sector == "11" or sector == "02" or sector == "20" :
            relative_energy = self.get_relative_energy()
            print "\n printing electronic state composition in table format"
            print "% 8.0f  %6s  % 2d  |" % (relative_energy*219474.631280634, self.label, self.index),
            count = 35
            for i_h in range_h :
                for i_p in range_p :
                    for p in range(len(dets)) :
                        d = dets[p]
                        if i_h == d[0] and i_p == d[1] and s_h == d[5] and s_p == d[6] and d[2] >= threshold :
                            print " % 3d % 3d |" % (i_h,i_p),
                            count += 10
            print "\n","-"*count
            print "% 8.0f  %6s  % 2d  |" % (relative_energy*219474.631280634, self.label, self.index),
            for i_h in range_h :
                for i_p in range_p :
                    for p in range(len(dets)) :
                        d = dets[p]
                        if i_h == d[0] and i_p == d[1] and s_h == d[5] and s_p == d[6] and d[2] >= threshold :
                            print "  % 5.1f  |" % (d[2]*100),
            print "\n"


    def print_list_and_table(self,sector,threshold,max_states,range_h=[],range_p=[],unique=True) :
        if self.index > max_states:
           return
        print "\n electronic state #",self.index," in symmetry ",self.label," energy: ",self.energy
        total_w = 0.0

        if not unique:
            dets = self.get_non_unique_dets(sector)
        else:
            dets = self.get_unique_dets(sector)

        for p in range(len(dets)) :
            d = dets[p]
            if d[2] >= threshold :
                self.print_determinant(d)
                total_w += d[2]

        print "    ----\n   % 5.1f\n" % (total_w*100)

        #
        # additional step, print in a table format
        #
        if sector == "11" :
            relative_energy = self.get_relative_energy()
            print "% 8.0f  %6s  % 2d  |" % (relative_energy*219474.631280634, self.label, self.index),
            count = 30
            for i_h in range(1,(self.max_index_h[sector])+1) :
                for i_p in range(1,(self.max_index_p[sector])+1) :
                    for p in range(len(dets)) :
                        d = dets[p]
                        if i_h == d[0] and i_p == d[1] and d[2] >= threshold :
                            print " % 3d % 3d |" % (i_h,i_p),
                            count += 10
            print "\n","-"*count
            print "% 8.0f  %6s  % 2d  |" % (relative_energy*219474.631280634, self.label, self.index),
            for i_h in range(1,(self.max_index_h[sector])+1) :
                for i_p in range(1,(self.max_index_p[sector])+1) :
                    for p in range(len(dets)) :
                        d = dets[p]
                        if i_h == d[0] and i_p == d[1] and d[2] >= threshold :
                            print "  % 5.1f  |" % (d[2]*100),
            print "\n"

