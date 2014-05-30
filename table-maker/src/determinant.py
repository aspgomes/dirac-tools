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


class determinant:
    """
    class to collect information on the determinants that make up a given electonic
    stat

    ca (cb) is the real (imaginary) part of the ci coefficient for the state
    from is the energy of the spinor the electron is excited from (sectors with e- removal)
    to   is the energy of the spinor  the electron is excited to  (sectors with e- addition)
    """

    def __init__ (self) :
        self.ca     = 0.0 
        self.cb     = 0.0 
        self.from_ener  = -1e+10 
        self.to_ener    =  1e+10 
        self.from_index = 0  
        self.to_index   = 0  
        self.weight     = 0
        self.from_symmetry = ""  
        self.to_symmetry   = ""

    def print_det (self):
        print "    Determinant:"
        if self.cb != 0.0 :
           print "\t c_a: % 6.4f,  c_b: % 6.4f (% 4.1f pc )"%(self.ca,self.cb,self.weight*100)
        else :
           print "\t c_a: % 6.4f (% 4.1f pc )"%(self.ca,self.weight*100)
        h_i = self.from_index
        p_i = self.to_index
        h_s = self.from_symmetry
        p_s = self.to_symmetry
        h_e = self.from_ener 
        p_e = self.to_ener
        if h_i != 0 :
           print "\t particle removed from spinor %4d (% 6.4f, %s) "%(h_i,h_e,h_s)
        if p_i != 0 :
           print "\t particle added   to   spinor %4d (% 6.4f, %s) "%(p_i,p_e,p_s)
        print " "

    def set_coef_r(self,c) :
        self.ca = float(c)
       
    def set_coef_i(self,c) :
        self.cb = float(c)
    
    def set_from_ener(self,e) :
        self.from_ener = e

    def set_to_ener(self,e) :
        self.to_ener = e
    
    def set_from_symmetry(self,s) :
        self.from_symmetry = s

    def set_to_symmetry(self,s) :
        self.to_symmetry = s 

    def set_from_index(self,i) :
        self.from_index = int(i)

    def set_to_index(self,i) :
        self.to_index = int(i)

    def set_weight(self) :
        self.weight = self.ca*self.ca + self.cb*self.cb

    def get_coef_r(self) :
        return self.ca
       
    def get_coef_i(self) :
        return self.cb
    
    def get_from_ener(self) :
        return self.from_ener

    def get_to_ener(self) :
        return self.to_ener
    
    def get_from_index(self) :
        return self.from_index

    def get_to_index(self) :
        return self.to_index

    def get_from_symmetry(self) :
        return self.from_symmetry

    def get_to_symmetry(self) :
        return self.to_symmetry

    def get_weight(self) :
        return self.weight 
