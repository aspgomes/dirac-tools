#!/usr/bin/env python

#
# makes tables with electronic states' composistion, in terms of spinors 
# from fock-space coupled cluster or tddft outputs 
#

import reader

#
# program
# 

my_sector = "11"
res = reader.tdrsp_results()
myfile = "test_tddft.out"
#

res.process_output(myfile)
print "Excited states (sector ",my_sector,")\n"
for s in res.states :
    sector = s.get_sector()
    if sector == my_sector :
        s.print_list_and_table(sector,0.01,5)



