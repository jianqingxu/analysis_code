#!/usr/bin/env python

## @file   evaluate_grafting.py
## @brief  analysis to evalute the antibody CDR grafting results
## @author Jianqing Xu
## @date   Last Modified: Apr.11th, 2013 by JQX

import optparse
import os, sys

from rosetta import *
init()


'''
	  CA   C   N
	 /  \ ? \ /
	C    N   CA
'''

_models_to_check_ = ["1.pdb", "2.pdb", "3.pdb", "4.pdb" ]


_grafting_list_ = ['L1', 'L2', 'L3', 'H1', 'H2'] # H3 is not on the list

_things_to_check_ = ['N-C', 'CA-N-C', 'N-C-CA', 'd_rmsd', 'd_score'  ]

_script_path_ = os.path.dirname( os.path.realpath(__file__) )

pose = Pose()



def load_models:
	pose_from_pdb(pose, file)

def check_N-C:

def check_CA-N:

def check_N-C-CA:

def check_d_rmsd:

def check_d_score:


def check_things:
	for thing in _things_to_check_:
		load_models()
		check_N-C()
		check_CA-N()
		check_N-C-CA()
		check_d_rmsd()
		check_d_score()

def output_results()


def main(args):
    parser = optparse.OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option('--profit',
			action="store", default='profit',
			help="Specify path+name for 'ProFIt' executable. Default is profit.",)
	
	parser.add_option('--benchmark_pdb_list',
			action="store", default=_script_path_+"Jianqing_List_for_RMSD_with_H3_Seq_NonRedudentH3",
			help="Specify path+name for the list of benchmark. Default is ./Jianqing_List_for_RMSD_with_H3_Seq_NonRedudentH3" ,)

    (options, args) = parser.parse_args(args=args[1:])
    global Options; Options = options

	
	check_things()
	
	output_results()


if __name__ == "__main__": 
    main(sys.argv)



