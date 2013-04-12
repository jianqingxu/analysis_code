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

_grafting_CDRs_ = {'L1' : {'ch_id':'L', 'nter': 24, 'cter': 34}, 
		   'L2' : {'ch_id':'L', 'nter': 50, 'cter': 56}, 
		   'L3' : {'ch_id':'L', 'nter': 89, 'cter': 97}, 
		   'H1' : {'ch_id':'H', 'nter': 26, 'cter': 35}, 
		   'H2' : {'ch_id':'H', 'nter': 50, 'cter': 65} 
		   } # H3 is not on the list
'''




_grafting_stems_ = {'L1_stem' : {'ch_id':'L', 
                                 'nter' : {'start':20, 'end':23}, 
                                 'cter' : {'start':35, 'end':38}
			        }, 
		    'L2_stem' : {'ch_id':'L', 
			         'nter' : {'start':46, 'end':49}, 
			         'cter' : {'start':57, 'end':60} 
			        }, 
		    'L3_stem' : {'ch_id':'L', 
				 'nter' : {'start':85, 'end':88}, 
				 'cter':  {'start':98, 'end':101}
				}, 
		    'H1_stem' : {'ch_id':'H', 
				 'nter' : {'start':22, 'end':25}, 
				 'cter' : {'start':36, 'end':39}
				}, 
		    'H2_stem' : {'ch_id':'H', 
				  'nter': {'start':46, 'end':49}, 
				  'cter': {'start':66, 'end':69} 
				} 
		   } # H3 is not on the list
_stem_size_ = 4

_script_path_ = os.path.dirname( os.path.realpath(__file__) )

_models_to_check_ = {'Grafted'  : "/output/details/1.pdb", 
		     'Grafted_Min' : "/output/details/2.pdb",
		     'Grafted_Min_Opt' : "/output/details/3.pdb",
		     'Grafted_Min_Opt_Min': "/output/details/4.pdb",
		     'Grafted_Min_Opt_Min_Relaxed': "/output/grafted.relaxed.pdb",
		     'Crystal': ".renum.best_packed.pdb"
		    }




#def check_N_C_bond():
#def check_CA_N_C():
#def check_N_C_CA():
#def check_d_rmsd():
#def check_d_score():


#def check_things(pose):
	#check_N_C(pose)
	#check_CA_N_C()
	#check_N_C_CA()
	#check_d_rmsd()
	#check_d_score()

#def output_results():


def main(args):

    parser = optparse.OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option('--profit',
	action="store", default='profit',
	help="Specify path+name for 'ProFIt' executable. Default is profit.",)

    parser.add_option('--benchmark_pdb_list',
	action="store", default=_script_path_+"/Jianqing_List_for_RMSD_with_H3_Seq_NonRedudentH3",
	help="Specify path+name for the list of benchmark. Default is ./Jianqing_List_for_RMSD_with_H3_Seq_NonRedudentH3" ,)

    (options, args) = parser.parse_args(args=args[1:])
    global Options; Options = options


    pose = Pose()

    print "Targets List: "+Options.benchmark_pdb_list
    targets_list_file = open( Options.benchmark_pdb_list, 'r' )

    # loop over all the 53 PDB targets in the Benchmark
    for line in targets_list_file:  
	target_name = line.split()[0]
	print 
	print "****************** Working on Target " + target_name + "  *************************************" 

	# loop over all the mdoels that should be inspected
	for model in _models_to_check_:  
	    if model ==  "Crystal":
		file_name = _script_path_ + "/" + target_name + "/" + target_name + _models_to_check_[model]
	    else:
		file_name = _script_path_ + "/" + target_name + _models_to_check_[model]

	    print "          " + file_name
	    pose.clear(); pose_from_pdb(pose, file_name)

	    # loop over L1, L2, L3, H1, H2
	    for graft in _grafting_stems_:
		print _grafting_stems_[graft]['nter']

		ch_id   = _grafting_stems_[graft]['ch_id']
		for _grafting_stems_[graft]['nter']['start']

		pose_num = pose.pdb_info().pdb2pose(
		#check_things(pose)
	
		#output_results()


if __name__ == "__main__": 
    main(sys.argv)



