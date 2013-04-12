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
     N    C    CA  
      \  /  ? /  \ 
       CA    N    C

_grafting_CDRs_ = {'L1' : {'ch_id':'L', 'nter': 24, 'cter': 34}, 
		   'L2' : {'ch_id':'L', 'nter': 50, 'cter': 56}, 
		   'L3' : {'ch_id':'L', 'nter': 89, 'cter': 97}, 
		   'H1' : {'ch_id':'H', 'nter': 26, 'cter': 35}, 
		   'H2' : {'ch_id':'H', 'nter': 50, 'cter': 65} 
		   } # H3 is not on the list
'''




_grafting_stems_ = {'L1_stem' : {'ch_id':'L', 
                                 'nter' : [20, 21, 22, 23], 
                                 'cter' : [35, 36, 37, 38]
			        }, 
		    'L2_stem' : {'ch_id':'L', 
			         'nter' : [46, 47, 48, 49], 
			         'cter' : [57, 58, 59, 60] 
			        }, 
		    'L3_stem' : {'ch_id':'L', 
				 'nter' : [85, 86, 87, 88], 
				 'cter':  [98, 99, 100, 101]
				}, 
		    'H1_stem' : {'ch_id':'H', 
				 'nter' : [22, 23, 24, 25], 
				 'cter' : [36, 37, 38, 39]
				}, 
		    'H2_stem' : {'ch_id':'H', 
				  'nter': [46, 47, 48, 49], 
				  'cter': [66, 67, 68, 69] 
				} 
		   } # H3 is not on the list

_script_path_ = os.path.dirname( os.path.realpath(__file__) )

_models_to_check_ = {'Grafted'  : "/output/details/1.pdb", 
		     'Grafted_Min' : "/output/details/2.pdb",
		     'Grafted_Min_Opt' : "/output/details/3.pdb",
		     'Grafted_Min_Opt_Min': "/output/details/4.pdb",
		     'Grafted_Min_Opt_Min_Relaxed': "/output/grafted.relaxed.pdb",
		     'Crystal': ".renum.best_packed.pdb"
		    }



#def check_C_N_bond():
#def check_CA_C_N_angle():
#def check_C_N_CA_angle():
#def check_d_rmsd():
#def check_d_score():


#def check_things(pose):
	#check_C_N_bond(pose)
	#check_CA_C_N()
	#check_C_N_CA()
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
	    # Four Stem Residues A-B-C-D
	    for stems in _grafting_stems_:
		print _grafting_stems_[stems]['nter']

		ch_id    = _grafting_stems_[stems]['ch_id']

		B_pose_num = pose.pdb_info().pdb2pose(ch_id, _grafting_stems_[stems]['nter'][1])
		C_pose_num = pose.pdb_info().pdb2pose(ch_id, _grafting_stems_[stems]['nter'][2])
		#check_things(pose)
	
		#output_results()


if __name__ == "__main__": 
    main(sys.argv)



