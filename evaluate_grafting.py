#!/usr/bin/env python

## @file   evaluate_grafting.py
## @brief  analysis to evalute the antibody CDR grafting results
## @author Jianqing Xu
## @date   Last Modified: Apr.11th, 2013 by JQX
#

import optparse
import os, sys
import json

from rosetta import *
init()


'''
     A-B-|-C-D  4 residues for a certain stem 

     B and C connecting in the following way

     ....N    C    CA  
          \  /  ? /  \ 
           CA    N    C....

_grafting_CDRs_ = {
               "L1" : {"ch_id":"L", 
                        "nter": 24, 
                        "cter": 34
                      }, 
		       "L2" : {"ch_id":"L", 
		                "nter": 50, 
		                "cter": 56
		              }, 
		       "L3" : {"ch_id":"L", 
		                "nter": 89, 
		                "cter": 97
		              }, 
		       "H1" : {"ch_id":"H", 
		                "nter": 26, 
		                "cter": 35
		              }, 
		       "H2" : {"ch_id":"H", 
		                "nter": 50, 
		                "cter": 65
		              } 
		   } # H3 is not on the list
'''



_grafting_stems_ ={
			"L1_stem" : {"ch_id":"L",
				 "nter" : { "pdb_nums" : [20, 21, 22, 23], 
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					  },
				 "cter" : { "pdb_num" : [35, 36, 37, 38],
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					  }
			        }, 
		    "L2_stem" : {"ch_id":"L", 
				 "nter" : { "pdb_nums" : [46, 47, 48, 49], 
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					    },
				 "cter" : { "pdb_nums" : [57, 58, 59, 60],
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					  }
			        }, 
		    "L3_stem" : {"ch_id":"L", 
				 "nter" : { "pdb_nums" : [85, 86, 87, 88], 
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					    },
				 "cter":  { "pdb_nums" : [98, 99, 100, 101],
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					    }
				}, 
		    "H1_stem" : {"ch_id":"H", 
				 "nter" : { "pdb_nums" : [22, 23, 24, 25], 
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					    },
				 "cter" : { "pdb_nums" : [36, 37, 38, 39],
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					    }
				}, 
		    "H2_stem" : {"ch_id":"H", 
				  "nter": { "pdb_nums" : [46, 47, 48, 49], 
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					    },
				  "cter": { "pdb_nums" :  [66, 67, 68, 69],
				            "C_N_bond" : 0,
					    "CA_C_N_angle" : 0,
					    "C_N_CA_angle" : 0,
					    "d_rmsd" : 0,
					    "d_score" :0
					    }
				} 
		   } # H3 is not on the list




_script_path_ = os.path.dirname( os.path.realpath(__file__) )


_scorefxn_=create_score_function_ws_patch('standard', 'score12')



_models_to_check_ ={
			'Grafted'  : "/output/details/1.pdb",
			'Grafted_Min' : "/output/details/2.pdb",
			'Grafted_Min_Opt' : "/output/details/3.pdb",
			'Grafted_Min_Opt_Min': "/output/details/4.pdb",
			'Grafted_Min_Opt_Min_Relaxed': "/output/grafted.relaxed.pdb",
		}



_profit_templates_ ={
"L":'''reference %s
mobile %s 
ATOM C,CA,N,O 
ZONE L5-L23 
ZONE L35-L49 
ZONE L57-L88 
ZONE L98-L100 
fit 
write %s.L_fitted.pdb 
quit''',
"H":'''reference %s 
mobile %s 
ATOM C,CA,N,O 
ZONE H5-H25 
ZONE H36-H49 
ZONE H66-H94 
ZONE H103-H105 
fit 
write %s.H_fitted.pdb
quit''' 
}


def align_to_native_by_framework(pose, native_pose, chain):

	native_path_name = native_pose.pdb_info().name()
	model_path_name  = pose.pdb_info().name()
	
	native_name = native_path_name[native_path_name.rfind("/")+1: ]
	model_name  =  model_path_name[model_path_name.rfind("/")+1: ]
	

	native_path = native_path_name[0: native_path_name.rfind("/")+1 ]
	model_path  =  model_path_name[0: model_path_name.rfind("/")+1 ]
	
	profit_in_file  =  native_path + "profit" + chain + ".in"
	profit_out_file =  native_path + "profit" + chain + ".out"

	with file(profit_in_file, 'w') as f:
		f.write(_profit_templates_[chain] % (native_path_name, model_path_name, native_path+model_name[:-4]) )


	#f_name = '\ '.join(f_name.split())
	sys.exit()
	'''
	commandline = '%s < %s.in > %s.out' % (Options.profit, f_name, f_name)
	res, output = commands.getstatusoutput(commandline);
	if res: print commandline, output; sys.exit(1)

	pathPrefix = '\ '.join(vars()['prefix'].split())
        pathPrefix = pathPrefix[:-1] if pathPrefix.endswith('/') else pathPrefix
	res, output = commands.getstatusoutput('cat {0}/fitted.L.pdb {0}/fitted.H.pdb > {0}/FR.pdb'.format(pathPrefix))

        #res, output = commands.getstatusoutput('cat %(prefix)s/fitted.L.pdb %(prefix)s/fitted.H.pdb > %(prefix)s/FR.pdb' % vars())
	if res: print output;  sys.exit(1)
	'''



# checking 1.pdb or 2.pdb or 3.pdb or 4.pdb in pose
def check_things(pose, native_pose):

	gf_st =  _grafting_stems_
	terminus=['nter', 'cter']

	# loop over L1, L2, L3, H1, H2
    # Four Stem Residues A-B-|-C-D
	for stems in gf_st:

		# loop over 'nter' and 'cter'
		for ter in terminus:

			ch_id = gf_st[stems]['ch_id']

			# define    B-|-C
			A_pose_num = pose.pdb_info().pdb2pose(ch_id, gf_st[stems][ter]["pdb_nums"][0])
			B_pose_num = pose.pdb_info().pdb2pose(ch_id, gf_st[stems][ter]["pdb_nums"][1])
			C_pose_num = pose.pdb_info().pdb2pose(ch_id, gf_st[stems][ter]["pdb_nums"][2])
			D_pose_num = pose.pdb_info().pdb2pose(ch_id, gf_st[stems][ter]["pdb_nums"][3])

			B_N  = AtomID(1, B_pose_num)  # in Rosetta, index starts from 1
			B_CA = AtomID(2, B_pose_num)
			B_C  = AtomID(3, B_pose_num)

			C_N  = AtomID(1, C_pose_num)
			C_CA = AtomID(2, C_pose_num)
			C_C  = AtomID(3, C_pose_num)

			#check_C_N_bond
			gf_st[stems][ter]["C_N_bond"] = pose.conformation().bond_length(B_C, C_N)

			#check_CA_C_N()
			gf_st[stems][ter]["CA_C_N_angle"] = pose.conformation().bond_angle(B_CA, B_C, C_N)

			#check_C_N_CA()
			gf_st[stems][ter]["C_N_CA_angle"] = pose.conformation().bond_angle(B_C, C_N, C_CA)

			#check_d_rmsd()
			print stems[0]
			if stems[0]== "H":
				align_to_native_by_framework(pose, native_pose, "H")
			if stems[0]== "L":
				align_to_native_by_framework(pose, native_pose, "L")


			#check_d_score()
			model_stem_score=0.0
			model_stem_score+=pose.energies().residue_total_energy(A_pose_num)
			model_stem_score+=pose.energies().residue_total_energy(B_pose_num)
			model_stem_score+=pose.energies().residue_total_energy(C_pose_num)
			model_stem_score+=pose.energies().residue_total_energy(D_pose_num)

			native_stem_score=0.0
			native_stem_score+=native_pose.energies().residue_total_energy(A_pose_num)
			native_stem_score+=native_pose.energies().residue_total_energy(B_pose_num)
			native_stem_score+=native_pose.energies().residue_total_energy(C_pose_num)
			native_stem_score+=native_pose.energies().residue_total_energy(D_pose_num)

			gf_st[stems][ter]["d_score"]=model_stem_score-native_stem_score

	return gf_st


#def output_results(results):


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


    # initialize the pose and scoring function in Rosetta
    pose = Pose()
    native_pose = Pose()

    print "Targets List: "+Options.benchmark_pdb_list
    targets_list_file = open( Options.benchmark_pdb_list, 'r' )

    # loop over all the 53 PDB targets in the Benchmark
    for line in targets_list_file:
		target_name = line.split()[0]
		print 
		print "****************** Working on Target " + target_name + "  *************************************"
		native_file_path = _script_path_ + "/" + target_name + "/"
		native_file_name = native_file_path + target_name + ".renum.best_packed.pdb"

		native_pose.clear(); pose_from_pdb(native_pose, native_file_name); _scorefxn_(native_pose)

		# loop over all the mdoels that should be inspected
		# like 1.pdb, 2.pdb, 3.pdb, 4.pdb, and grafted.relaxed.pdb
		results=[]
		for model in _models_to_check_ :
		
			# read the model 1.pdb, 2.pdb, 3.pdb, 4.pdb and grafted.relaxed.pdb
			file_name = _script_path_ + "/" + target_name + _models_to_check_[model]
			print "          ###############" + file_name
			pose.clear(); pose_from_pdb(pose, file_name); _scorefxn_(pose)

	    	#do all the measurements
			gf_st =  check_things(pose, native_pose)

			results.append(gf_st)


		#output_results(results)






if __name__ == "__main__": 
    main(sys.argv)



