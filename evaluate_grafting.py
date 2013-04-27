#!/usr/bin/env python

## @file   evaluate_grafting.py
## @brief  analysis to evalute the antibody CDR grafting Results
## @author Jianqing Xu
## @date   Last Modified: Apr.11th, 2013 by JQX
#

import optparse
import os, sys
import json
import commands
import copy
import math

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
				 "cter" : { "pdb_nums" : [35, 36, 37, 38],
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


_scorefxn_=create_score_function_ws_patch("standard", "score12")



_models_to_check_ ={
			"Grafted"  : "/output/details/1.pdb",
			"Grafted_Min" : "/output/details/2.pdb",
			"Grafted_Min_Opt" : "/output/details/3.pdb",
			"Grafted_Min_Opt_Min": "/output/details/4.pdb",
			"Grafted_Min_Opt_Min_Relaxed": "/output/grafted.relaxed.pdb",
		}


_profit_templates_ ={
"L1_stem":'''reference %s
mobile %s 
ATOM C,CA,N,O 
ZONE L16-L19 
ZONE L39-L42 
fit 
write %s
quit''',
"L2_stem":'''reference %s
mobile %s 
ATOM C,CA,N,O 
ZONE L42-L45 
ZONE L61-L64
fit 
write %s
quit''',
"L3_stem":'''reference %s
mobile %s 
ATOM C,CA,N,O 
ZONE L81-L84 
ZONE L102-L105
fit 
write %s
quit''',
"H1_stem":'''reference %s
mobile %s 
ATOM C,CA,N,O 
ZONE H18-H21 
ZONE H40-H43
fit 
write %s
quit''',
"H2_stem":'''reference %s
mobile %s 
ATOM C,CA,N,O 
ZONE H42-H45 
ZONE H70-H73
fit 
write %s
quit''',
}


def align_to_native_by_framework(pose, native_pose, stems):

	native_path_name = native_pose.pdb_info().name()
	model_path_name  = pose.pdb_info().name()
	
	native_name = native_path_name[native_path_name.rfind("/")+1: ]
	model_name  =  model_path_name[model_path_name.rfind("/")+1: ]

	native_path = native_path_name[0: native_path_name.rfind("/")+1 ]
	model_path  =  model_path_name[0: model_path_name.rfind("/")+1 ]

	analysis_path = native_path+"analysis/"
	if not os.path.exists( analysis_path ):
		os.makedirs(analysis_path)
	
	profit_in_file  =  analysis_path + "profit_" + stems + "_" + model_name[:-4] + ".in"
	profit_out_file =  analysis_path + "profit_" + stems + "_" + model_name[:-4] + ".out"

	aligned_model_path_name = analysis_path + model_name[:-4] + "_" + stems + "_fitted.pdb"

	with file(profit_in_file, 'w') as f:
		f.write(   _profit_templates_[stems] % (native_path_name, model_path_name, aligned_model_path_name)   )


	commandline = '%s < %s > %s' % (Options.profit, profit_in_file, profit_out_file)
	status, output = commands.getstatusoutput(commandline);

	aligned_pose = Pose()
	pose_from_pdb(aligned_pose, aligned_model_path_name); _scorefxn_(aligned_pose)

	return aligned_pose


# checking 1.pdb or 2.pdb or 3.pdb or 4.pdb in pose
def check_things(pose, native_pose):

	gf_st =  copy.deepcopy(_grafting_stems_)

	# loop over L1_stem, L2_stem, L3_stem, H1_stem, H2_stem
	# Four Stem Residues A-B-|-C-D
	for stems in gf_st:


		aligned_pose = align_to_native_by_framework(pose, native_pose, stems)
		print "                            &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
		print "                            &&&&&&&&&&&&&&&&&&&&&&&&&&&  Finish Reading ....  ", aligned_pose.pdb_info().name()

		# loop over "nter" and "cter"
		terminus=["nter", "cter"]
		for ter in terminus:

			ch_id = gf_st[stems]["ch_id"]

			# define    B-|-C
			A_pose_num = pose.pdb_info().pdb2pose(ch_id, gf_st[stems][ter]["pdb_nums"][0])
			B_pose_num = pose.pdb_info().pdb2pose(ch_id, gf_st[stems][ter]["pdb_nums"][1])
			C_pose_num = pose.pdb_info().pdb2pose(ch_id, gf_st[stems][ter]["pdb_nums"][2])
			D_pose_num = pose.pdb_info().pdb2pose(ch_id, gf_st[stems][ter]["pdb_nums"][3])


			B_n  = AtomID(1, B_pose_num)  # in Rosetta, index starts from 1
			B_ca = AtomID(2, B_pose_num)
			B_c  = AtomID(3, B_pose_num)

			C_n  = AtomID(1, C_pose_num)
			C_ca = AtomID(2, C_pose_num)
			C_c  = AtomID(3, C_pose_num)

			#check_C_N_bond
			gf_st[stems][ter]["C_N_bond"] = pose.conformation().bond_length(B_c, C_n)

			#check_CA_C_N()
			gf_st[stems][ter]["CA_C_N_angle"] = math.degrees( pose.conformation().bond_angle(B_ca, B_c, C_n) )

			#check_C_N_CA()
			gf_st[stems][ter]["C_N_CA_angle"] = math.degrees( pose.conformation().bond_angle(B_c, C_n, C_ca) )

			#check_d_rmsd()
			loop = Loop( A_pose_num, D_pose_num, B_pose_num )
			loops= Loops(); loops.add_loop(loop)
			rms = loop_rmsd(aligned_pose, native_pose, loops, True)
			gf_st[stems][ter]["d_rmsd"] = rms

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





_results_={
	    "C_N_bond"     : [],
            "CA_C_N_angle" : [],
            "C_N_CA_angle" : [],
            "d_rmsd"       : [],
            "d_score"      : []
            }


def get_all_results(All_Targets):

    models_results = {}

    # loop over l1.pdb, l2.pdb, l3.pdb, l4.pdb ...
    for model in sorted(_models_to_check_):

	results = copy.deepcopy(  _results_  )

	#loop over "C_N_bond", "CA_C_N_angle", "C_N_CA_angle", "d_rmsd", "d_score"
	for result in results:

	    for target in All_Targets:
		for stem in All_Targets[target][model]:
		    for ter in ["nter", "cter"]:
			results[result].append  ( All_Targets[target][model][stem][ter][result] )

	models_results[model] = copy.deepcopy(results)

    for model in sorted(models_results):
	print model
	print models_results[model]
	print ""
	print ""

    return models_results


def histogram(data, bin_size):

    max_bound = float(max(data))
    min_bound = float(min(data))
    num_bins = int(  (max_bound-min_bound)/float(bin_size)  ) + 1

    hist = [0.0 for x in range(num_bins)]

    for num in data:
	bin_location=int(  (float(num) - float(min_bound))/float(bin_size)  ) 
	hist[ bin_location ] += 1.0

    x_values = [0.0 for x in range(num_bins)] 
    for i in range(0, num_bins):
	hist[i] /= float(len(data)) 
	x_values[i] = min_bound + float(i)*float(bin_size)+0.5*float(bin_size)



    histogram_xy = {"x_values": copy.deepcopy(x_values), 
		    "y_values": copy.deepcopy(hist)  }

    return histogram_xy



_bin_creteria_={ "C_N_bond"     : 0.01 ,
                 "CA_C_N_angle" : 0.1 ,
                 "C_N_CA_angle" : 0.1 ,
                 "d_rmsd"       : 0.1 ,
                 "d_score"      : 0.1 
	}

def calculate_distribution(models_results):
    models_distribution = { }

    # loop over 1.pdb, 2.pdb, 3.pdb ...
    for model in sorted(models_results):
	#print model
	distribution = {}
	for result in sorted(models_results[model]):
	    #print result 
	    bin_size  = copy.deepcopy(  _bin_creteria_[result]  )
	    distribution[result] = copy.deepcopy( histogram(models_results[model][result], bin_size) )
	    #print models_results[model][result]
	    #print 
	    #print

	models_distribution[model] = copy.deepcopy(distribution)


    return models_distribution



def output_final_distribution_results( models_distribution, native_pose):
    native_path_name = native_pose.pdb_info().name()
    native_path = native_path_name[0: native_path_name.rfind("/") ]
    path = native_path[0: native_path.rfind("/") +1 ] + "analysis/"

    if not os.path.exists( path ):
	os.makedirs(path)



    for model in sorted(models_distribution):
	for result in sorted(models_distribution[model]):
	    fname = path+model+"_"+result+"_distribution"
	    f = open(fname,'w')

	    tot_num_x = len(models_distribution[model][result]["x_values"])
	    tot_num_y = len(models_distribution[model][result]["y_values"]) 
	    if tot_num_x != tot_num_y:
		print "The length in X and Y dismatches !!!!"
		sys.exit()

	    x_less = str(  models_distribution[model][result]["x_values"][0]- _bin_creteria_[result]  )
	    y_less = str(0.0)
	    f.write('%s  %s\n' % (x_less, y_less) )

	    for i in range( 0, tot_num_x  ) :
		str1=str(   str(models_distribution[model][result]["x_values"][i])  )
		str2=str(   str(models_distribution[model][result]["y_values"][i])  )
		f.write('%s  %s\n' % (str1, str2) )


	    x_more = str(  models_distribution[model][result]["x_values"][tot_num_x-1]+ _bin_creteria_[result]  )
	    y_more = str(0.0)
	    f.write('%s  %s\n' % (x_more, y_more) )

	    f.close()

	    
	    
	    
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
    All_Targets={}
    for line in targets_list_file:
		target_name = line.split()[0]
		print 
		print "********************************************************************************"
		print "*                                                                              *"
		print "*                  Working on Target " + target_name + "                                      *"
		print "*                                                                              *"
		print "********************************************************************************"
		native_file_path = _script_path_ + "/" + target_name + "/"
		native_file_name = native_file_path + target_name + ".renum.best_packed.pdb"

		native_pose.clear(); pose_from_pdb(native_pose, native_file_name); _scorefxn_(native_pose)

		# loop over all the mdoels that should be inspected
		# like 1.pdb, 2.pdb, 3.pdb, 4.pdb, and grafted.relaxed.pdb
		values = {}
		for model in sorted(_models_to_check_) :   # the sort here is to sort the order of keys
		
			# read the model 1.pdb, 2.pdb, 3.pdb, 4.pdb and grafted.relaxed.pdb
			file_name = _script_path_ + "/" + target_name + _models_to_check_[model]
			print 
			print "          ####################################################################################################"
			print "          ###############" + file_name
			pose.clear(); pose_from_pdb(pose, file_name); _scorefxn_(pose)

			#do all the measurements
			values[model] = copy.deepcopy(   check_things(pose, native_pose)    )

		All_Targets[target_name] = copy.deepcopy( values )

    models_results =   get_all_results(All_Targets)  

    models_distribution = calculate_distribution(models_results)

    output_final_distribution_results( models_distribution , native_pose )



if __name__ == "__main__": 
    main(sys.argv)



