#!/usr/bin/env python

## @file   pdb_numb_renum.py
## @brief  use the pdb_info of one PDB, to renumber the 
##         PDBs having the same sequence
## @author Jianqing Xu
## @date   Last Modified: Apr.16th, 2013 by JQX
#


from rosetta import *
init()



_script_path_ = os.path.dirname( os.path.realpath(__file__) )

_models_to_check_ = {'Grafted'  : "/output/details/1.pdb", 
		     'Grafted_Min' : "/output/details/2.pdb",
		     'Grafted_Min_Opt' : "/output/details/3.pdb",
		     'Grafted_Min_Opt_Min': "/output/details/4.pdb",
		     'Grafted_Min_Opt_Min_Relaxed': "/output/grafted.relaxed.pdb",
		    }

def get_uniq_chain_ids(pose):

    uniq_chain_ids=[ pose.pdb_info().chain(1) ] # the 1st residue chain ID

    for i in range(2, pose.total_residue()+1):  # loop the rest residue's chain ID
	count = 0
	for ch_id in uniq_chain_ids:
	    if pose.pdb_info().chain(i) == uniq_chain_ids[ch_id]:
		count+=1;
	if count ==0:
	    uniq_chain_ids.append( pose.pdb_info().chain(i) )

    return uniq_chain_ids


def main(args):

    parser = optparse.OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)


    parser.add_option('--pdb_list',
	action="store", default=_script_path_+"/list_to_change",
	help="Specify path+name for the list of benchmark. Default is ./list_to_change" ,)

    parser.add_option('--native',
	action="store", default=_script_path_+'native',
	help="Specify path+name for the list of benchmark. Default is ./native.pdb" ,)

    (options, args) = parser.parse_args(args=args[1:])
    global Options; Options = options


    pose = Pose()

    native_pose = Pose()
    pose_from_pdb(native_pose, Options.native)


    uniq_chain_ids=get_uniq_chain_ids(native_pose)

    print "Targets List: "+Options.pdb_list
    targets_list_file = open( Options.pdb_list, 'r' )

    # loop over all the 53 PDB targets in the Benchmark
    for line in targets_list_file:  
	target_name = line.split()[0]
	print 
	print "****************** Working on Target " + target_name + "  *************************************" 

	n_chains=get_uniq_chain_ids(target_name)
	if uniq_chain_ids !=


	# loop over all the mdoels that should be inspected


	for model in _models_to_check_:  
	    file_name = _script_path_ + "/" + target_name + _models_to_check_[model]
	    print "          " + file_name
	    pose.clear(); 
	    pose_from_pdb(pose, file_name)

	    #align to the native structure
	    align_to_native_by_framework(pose, native_pose)

	    #do all the measurement
	    check_things(pose,native_pose, scorefxn)
	
	    #output_results()


if __name__ == "__main__": 
    main(sys.argv)



