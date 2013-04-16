#!/usr/bin/env python

## @file   pdb_numb_renum.py
## @brief  use the pdb_info of one PDB, to renumber the 
##         PDBs having the same sequence
## @author Jianqing Xu
## @date   Last Modified: Apr.16th, 2013 by JQX
#

import sys
from rosetta import *
init()



_script_path_ = os.path.dirname( os.path.realpath(__file__) )



def get_uniq_chain_ids(pose):
    chains=[]
    for i in range(1, pose.total_residue()+1):
	chains.append( pose.pdb_info().chain(i) )

    from collections import OrderedDict
    return OrderedDict.fromkeys(chains).keys()


def main(args):

    parser = optparse.OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)


    parser.add_option('--pdb_list',
	action="store", default=_script_path_+"/list_to_change",
	help="Specify path+name for the list of benchmark. Default is ./list_to_change" ,)

    parser.add_option('--native',
	action="store", default=_script_path_+'/native',
	help="Specify path+name for the list of benchmark. Default is ./native.pdb" ,)

    (options, args) = parser.parse_args(args=args[1:])
    global Options; Options = options





    try:
	targets_list_file = open( Options.pdb_list, 'r' )
    except IOError as e: 
	print (  "({})".format(e)         )
	sys.exit()
    else:
	print "Targets List: "+Options.pdb_list
	targets_list_file.close()


    native_pose = Pose()
    pose_from_pdb(native_pose, Options.native)
    native_total_sequence = pose.sequence() # all the AAs in the native pose



    pose = Pose()

    # loop over all the mdoels that should be inspected
    for line in targets_list_file:  
	target_name = line.split()[0]
	print 
	print "****************** Working on Target " + target_name + "  *************************************" 

	pose.clear()
	pose_from_pdb(pose, _script_path_+target_name)

	model_uniq_chains==get_uniq_chain_ids(pose)

	for chain_id in model_uniq_chains:
	    seq = pose.chain_sequence(chain_id)
	    try:
		start_n = native_total_sequence.index(seq)
	    except ValueError as e:
		print (      "({})".format(e)           )
		print "Could not find the sequence in native pose"
		print "query_sequence:  "+seq
		print "native_sequence: "+native_total_sequence
		sys.exit()
	    else:
		end_n = start_n + len(seq) - 1 

	    truncated_pose = Pose( pose, start_id+1, end_id+1); # Rosetta style, starting from index 1
	    truncated_pose.dump_pdb(chain_id+".pdb")
	    truncated_pose.clear()


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



