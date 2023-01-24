import argparse, os, sys
import numpy as np
# You can choose to write classes in other python files
# and import them here.
import MEMEoops as model

# This is the main function provided for you.
# Define additional functions to implement MEME
def main(args):
    # Parse input arguments
    # It is generally good practice to validate the input arguments, e.g.,
    # verify that the input and output filenames are provided
    seq_file_path = args.sequences_filename
    W = args.width
    model_file_path = args.model
    position_file_path = args.positions
    subseq_file_path = args.subseqs

    # Where you run your code.
    
    ### create a object of MEMEoops class + starting position pre
    o = model.MEMEoops(seq_file_path,
                     W,
                     output_paths = [model_file_path,position_file_path,subseq_file_path])
    o.enum_candidate()
    o.picking_starting_point()
    
    ### Starting
    stop_thres = 0.001
    em_progress = stop_thres + 1
    cur_pwm = o.starting_pwm
    #
    i = 0
    all_likelihood = []
    print("Starting EM iteration")
    while em_progress >= stop_thres:
        i += 1
        old_liklihood = o.cur_likelihood
        # E step
        o.E(cur_pwm)
        cur_z = o.cur_z
        # M step
        o.M(cur_z)
        cur_pwm = o.cur_pwm
        # 
        all_likelihood.append(o.cur_likelihood)
        new_liklihood = o.cur_likelihood
        em_progress = new_liklihood - old_liklihood
    #
    print("Finished.")
    print("Stop threshold: change in log likelihood smaller than ",stop_thres,")",sep = "")
    print("Totoal iter:", i)
    o.outputfile()
    ### END ###
        
# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module
if __name__ == "__main__":
    # Note: this example shows named command line arguments.  See the argparse
    # documentation for positional arguments and other examples.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sequences_filename',
                        help='sequences file path.',
                        type=str)
    parser.add_argument('--width',
                        help='width of the motif.',
                        type=int,
                        default=6)
    parser.add_argument('--model',
                        help='model output file path.',
                        type=str,
                        default='model.txt')
    parser.add_argument('--positions',
                        help='position output file path.',
                        type=str,
                        default='positions.txt')
    parser.add_argument('--subseqs',
                        help='subsequence output file path.',
                        type=str,
                        default='subseqs.txt')

    args = parser.parse_args()
    # Note: this simply calls the main function above, which we could have
    # given any name
    main(args)
