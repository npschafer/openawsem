#!/usr/bin/env python3
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="AWSEM command line tool")
    subparsers = parser.add_subparsers(dest='subcommand', required=True, help='AWSEM subcommands')
       
    #parser.add_argument_group("Project", "AWSEM project management commands")
    subparsers.add_parser('create', help='Creates a new AWSEM project', add_help=False,)
    subparsers.add_parser('run', help='Runs the AWSEM project', add_help=False)
    subparsers.add_parser('analyze', help='Analyzes the AWSEM project', add_help=False)

    #parser.add_argument_group("Helper", "Useful functions for AWSEM")
    subparsers.add_parser('pdb2gro', help='Converts a pdb file to a gro file', add_help=False,)
    subparsers.add_parser('align_fragments', help='Aligns fragments to a fasta file', add_help=False)

    args, remaining_args = parser.parse_known_args()
    if len(remaining_args)==0:
        remaining_args = ['--help']
    
    original_argv = sys.argv[:] 
    sys.argv[0] = f"{original_argv[0]} {args.subcommand}"

    print(sys.argv)
    if args.subcommand == 'create':
        from openawsem.scripts import mm_create_project
        mm_create_project.main(remaining_args)
    elif args.subcommand == 'run':
        from openawsem.scripts import mm_run
        mm_run.main(remaining_args)
    elif args.subcommand == 'analyze':
        from openawsem.scripts import mm_analyze
        mm_analyze.main(remaining_args)
    elif args.subcommand == 'pdb2gro':
        from openawsem.helperFunctions import Pdb2Gro
        Pdb2Gro.main(remaining_args)
    elif args.subcommand == 'align_fragments':
        from openawsem.helperFunctions import align_fragments
        align_fragments.main(remaining_args)
    else:
        parser.print_help()

    sys.argv = original_argv

if __name__ == '__main__':
    main()
