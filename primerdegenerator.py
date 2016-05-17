#!/usr/bin/env python
from SPAdesPipeline.OLCspades.accessoryFunctions import *

__author__ = 'adamkoziol'


class DegeneratePrimers(object):
    def objectifier(self):
        import primerobject
        # Initialise the primer object
        self.runmetadata = primerobject.PrimerObject(self)
        if self.batch:
            printtime('Performing batch analyses', self.starttime)
            self.runmetadata.batch()
        else:
            printtime('Performing single analysis', self.starttime)
            self.runmetadata.single()

    def degenerate(self):
        """
        Creates all the possible primers from sequences with ambiguous characters
        """
        import itertools
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        from Bio.Seq import Seq
        # Create the output file path/name
        filename = '{}unambiguousprimers.fasta'.format(self.path)
        with open(filename, 'wb') as unambiguous:
            for sample in self.runmetadata.samples:
                # Create a list of all possible forward and reverse primers
                ambiguouslist = self.extend_ambiguous_dna(sample.forward, sample.reverse)
                # Create a list of tuples of all possible primer pairs
                ambiguouspairs = list(itertools.product(*ambiguouslist))
                for iterator, pair in enumerate(ambiguouspairs):
                    # Zip together the forward and reverse primers in :pair with 'F' and 'R', respectively
                    for item in zip(pair, ['F', 'R']):
                        # The definition line will be something like: stx1a_F_0
                        definitionline = '{}_{}_{}'.format(sample.name, item[1], iterator)
                        # Convert the string of the sequence to a Seq
                        sequence = Seq(item[0])
                        # Create a sequence record using BioPython
                        fasta = SeqRecord(sequence,
                                          # Without this, the header will be improperly formatted
                                          description='',
                                          # Use >:definitionline as the header
                                          id=definitionline)
                        # Use the SeqIO module to properly format the new sequence record
                        SeqIO.write(fasta, unambiguous, "fasta")

    @staticmethod
    def extend_ambiguous_dna(forward, reverse):
        """
        Return list of all possible sequences given an ambiguous DNA input
        from: https://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence
        """
        from Bio import Seq
        from itertools import product
        d = Seq.IUPAC.IUPACData.ambiguous_dna_values
        return [list(map("".join, product(*map(d.get, forward)))), list(map("".join, product(*map(d.get, reverse))))]

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args:
        :param pipelinecommit:
        :param startingtime:
        :param scriptpath:
        """
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.starttime = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        self.args = args
        # Forward and reverse primers (if supplied)
        self.forward = args.forwardprimer
        self.reverse = args.reverseprimer
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Output location is not a valid directory {0!r:s}'.format(self.path)
        self.batch = False
        self.primerfile = ''
        # If primer sequences are supplied, don't look for a primer file. If one or more of the primers are missing,
        # then attempt to try the primer file
        if not self.forward or not self.reverse:
            print 'One or more of the primers was not provided. Attempting to find and use a supplied primer file.'
            self.primerfile = args.primerfile
            assert os.path.isfile(self.primerfile), u'Cannot file primer file {0!r:s}'.format(self.primerfile)
            self.batch = True
        # If both primers are provided, set them to uppercase
        if self.forward and self.reverse:
            self.forward = self.forward.upper()
            self.reverse = self.reverse.upper()
        #
        self.runmetadata = MetadataObject()
        # Create the objects
        self.objectifier()
        # Run the degeneration
        self.degenerate()


if __name__ == '__main__':
    import subprocess
    import time
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser

    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Parser for arguments
    parser = ArgumentParser(description='Output all possible primer pairs from primers with degenerate bases.'
                                        'You must provide either a file of all the primer pairs for a batch analysis'
                                        'OR the forward and reverse primers for a single analysis')
    # parser.add_argument('-v', '--version',
    #                     version='%(prog)s commit {}'.format(commit))
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-p', '--primerfile',
                        help='Specify path and name of text file with degenerate primers. The format of the file '
                             'should be: "primername1,forwardprimer1,reverseprimer1\n'
                             'primername2,forwardprimer2,reverseprimer2\n"')
    parser.add_argument('-f', '--forwardprimer',
                        help='Sequence of the forward primer')
    parser.add_argument('-r', '--reverseprimer',
                        help='Sequence of the reverse primer')

    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    DegeneratePrimers(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'
