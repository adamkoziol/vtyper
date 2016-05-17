#!/usr/bin/env python
from SPAdesPipeline.OLCspades.accessoryFunctions import *

__author__ = 'adamkoziol'


class PrimerObject(object):
    def single(self):
        """
        Creates an object with a sample name as well as the input primer sequences
        """
        # Initialise the metadata object and set the name
        metadata = MetadataObject()
        metadata.name = 'singleprimerset'
        # Initialise the forward and reverse attributes
        metadata.forward = self.forward
        metadata.reverse = self.reverse
        self.samples.append(metadata)

    def batch(self):
        """
        Creates objects with sample name, forward and reverse primer input sequences for each row in the primer file
        """
        from csv import DictReader
        # Load the primers from the file into a dictionary
        primerdict = DictReader(open(self.primerfile), fieldnames=self.primerheadings)
        for row in primerdict:
            # Initialise the metadata object and set the name
            metadata = MetadataObject()
            metadata.name = row['name']
            # Initialise the forward and reverse attributes
            metadata.forward = row['forward']
            metadata.reverse = row['reverse']
            self.samples.append(metadata)

    def __init__(self, inputobject):
        self.starttime = inputobject.starttime
        self.forward = inputobject.forward
        self.reverse = inputobject.reverse
        self.primerfile = inputobject.primerfile
        self.primerheadings = ['name', 'forward', 'reverse']
        self.samples = []
