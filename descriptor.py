#!/usr/bin/env python
#
"""
Code to generate basic DECam observing scripts 

""" 

from __future__ import (division, print_function)
import astropy.units as u

class Target(object):
    """An object that holds basic information on the observing target
    
    :param name:
        Target name
    
    :param ra:
        Target Right Ascension, astropy Angle
    
    :param dec:
        Target Declination, astropy Angle
    
    """
    
    def __init__(self, name, ra, dec):
        self.name = name
        self.ra = ra
        self.dec = dec
    
class Sequence(object):
    """ A sequence of target exposures to be taken 
    
    :param name:
        Sequence name
    
    :param target:
        Target object, containing information on the name and location of the target.
    
    :param exptype:
        Exposure type ['object' | 'dark' | 'dome flat' | 'sky flat' | 'zero']
    
    :param exptime:
        Exposure time, astropy Quantity
    
    :param band:
        Loaded filter cartridge ['u' | 'g' | 'r' | 'i' | 'z' | 'Y' | 'VR' | 'block' | 'pinhole' | 'None'] 
    
    :param seqtype:
        Sequence type ['single' | 'dither']
    
    :param count: (optional)
        Number of repeat exposures (default: ``1``)
    
    :param ditherpattern: (optional)
        Include the following dither pattern in the observing sequence ['line' | 'rectangle' | 'center+rectangle'] (default: ``'line'``)
    
    :param offset_ra: (optional)
        RA offset in dithering, astropy Quantity (default: ``90*u.arcsec``)
    
    :param offset_dec: (optional)
        Dec offset in dithering, astropy Quantity (default: ``90*u.arcsec``)
    
    :param num: (optional)
        Number of dithered exposures (default: ``3``)
    
    """
    
    def __init__(self, name, target, exptype, exptime, band, seqtype, count=1, ditherpattern='line', offset_ra=90*u.arcsec, offset_dec=90*u.arcsec, num=1):
        self.name = name
        self.target = target
        self.exptype = exptype
        self.exptime = exptime
        self.band = band
        self.seqtype = seqtype
        self.count = count
        
        # dithering parameters
        self.ditherpattern = ditherpattern
        self.offset_ra = offset_ra
        self.offset_dec = offset_dec
        self.num = num
        
        # make sure number of dithers is correct for a given ditherpattern
        self.num = self._numdither()
        
        # create list of exposures
        self.exposures = []
        self._make_exposures()
        
    def _get_position(self, index):
        """Return position (RA, dec) in degrees for a position in the observing sequence"""
        
        # define prefixes for calculating position offsets
        prefix_ra = {'line': range(self.num), 'rectangle': [0, 1, 1, 0], 'center+rectangle': [0, 1, -1, -1, 1]}
        prefix_dec = {'line': range(self.num), 'rectangle': [0, 0, 1, 1], 'center+rectangle': [0, 1, 1, -1, -1]}
        
        # find current prefix
        pra = prefix_ra[self.ditherpattern][index]
        pdec = prefix_dec[self.ditherpattern][index]
        
        # offset positions
        ra = self.target.ra + pra * self.offset_ra
        dec = self.target.dec + pdec * self.offset_dec
        
        return(ra, dec)
        
    def _make_exposures(self):
        """Populate list of exposures with the relevant details"""
        
        for i in range(self.num):
            
            ra, dec = self._get_position(i)
            
            exp = { 'seqid': self.name, 'seqnum': i+1, 'seqtot': self.num, 'expType': self.exptype, 'object': self.target.name, 'expTime': self.exptime, 'filter': self.band, 'RA': ra.deg, 'dec': dec.deg }
            
            if self.count>1:
                exp['count'] = self.count
            
            self.exposures.append(exp)
        
    def _numdither(self):
        """Return number of exposures for a given ditherpattern (different from 'line')"""
        
        dithermap = {'line': self.num, 'rectangle': 4, 'center+rectangle': 5}

        # check that valid ditherpattern is passed
        if self.ditherpattern not in dithermap.keys():
            raise ValueError('Dither pattern should be one of the following: {0}'.format(dithermap.keys()))
        else:
            return dithermap[self.ditherpattern]
        
    
class Script(object):
    """ Script holding multiple observing sequences, can print out to json file
    
    :param name:
        Filename for resulting script
    
    :param sequences: (optional)
        List of sequences in this script (default: ``[]``)
        
    :param dr: (optional)
        Output directory (default: ``"./"``)
    
    """
    
    def __init__(self, name, sequences=[], dr='./'):
        self.name = dr + name + '.json'
        self.sequences = sequences
        
    def add_sequence(self, sequence):
        """Add an exposure sequence to self.sequences"""
        
        self.sequences.extend(sequence)
    
    def add_sequences(self, sequences):
        """Add a list of exposure sequences to self.sequences"""
        
        for s in sequences:
            self.add_sequence(s)
        
    def write_json(self):
        """Write json script of provided observing sequences for DECam"""
        
        jf = open(self.name, 'w')
        jf.write('[' + '\n')
        
        for s in self.sequences:
            self.write_sequence(jf, s)
        
        # remove last comma
        jf.seek(-2,2)
        
        # close
        jf.write('\n' + ']')
        jf.close()
    
    def write_sequence(self, jf, s):
        """Write an exposure sequence to a json file"""
        
        jf.write(' {' + '\n')
        
        jf.write('  "seqid": "{0}", \n'.format(s['seqid']) )
        jf.write('  "seqnum": {0:d}, \n'.format(s['seqnum']) )
        jf.write('  "seqtot": {0:d}, \n'.format(s['seqtot']) )
        jf.write('  "expType": "{0}", \n'.format(s['expType']) )
        jf.write('  "object": "{0}", \n'.format(s['object']) )
        jf.write('  "expTime": {0:d}, \n'.format(s['expTime']) )
        jf.write('  "filter": "{0}", \n'.format(s['filter']) )
        jf.write('  "RA": {0:.6f}, \n'.format(s['RA']) )
        jf.write('  "dec": {0:.6f} \n'.format(s['dec']) )
        
        jf.write(' },' + '\n')
