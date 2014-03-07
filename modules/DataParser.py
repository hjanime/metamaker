#!/usr/bin/python2.7
"""
DataParser is a part of MetaMaker, used to extrapolate approximate sequencing 
technology statistics from real datasets.
"""

import sys
import numpy
import logging
import threading
import curses.ascii
from Bio import SeqIO

class MmpIO(object):
    """
    I/O handler for the super simple "MetaMaker Profile" file type.
    """
    
    @staticmethod
    def read(filename):
        """
        Reads a MetaMaker profile and returns a dictionary of values.
        """
        output = {}
        try:
            with open(filename, 'r') as f:
                for row in f:
                    if not row or row[0] == '#':
                        continue
                    key, value = row.strip().split(':')
                    key = key.strip().lower()
                    if key == 'key':
                        output['key'] = value.strip()
                    elif key == 'default reads':
                        output['reads'] = float(value)
                    elif key == 'read length mean':
                        if not 'length' in output:
                            output['length'] = {}
                        output['length']['mean'] = float(value)
                    elif key == 'read length var':
                        if not 'length' in output:
                            output['length'] = {}
                        output['length']['var'] = float(value)
                    elif key == 'quality mean':
                        if not 'quality' in output:
                            output['quality'] = {}
                        value = map(float, value.strip()[1:-1].split(', '))
                        output['quality']['mean'] = value
                    elif key == 'quality var':
                        if not 'quality' in output:
                            output['quality'] = {}
                        value = map(float, value.strip()[1:-1].split(', '))
                        output['quality']['var'] = value
        except Exception as e:
            print e
            pass
        
        return output
    
    @staticmethod
    def write(filename, data):
        """
        Writes a MetaMaker profile, given a formatted dictionary of values.
        """
        try:
            with open(filename, 'w') as f:
                f.write('Key: %s\n'              % data['key'])
                f.write('Default Reads: %s\n'    % data['reads'])
                f.write('Read Length Mean: %s\n' % float(data['length']['mean']))
                f.write('Read Length Var: %s\n'  % float(data['length']['var']))
                f.write('Quality Mean: %s\n'     % list(data['quality']['mean']))
                f.write('Quality Var: %s\n'      % list(data['quality']['var']))
        except Exception as e:
            print e
            pass
            

class DataParser( threading.Thread ):
    """
    DataParser parses fastq sequence files and uses least squares approximation 
    to fit approximation functions to the read length and quality statistics.
    """
    
    def __init__(self, infiles):
        """
        Init function, stores infiles, and initializes thread overhead and 
        globals.
        """
        threading.Thread.__init__(self)
        self.infiles  = infiles
        self.quality  = {}
        self.length   = {'mean':0, 'variance':0}
        self.variance = numpy.array([])
        self.count    = numpy.array([])
        self.settings = {'log':'DataParser',
                         'progress':False,
                         'output':None,
                         'profile dir':'../profiles',
                         'plot':False,
                        }
    
    def _parse_fastq(self, infile):
        """
        Parses a single fastq file and adds the results to the global 
        variables.
        """
        if self.settings['progress']:
            records = 0
            for record in SeqIO.parse(infile, 'fastq'):
                records += 1
            self.log.info("%s has %i records" % (filename, records))
        
        
        for i, record in enumerate(SeqIO.parse(infile, 'fastq')):
            
            qual = numpy.array(record.letter_annotations['phred_quality'], float)
            
            self.length['mean']     += len(record.seq)
            self.length['variance'] += len(record.seq)**2
            if len(record.seq) < self.length.get('min', 1e12):
                self.length['min'] = len(record.seq)
            if len(record.seq) > self.length.get('max', 0):
                self.length['max'] = len(record.seq)
            
            if 'mean' not in self.quality:
                self.quality['mean'] = numpy.array(qual)
                self.quality['max']  = qual
                self.quality['min']  = qual
                self.variance = qual**2
            else:
                for p, q in enumerate(qual):
                    if p >= len(self.quality['mean']):
                        self.quality['mean'] = numpy.append(
                                                    self.quality['mean'], q)
                        self.quality['max']  = numpy.append(
                                                    self.quality['max'],  q)
                        self.quality['min']  = numpy.append(
                                                    self.quality['min'],  q)
                        self.variance = numpy.append( self.variance, q**2 )
                    else:
                        self.quality['mean'][p] += q
                        self.variance[p] += q**2
                        if q > self.quality['max'][p]:
                            self.quality['max'][p] = q
                        if q < self.quality['min'][p]:
                            self.quality['min'][p] = q
            
            # counter to keep track of how many values are stored for each 
            # nucleotide position
            for p, q in enumerate(qual):
                if p >= len(self.count):
                    self.count = numpy.append(self.count, 1.0)
                else:
                    self.count[p] += 1.0
            
            if i and self.settings.get('progress'):
                if not self.settings['progress'] % i:
                    self.log.info("  parsed %i/%i records" % (i, records))
            if i > 999:
                break
    
    def run(self):
        """
        Starts the log and the starts parsing the input data.
        """
        self.log = logging.getLogger( self.settings['log'] )
        
        for infile in self.infiles:
            filename = infile.split('/')[-1]
            self.log.info("Parsing input file '%s'" % filename)
            
            self._parse_fastq(infile)
        
        self.statistics()
    
    def set(self, key, value):
        """
        Changes the settings value of 'key' to 'value'.
        """
        if key in self.settings:
            self.settings[key] = value
        else:
            raise Exception("Settings has no key '%s'." % key)
    
    def statistics(self):
        """
        Calculates and displays statistics from the global variables.
        """
        
        # convenience variables
        tot = float(self.count[0]) # total number of reads
        
        mean_reads    = round(tot / len(self.infiles))
        length_mean   = self.length['mean'] / tot
        length_var    = self.length['variance'] / tot - length_mean**2
        quality_means = self.quality['mean'] / self.count
        quality_vars  = self.variance / self.count - quality_means**2
        
        if self.settings['output']:
            output = self.settings['output']
        else: 
            output = self.infiles[0].split('/')[-1].split('.')[0]
        
        # least squares approximation coefficients
        
        m = numpy.polyfit(range(1,len(quality_means)+1), quality_means, 4)
        v = numpy.polyfit(range(1,len(quality_vars )+1), quality_vars,  2)
        
        # write results to log and file
        
        self.log.info("Key: %s" % output)
        self.log.info("Default reads: %i" % mean_reads)
        self.log.info("Read Length Mean: %i" % length_mean)
        self.log.info("Read Length Var: %3.2f" % length_var)
        self.log.info("Quality Mean: %s" % m)
        self.log.info("Quality Var:  %s" % v)
        
        outfile = '%s/%s.mmp' % (self.settings['profile dir'], output)
        MmpIO.write(outfile, {'key':output,
                              'reads':mean_reads,
                              'length':{'mean':length_mean,
                                        'var':length_var,},
                              'quality':{'mean':m,
                                         'var':v,},
                              })
        
        # plot stuff.
        if self.settings['plot']:
            try: 
                from matplotlib import pyplot
                fm = numpy.poly1d(m)
                fv = numpy.poly1d(v)
                x  = numpy.linspace(0, len(self.count)+1, len(self.count))
                
                pyplot.plot(x, quality_means, 'b.',
                            x, self.quality['max'], 'g.',
                            x, self.quality['min'], 'r.',
                            x, fm(x),         'r-',
                            x, fv(x),         'g--',
                            )
                pyplot.savefig('%s.png' % outfile[:-4], bbox_inches='tight')
                pyplot.show()
            except Exception as e:
                print e
                self.log.error(e)

if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser( description = __doc__ )
    parser.add_argument("input", nargs='+', help="input file(s) in fastq format.")
    parser.add_argument("-v", "--verbose", action='count', default=0, 
                        help="input file(s) in fastq format.")
    parser.add_argument("-d", "--profile_dir", default="../profiles",
                        help="Directory to store output MetaMaker profile.")
    parser.add_argument("-p", "--plot", default=False, action='store_true',
                        help=("Display a plot of the input data and functions. "
                              "(Requires matplotlib)") ) 
    parser.add_argument("-o", "--output", default=None, help="output filename")

    args = parser.parse_args()
    
    # Start Logging
    
    level = 50-args.verbose*10
    level = 10 if level < 10 else level
    
    log = logging.getLogger( "DataParser" )
    log.setLevel( level )
    
    formatter = logging.Formatter( ('%(asctime)s | %(name)s '
                                    '%(levelname)s: %(message)s') )
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel( level )
    console_handler.setFormatter(formatter)
    
    log.addHandler(console_handler)
    
    # Start the main script
    
    app = DataParser( args.input )
    app.set('output',      args.output)
    app.set('profile dir', args.profile_dir)
    app.set('plot',        args.plot)
    app.run()
    
    
    