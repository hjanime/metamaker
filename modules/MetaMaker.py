#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
MetaMaker is a basic tool for simulating metagenomic datasets. It downloads a 
number of genomes, splits them into reads and creates a fastq output file 
simulating a sequencing run.
"""


import os
import sys
import time
import numpy
import random
import logging
import threading
import curses.ascii
from Bio import Entrez, SeqIO
from DataParser import MmpIO

# Please add your own e-mail address here!
# It makes the people at Entrez super happy!
Entrez.email = "MetaMaker@slu.se"

class MetaMaker( threading.Thread ):
    """
    Viral Metagenomic dataset simulator.
    """
    
    def __init__(self, outfile = "output", num_genomes = 10):
        """
        Reads arguments and sets default settings.
        """
        threading.Thread.__init__(self)
        self.settings = {'num_genomes':num_genomes,
                         'outfile':"%s.fastq" % outfile,
                         'keyfile':None,
                         'taxa':'viruses',
                         'reads':1000,
                         'basepairs':200000,
                         'read length':200,
                         'length var':0,
                         'quality mean':[25],
                         'quality var':[10],
                         'distribution':'uniform',
                         'progress':False,
                         'log':'MetaMaker',
                         'template':None,
                         'template_dir':'../profiles',
                         }
        self.quality_cache = []
        self.variance_cache = []
        self._progress = 0
    
    def _get_tax_from_id(self, project_id):
        """
        Attempts to get taxonomic id from NCBI, given a sequencing project id.
        """
        project_term = "%s[BioProject]" % project_id
        tax_id = None
        self.log.debug('Asking NCBI for "%s" to get tax id.' % project_term)
        for retries in xrange(5):
            try:
                nuc_data = Entrez.read(Entrez.esearch("nucleotide", project_term))
                nuc_id = nuc_data['IdList'][0]
                summary = Entrez.read(Entrez.esummary(db="nucleotide", id=nuc_id))
                tax_id = summary[0]['TaxId']
                break
            except Exception as e:
                pass
        if tax_id:
            self.log.debug('Tax id found (%s)' % tax_id)
        else:
            self.log.debug('Tax id not found')
        
        return tax_id
    
    def _list_ncbi(self, max = 10000):
        """
        Lists (searches) NCBI entries for the specified taxa.
        """
        self.log.info('Getting list of %s from NCBI' % self.settings['taxa'])
        term = "%s[Organism Kingdom]" % self.settings['taxa']
        handle = Entrez.esearch("genome", term = term, retmax = max)
        results = Entrez.read(handle)
        self.log.info(' + Found %i %s' % (len(results['IdList']), 
                                          self.settings['taxa']))
        return results['IdList']
    
    def _list(self):
        """
        Wrapper function in case more sources are added.
        """
        id_list  = []
        id_list += self._list_ncbi()
        return id_list
    
    def _load_template(self, template):
        """
        Loads the run values from a MetaMaker profile template into the system.
        """
        template_dir = self.settings['template_dir']
        templates    = self.get_templates(None, template_dir)
        
        if template in templates:
            template_data = templates[template]
            
            self.settings['reads']       = template_data['reads']
            self.settings['read length'] = template_data['length']['mean']
            self.settings['length var']  = template_data['length']['var']
            self.settings['quality mean']= template_data['quality']['mean']
            self.settings['quality var'] = template_data['quality']['var']
            
            self.log.info("Using template '%s'" % self.settings['template'])
            self.log.info(" + Number of reads: %.1e" % self.settings['reads'])
            self.log.info(" + Read length    : %i±%i Bp" % (
                                                self.settings['read length'],
                                     numpy.sqrt(self.settings['length var']),))
        else:
            self.log.warning("Unknown template '%s', ignoring." % 
                                               self.settings['template'])
    
    def _make_dataset(self):
        """
        Creates the metadata for the project.
        """
        dataset = []
        avg_reads = self.settings['reads']/self.settings['num_genomes']
        
        if self.settings['distribution'].lower() == 'exponential':
            n = self.settings['reads']**(1.0/(self.settings['num_genomes']-1))
        
        ids = self._list()
        last = 0
        i = 0
        self.log.info("Making dataset")
        while i < self.settings['num_genomes']:
            genome_id = random.choice(ids)
            
            new = True
            for prev in dataset:
                if prev['genome_id'] == genome_id:
                    new = False
                    break
            if not new:
                continue
            
            summary = Entrez.read(Entrez.esummary(db="genome", id=genome_id))[0]
            
            tax_id = self._get_tax_from_id(summary['ProjectID'])
            if not tax_id:
                tax_id = '-'
                if self.settings['keyfile']:
                    continue
            
            data = {'genome_id':genome_id, 'def':summary['DefLine'], 
                    'project':summary['ProjectID'], 
                    'tax_id':tax_id}
            
            if self.settings['distribution'].lower() == 'uniform':
                data['reads'] = avg_reads
            elif self.settings['distribution'].lower() == 'exponential':
                data['reads'] = max(1, int(round(n**i - last)))
                last += data['reads']
            else:
                self.log.warning("WARNING: couldn't understand distribution" \
                                 " '%s', Defaulting to: Uniform" % \
                                 self.settings['distribution'])
                data['reads'] = avg_reads
            dataset += [data]
            i += 1
        return dataset
    
    def _make_read(self, seq):
        """
        Extracts a single read from a sequence and returns the read sequence as
        well as the position metadata.
        """
        
        # I seem to remember some early SOLiD reads being as short as 27bp,
        # so I'm putting that as a hard limit for short reads. Feel free to
        # change it if you feel like it though!
        
        read_length = int(self.settings['read length'])
        stdev       = numpy.sqrt(self.settings['length var'])
        min_length  = max(int(read_length-stdev), 27)
        
        start = random.randint(0, max(0, read_length-min_length))
        pos = (start, min(read_length, start + read_length))
        self.log.debug('Extracting read between %i-%i (%ibp).' % \
                       (pos[0], pos[1], pos[1]-pos[0]))
        return seq[pos[0]:pos[1]], pos
    
    def _make_quality(self, seq):
        """
        Simulates read quality from an error function.
        Qualities are in Sanger Fastq format (Phred+33), i.e. quality is 
        represented by an integer from 0 to 93, represented by the ascii 
        characters 33-126. 
        Errors are represented as 10^-0.0 (random base) to 10^-9.3 (super 
        accurate).
        
        ref: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/?tool=pubmed
        
        This might be re-written in the future using Biopythons QualityIO,
        http://www.biopython.org/DIST/docs/api/Bio.SeqIO.QualityIO-module.html
        """
        
        output = ""
        for i, q in enumerate(seq):
            if len(self.quality_cache) <= i:
                f = numpy.poly1d(self.settings['quality mean'])
                self.quality_cache += [f(len(self.quality_cache))]
            if len(self.variance_cache) <= i:
                v = numpy.poly1d(self.settings['quality var'])
                self.variance_cache += [v(len(self.variance_cache))]
            
            quality = self.quality_cache[i]
            quality += numpy.random.normal(0, numpy.sqrt(self.variance_cache[i]))
            quality = min(93, max(int(quality), 0))
            output += "%c" % (33+quality)
            
        return output
    
    def _write_csv(self, dataset, separator = ','):
        """
        Writes a csv file 
        """
        self.log.info('Creating Key file')
        header = ['Genome ID', 'Tax ID', 'Definition', 'Project', 'No. Reads']
        with open(self.settings['keyfile'], 'w') as key:
            key.write( "%s\n" % (separator.join(header)) )
            for i in dataset:
                data = [i['genome_id'], i['tax_id'], i['def'], 
                        i['project'],   i['reads']]
                key.write( "%s\n" % (separator.join(map(str,data))) )
        
    def progress(self):
        """
        Returns the progress of the current action.
        """
        return self._progress
    
    @staticmethod
    def get_templates(return_format = None, template_dir = '../profiles'):
        """
        Returns a list of allowed sequencing templates.
        This function uses MmpIO from DataParser to read templates from the
        template directory. It will generate an error is the template 
        directory isn't found.
        """
        templates = {}
        
        for template_file in os.listdir(template_dir):
            if template_file.split('.')[-1].lower() == 'mmp':
                data = MmpIO.read("%s/%s" % (template_dir, template_file))
                templates[data['key']] = data
        
        if return_format == "human":
            keys = templates.keys()
            if len(templates) < 2:
                return keys[0]
            return ", ".join(keys[:-1]) + " or " + keys[-1]
        if return_format == "keys":
            return templates.keys()
        return templates
    
    def run(self):
        """
        Starts the job of creating a metagenomic sample set.
        """
        self.log = logging.getLogger( self.settings['log'] )
        
        if self.settings['template']:
            self._load_template( self.settings['template'] )
        
        dataset = self._make_dataset()
        
        # Print debug information about the dataset
        self.log.debug('DATASET:')
        tot_reads = 0
        for i in dataset:
            self.log.debug("%i\t%s" % (i['reads'], i['def']))
            tot_reads += i['reads']
        self.log.debug("TOTAL READS: %i" % tot_reads)
        
        # Create the key file
        if self.settings['keyfile']:
            self._write_csv(dataset)
        
        # Start creating the fastq output file.
        out = open(self.settings['outfile'], 'w')
        for metadata in dataset:
            self._progress = 0.0
            
            self.log.info("* Parsing %s" % metadata['def'])
            self.log.info("  * Downloading")
            
            project_term = "%s[BioProject]" % metadata['project']
            for tries in xrange(5):
                try:
                    handle = Entrez.esearch("nucleotide", project_term)
                    nuc_id = Entrez.read(handle)['IdList'][0]
                    data = Entrez.efetch(db="nucleotide", id=nuc_id, 
                                         rettype="gb",    retmode="text")
                    break
                except:
                    self.log.info("    * Retrying")
                    pass
            
            self.log.info("  * Creating Reads" )
            
            for record in SeqIO.parse(data,"gb"):
                # TODO: make use of several records if present
                for i in xrange(int(metadata['reads'])):
                    seq = None
                    while not seq:
                        seq, meta = self._make_read(record.seq)
                        quality = self._make_quality(seq)
                    
                    # apply quality to read
                    seq = list(seq)
                    for j, q in enumerate(quality):
                        if numpy.random.random() < (10**-((ord(q)-33)/10.0)):
                            seq[j] = 'actg'[numpy.random.randint(4)]
                    seq = "".join(seq)
                    header = "@%s|ref:%s-%i|pos:%i-%i\n" % (record.id, 
                                                        metadata['genome_id'],
                                                        i, 
                                                        meta[0], meta[1])
                    out.write(header)
                    out.write("%s\n" % seq)
                    out.write("+\n%s\n" % quality)
                    self._progress = (i+1)/float(int(metadata['reads']))
                break
        
        out.close()
        self._progress = -1
        self.log.info("Finished. All went well!")
    
    def set(self, key, value):
        """
        Sets a value in the settings dictionary.
        """
        if key in self.settings:
            if key == 'keyfile' and not value.endswith('.csv'):
                self.settings[key] = "%s.csv" % value
            else:
                self.settings[key] = value
        else:
            raise Exception("Unknown key '%s' in settings." % key)
    

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser( description = __doc__,
                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("number_of_genomes", help="Number of genomes.", 
                        type=int)
    parser.add_argument("-o", "--output", help="Output fastq filename",
                        default="meta_output")
    parser.add_argument("-k", "--keyfile", help="key filename.",
                        default=None)
    parser.add_argument("-r", "--reads", help="Number of reads.",
                        default="50M")
    parser.add_argument("-b", "--basepairs", help="Total number of basepairs.",
                        default="10G")
    parser.add_argument("-l", "--readlength", help="Read length",
                        default="200")
    parser.add_argument("-d", "--distribution", 
                        help="Read distribution, 'uniform' or 'exponential'",
                        default="uniform")
    parser.add_argument("-e", "--error_function", nargs="+", type=float,
                        default = [25.0],
                        help="Factors for the error approximation equation.")
    parser.add_argument("-c", "--error_variance", nargs="+", type=float,
                        default = [0],
                        help=("Factors for the error variance approximation "
                              "equation.") )
    parser.add_argument("-p", "--progress", default=False,
                        help="Display progress information for long tasks.")
    parser.add_argument("-t", "--template", default=None,
                        help=("Sequencing template to use for read generation. "
                              "Overrides reads, basepairs, readlength and "
                              "error_function. Valid options are %s") % \
                              MetaMaker.get_templates('human'))
    parser.add_argument("-x", "--taxa", default="viruses",
                        help=("Taxonomic identifier of the species to "
                              "download."))
    parser.add_argument("-v", "--verbose", action = "count", default = 0, 
                        help="Set output Verbosity")
    
    args = parser.parse_args()
    
    for arg in ['reads', 'basepairs', 'readlength']:
        if eval("args.%s" % arg)[-1] in ['K', 'k']:
            exec("args.%s = int(args.%s[:-1])*1000" % (arg, arg))
        elif eval("args.%s" % arg)[-1] in ['M', 'm']:
            exec("args.%s = int(args.%s[:-1])*1000000" % (arg, arg))
        elif eval("args.%s" % arg)[-1] in ['G', 'g']:
            exec("args.%s = int(args.%s[:-1])*1000000000" % (arg, arg))
        else:
            exec("args.%s = int(args.%s)" % (arg, arg))

    level = 50-args.verbose*10
    level = 10 if level < 10 else level
    
    log = logging.getLogger( "MetaMaker" )
    log.setLevel( level )
    
    formatter = logging.Formatter( ('%(asctime)s | %(name)s '
                                    '%(levelname)s: %(message)s') )
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel( level )
    console_handler.setFormatter(formatter)
    
    log.addHandler(console_handler)
    
    app = MetaMaker( args.output, args.number_of_genomes )
    app.set('keyfile',      args.keyfile)
    app.set('taxa',         args.taxa)
    app.set('reads',        args.reads)
    app.set('basepairs',    args.basepairs)
    app.set('read length',  args.readlength)
    app.set('quality mean', args.error_function)
    app.set('distribution', args.distribution)
    app.set('template',     args.template)
    app.set('progress',     args.progress)
    
    app.run()