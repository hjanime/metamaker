#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Welcome to MetaMaker!

This program uses NCBI data to create simulated datasets from known genomes.
The data is intended for use in testing metagenomic classifiers (or similar), 
or for benchmarking bioinformatic systems. It is in no way intended to be used 
in live research, other than as a control.

To use the system, simply input the variables above and click "generate". The 
options are:

Number of Genomes: The number of (randomly selected) genomes to include in the 
                   final dataset.

Taxa: The taxonomic classification of the randomly selected genomes. Currently
      only viruses are supported.

Species distribution: The amounts of the different species to be included in 
                      the dataset. "Uniform" means that equal amounts of all 
                      species are included, while "Exponential" means that the
                      amounts will be taken from an exponential function.

Technology: The system will attempt to approximate the output from the chosen 
            sequencing platform. Note that the system uses VERY simplified 
            generation functions, and may in some cases differ significantly 
            from real data.

Output filename: The output filename to use. The filename will be appended with
                 ".fastq" and ".key" for the individual output files.

The system will produce two output files, a fastq file of simulated reads from 
a pool of randomly selected genomes downloaded from ENTREZ. The system uses 
Sanger quality values (Phred+33) according to
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/?tool=pubmed. The second 
file is a "key file", a csv-file with summary information about the dataset.

================================================================================

"""

import time
import logging
from Tkinter import *
from modules.MetaMaker import MetaMaker

class MetaMakerGUI(object):
    """
    TkInter wrapper class for using MetaMaker with a GUI.
    """
    
    __name__ = "MetaMakerGUI"
    __version__ = '0.1.0'
    
    class GUILog(logging.Handler):
        """
        Log handler to print log messages to TkInter GUI console
        """
        def __init__(self, console):
            """
            Init function, just stores the console and calls super.
            """
            logging.Handler.__init__(self)
            self.console = console
        
        def emit(self, message):
            """
            Appends log message to the console.
            """
            formattedMessage = "%s\n" % self.format(message)
            
            self.console.configure(state=NORMAL)
            self.console.insert(END, formattedMessage)
            self.console.configure(state=DISABLED)
            self.console.see(END)
        
    class RunTimeFormatter(logging.Formatter):
        """
        Logging formatter that adds 'relativeTime' to the logging record. This 
        is a string of hours, minutes, seconds, and milliseconds formatted as
        HH:MM:SS.mm to be easily read by humans.
        """
        
        def h_m_s(self, msecs):
            """
            Formats a millisecond time as HH:MM:SS.mm.
            """
            ms = msecs%1e3
            s = msecs/1e3%60
            m = msecs/6e4%60
            h = msecs/36e4%24
            return "%02i:%02i:%02i.%02.f" % (h, m, s, ms/10)
        
        def format(self, record):
            """
            Interrupts super to add relativeTime to the record.
            """
            record.relativeTime = self.h_m_s(record.relativeCreated)
            return super(type(self), self).format(record)
    
    def __init__(self):
        """
        Sets up the tkinter root, as well as some basic parameters
        """
        self.gui = Tk()
        self.gui.wm_title("%s v. %s" % (self.__name__, self.__version__))
        self.gui.wm_minsize(800, 600)
        self._create_layout()
        self._start_log()
    
    def _create_layout(self):
        """
        Creates the GUI layout.
        """
        
        # Options
        
        Label(self.gui, text="No Genomes:").grid(row=0, column=0)
        self.no_genomes = Entry(self.gui)
        self.no_genomes.insert(0, 10)
        self.no_genomes.grid(row=0, column=1)
        
        Label(self.gui, text="Taxa:").grid(row=1,column=0)
        self.taxa = StringVar(self.gui)
        self.taxa.set("Viruses")
        drop = OptionMenu(self.gui, self.taxa, "Viruses")
        drop.grid(row = 1, column = 1)
        
        Label(self.gui, text="Species Distribution:").grid(row=2,column=0)
        self.distribution = StringVar(self.gui)
        self.distribution.set("Uniform")
        drop = OptionMenu(self.gui, self.distribution, "Uniform", "Exponential")
        drop.grid(row = 2, column = 1)
        
        Label(self.gui, text="Technology:").grid(row=3,column=0)
        self.technology = StringVar(self.gui)
        templates = MetaMaker.get_templates('keys', 'profiles')
        self.technology.set(templates[0])
        drop = OptionMenu(self.gui, self.technology, *templates)
        drop.grid(row = 3, column = 1)
        
        Label(self.gui, text="Output filename:").grid(row=4, column=0)
        self.filename = Entry(self.gui)
        self.filename.insert(0, "MetaMaker_output")
        self.filename.grid(row=4, column=1)
        
        # Start button
        self.start = Button(self.gui, text="Generate", command=self.generate)
        self.start.grid(row=5, column=1)
        
        # Logging console
        self.console = Text(self.gui, height=26, width=113, name="console")
        self.console.grid(row=10, column=0, columnspan=3)
        self.console.insert(END, __doc__)
    
    def _start_log(self):
        """
        Sets up the log parameters and registers a log with the program name.
        """
        
        self.log = logging.getLogger( self.__name__ )
        self.log.setLevel( logging.INFO )
    
        formatter = self.RunTimeFormatter( ('%(relativeTime)s | %(name)s '
                                            '%(levelname)s: %(message)s') )
    
        console_handler = self.GUILog(self.console)
        console_handler.setLevel( logging.INFO )
        console_handler.setFormatter(formatter)
    
        self.log.addHandler(console_handler)
    
    def generate(self):
        """
        Check the input and start generating a dataset!
        """
        no_of_genomes = 0
        read_distribution = self.distribution.get()
        try:
            no_of_genomes = int(self.no_genomes.get())
        except:
            self.log.error("Number of genomes must be an integer!")
            return
        
        outfile = self.filename.get()
        if not outfile:
            self.log.error("No filename for output files!")
            return
            
        self.log.info("Running %s genomes" % no_of_genomes)
        
        app = MetaMaker( outfile, no_of_genomes )
        app.set('log', self.__name__)
        app.set('distribution', read_distribution)
        app.set('template', self.technology.get())
        app.set('taxa', self.taxa.get())
        app.set('template_dir', 'profiles')
        app.start()
    
    def run(self):
        """
        Starts the gui
        """
        self.gui.mainloop()
    

if __name__ == '__main__':
    
    main = MetaMakerGUI()
    main.run()