#!/usr/bin/env python3

# -*- coding: utf-8 -*-
#############################################################
## Application:               panGraphViewer
## Author:                    Yuxuan Yuan
## Email:                     yuxuan.yuan@outlook.com
## Last modification Date:    08/01/2022
## Copyright: Copyright (c) 2021-2022 the Ting-Fung Chan Lab
## CUHK, Hong Kong SAR, China
#############################################################

#======================== libraries ==========================
import os
import sys
import shlex
import codecs
import shutil
import datetime
import traceback
import threading
from configparser import ConfigParser, ConverterMapping
import webbrowser
import subprocess
from time import sleep
from natsort import natsorted
from networkx.drawing.layout import rescale_layout 
#from scripts.screenClipping import *
#from scripts.screenshot import *
from scripts.utilities import *
from scripts.vcf2rGFA import *
from scripts.gfa2rGFA import *
from scripts.NodesInfoBrowser import *

from scripts.panGraph import *

if sys.platform == 'darwin':
    from scripts.entireUIDesign_mac import Ui_panGraphViewer
    from scripts.UI_vcf2rGFA_mac import Ui_vcf2rGFA
    from scripts.About_mac import Ui_About
    from scripts.renameFastaHeaderUI_mac import *
    from scripts.GraphModificationUI_mac import *
    from scripts.NodeShapesUI_mac import *
else:
    from scripts.About import Ui_About
    from scripts.UI_vcf2rGFA import Ui_vcf2rGFA
    from scripts.entireUIDesign import Ui_panGraphViewer
    from scripts.renameFastaHeaderUI import *
    from scripts.GraphModificationUI import *
    from scripts.NodeShapesUI import *
try:
    from PyQt5.QtGui import*
    from PyQt5 import (
        sip,
        QtCore,
        QtGui,
        QtWidgets,
        QtNetwork,
        QtWebEngineCore,
        QtWebChannel,
        QtPrintSupport,
        QtWebEngineWidgets,
    )
    from PyQt5.QtWebEngineWidgets import QWebEngineView
    from PyQt5.QtWidgets import (
        QApplication,
        QMainWindow,
        QWidget,
        QInputDialog,
        QLineEdit,
        QFileDialog,
        QVBoxLayout
    )
    from PyQt5.QtCore import (
        QObject,
        QThread,
        pyqtSignal,
        pyqtSlot,
        QRunnable,
        QThreadPool,
    )

except ImportError:
    print("\nOops! Something is wrong when load the python libaries")
    print("Please follow the instruction to install the python dependencies and then try again!\n")

# Global variable
tool = os.path.realpath(__file__)
toolPath = os.path.dirname(tool)
manualPath = os.path.join(toolPath,'..', 'doc', 'Manual.html').replace('\\', '/')
configfile = os.path.join(toolPath, 'config_default.ini')
copied_config = os.path.join(toolPath, 'config.ini')
with open(copied_config, 'wb') as f:
    with open(configfile, 'rb') as c:
        for line in c:
            f.write(line)

#=============================================================
class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.
    Supported signals are:
    - finished: No data
    - error:`tuple` (exctype, value, traceback.format_exc() )
    - result: `object` data returned from processing, anything
    - progress: `tuple` indicating progress metadata
    '''
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(tuple)

class Worker(QRunnable):
    '''
    Worker thread
    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.
    '''
    def __init__(self, fn, *args):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn(*self.args)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done

"""
class options():
    def __init__(self, vcf, fasta, backbone, outdir, nthread, chr=None, prefix=None):
        self.vcf = vcf
        self.fasta = fasta
        self.backbone = backbone
        self.outdir = outdir
        self.chr = chr
        self.prefix = prefix
        self.nthread = nthread
"""

class NodesInfoBrowser(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.ui = Ui_NodesInfoCanvas()
        self.ui.setupUi(self)

class about(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.ui = Ui_About()
        self.ui.setupUi(self)
        self.ui.closeButton.clicked.connect(self.close)

class vcf_2_rGFA(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.ui = Ui_vcf2rGFA()
        self.ui.setupUi(self)
        self.fasta = None
        self.vcf = None
        self.backbone = None
        self.outDirPath = None
        self.cvt = None
        self.threads = 4
        self.chrs = None
        self.threadpool = QThreadPool()

        # =================== UI elements =====================
        self.ui.fastaSelect.clicked.connect(self.select_fasta)
        self.ui.fastaClear.clicked.connect(self.clear_fasta)
        self.ui.vcfSelect.clicked.connect(self.select_vcf)
        self.ui.vcfClear.clicked.connect(self.clear_vcf)
        self.ui.outSelect.clicked.connect(self.outDir_select)
        self.ui.outClear.clicked.connect(self.outDir_clear)
        self.ui.covert.clicked.connect(self.convertWorker)

    def select_fasta(self):  # load GFA file
        if sys.platform == 'win32':
            fasta = QFileDialog.getOpenFileName(self, 'Select the backbone genome assembly', '',
                                                'fasta(*fa *fasta *fna)')
            self.fasta = codecs.decode(str(fasta)[1:-1].split(',')[0][1:-1], 'unicode_escape')
        else:
            self.fasta = \
            QFileDialog.getOpenFileName(self, 'Select the backbone genome assembly', '', 'fasta(*fa *fasta *fna)')[0]
        self.fasta = os.path.abspath(self.fasta)
        
        if checkFile(self.fasta) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the input FASTA file. Please check !',
            QtWidgets.QMessageBox.Ok)
            self.ui.fastaPath.clear()
            self.fasta = None
        else:
            self.ui.fastaPath.setText(self.fasta)

    def clear_fasta(self):
        self.ui.fastaPath.clear()
        self.fasta = None

    def select_vcf(self):
        if sys.platform == 'win32':
            vcf = QFileDialog.getOpenFileName(self, 'Select the vcf file', '', 'vcf(*vcf)')
            self.vcf = codecs.decode(str(vcf)[1:-1].split(',')[0][1:-1], 'unicode_escape')
        else:
            self.vcf = QFileDialog.getOpenFileName(self, 'Select the vcf file', '', 'vcf(*vcf)')[0]
        self.vcf = os.path.abspath(self.vcf)
        if checkFile(self.vcf) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the input VCF file. Please check !',
            QtWidgets.QMessageBox.Ok) 
            self.ui.vcfPath.clear()
            self.vcf = None
        else:
            self.ui.vcfPath.setText(self.vcf)           

    def clear_vcf(self):
        self.ui.vcfPath.clear()
        self.vcf = None

    def outDir_select(self):
        self.outDirPath = QFileDialog.getExistingDirectory()
        self.outDirPath = os.path.abspath(self.outDirPath)
        if checkDir(self.outDirPath) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the outDir. Please check !',
            QtWidgets.QMessageBox.Ok)
            self.ui.outDirPath.clear()
            self.outDirPath = None
        else:
            self.ui.outDirPath.setText(self.outDirPath)

    def outDir_clear(self):
        self.ui.outDirPath.clear()
        self.outDirPath = None

    def checkParseElem(self):
        if self.fasta == None or len(self.fasta) == 0:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select a FASTA file',
                                           QtWidgets.QMessageBox.Ok)
            return -1
        if self.vcf == None or len(self.vcf) == 0:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select a VCF file',
                                           QtWidgets.QMessageBox.Ok)
            return -2
        if self.outDirPath == None or len(self.outDirPath) == 0:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please specify the output directory',
                                           QtWidgets.QMessageBox.Ok)
            return -3
        return 0

    def showProgress(self):
        pass

    def getBB(self):
        self.backbone = self.ui.bbName.text().rstrip(' ')

    def getThreads(self):
        self.threads = self.ui.threads.text().rstrip(' ')

    def getChrs(self):
        #self.chrs = self.ui.chrTextEdit.toPlainText().strip().split('\n')
        manualInput = self.ui.chrTextEdit.toPlainText().split('\n')
        self.chrs = ''
        for i in manualInput:
            if i.strip() != '':
                self.chrs = self.chrs + i.strip() + ' '

    def checkElements(self):
        if self.checkParseElem() == 0:
            self.backbone = None
            self.getBB()
            if self.backbone == '' or self.backbone == None:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Please input a backbone name',
                                               QtWidgets.QMessageBox.Ok)
                return -1
            self.getThreads()

            try:
                self.threads = int(self.threads)
            except ValueError:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'The input of threads must be an integer !',
                                               QtWidgets.QMessageBox.Ok)
                return -2

            self.getChrs()
            
            if 'All' in self.chrs.strip().split(' ') and len(self.chrs.strip().split(' ')) >1:
                QtWidgets.QMessageBox.question(self, 'Warning !', "'All' and other chromosome names cannot be input together. Please check !",
                                               QtWidgets.QMessageBox.Ok)
                return -3              
            if self.chrs.strip() == '':
                QtWidgets.QMessageBox.question(self, 'Warning !', "Please input a chromosome name(s) or 'All' to convert !",
                                                QtWidgets.QMessageBox.Ok)
                return -4
            
            #self.chrs = None if self.chrs == ['All'] else self.chrs
            self.chrs = None if self.chrs.strip() == 'All' else self.chrs
            self.Date_1 = datetime.datetime.today()
            self.ui.convertStatus.setStyleSheet('color: blue')
            self.ui.convertStatus.setText('Running ...')
            return 0

    def convert(self):
        path = os.path.join(toolPath, 'scripts')
        script = os.path.join(path, 'vcf2rGFA.py')
        if self.chrs != None:
            cmd = f'"{sys.executable}" "{script}" -f "{self.fasta}" -v "{self.vcf}" -n "{self.threads}" -c "{self.chrs}" -o "{self.outDirPath}" -b "{self.backbone}"'
        else:
            cmd = f'"{sys.executable}" "{script}" -f "{self.fasta}" -v "{self.vcf}" -n "{self.threads}" -o "{self.outDirPath}" -b "{self.backbone}"'
        
        if sys.platform == 'win32':
            self.retcode = subprocess.call(shlex.split(cmd), shell = True)      
        else:
            self.retcode = subprocess.call(shlex.split(cmd), shell = False)
        """
        args = options(self.vcf, self.fasta, self.backbone, self.outDirPath, self.threads, self.chrs)
        try:
            self.retcode = VCF2rGFA(args).run()  
        except SystemExit:
            self.retcode = VCF2rGFA(args)
        """

    def endConvert(self):
        self.Date_2 = datetime.datetime.today()
        Time_duration = (self.Date_2 - self.Date_1).total_seconds()
        self.ui.convertStatus.setStyleSheet('color: green')
        self.ui.convertStatus.setText('Finished in %.2fs !' % Time_duration)
        try:
            if self.retcode == 2:
                QtWidgets.QMessageBox.question(self, 'Error!', 'Oops ! Unmathced FASTA file is given. Please check the log info. Abort !',
                QtWidgets.QMessageBox.Ok) 

            if self.retcode == 4: 
                QtWidgets.QMessageBox.question(self, 'Warning!', f"Oops ! Illegal chromosome name(s) in '{self.chrs}'. Ignored !",
                QtWidgets.QMessageBox.Ok)  
    
            if self.retcode == 1: 
                QtWidgets.QMessageBox.question(self, 'Error!', 'Oops ! No sample columns found in the input VCF. Abort !',
                QtWidgets.QMessageBox.Ok)  
                                           
        except: 
            pass

    def convertWorker(self):
        if self.checkElements() == 0:
            self.worker = Worker(self.convert)
            self.threadpool.start(self.worker)
            self.worker.signals.finished.connect(self.endConvert)
            self.worker.signals.progress.connect(self.showProgress)

class renameFaHeader(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.ui = Ui_header()
        self.ui.setupUi(self)
        self.fasta = None
        self.sampleName = None
        self.delim = "||"
        self.outDirPath = None
        self.threadpool = QThreadPool()

        self.ui.fastaSelect.clicked.connect(self.select_fasta)
        self.ui.fastaClear.clicked.connect(self.clear_fasta)
        self.ui.outDirSelect.clicked.connect(self.outDir_select)
        self.ui.outDirClear.clicked.connect(self.outDir_clear)
        self.ui.startRename.clicked.connect(self.renameWorker)

    def select_fasta(self):
        if sys.platform == 'win32':
            fasta = QFileDialog.getOpenFileName(self, 'Select the backbone genome assembly', '',
                                                'fasta(*fa *fasta *fna)')
            self.fasta = codecs.decode(str(fasta)[1:-1].split(',')[0][1:-1], 'unicode_escape')
        else:
            self.fasta = \
            QFileDialog.getOpenFileName(self, 'Select the backbone genome assembly', '', 'fasta(*fa *fasta *fna)')[0]
        self.fasta = os.path.abspath(self.fasta)
        
        if checkFile(self.fasta) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the input FASTA file. Please check !',
            QtWidgets.QMessageBox.Ok)
            self.ui.self.ui.fastaPath.clear()
            self.fasta = None
        else:
            self.ui.fastaPath.setText(self.fasta)

    def clear_fasta(self):
        self.ui.fastaPath.clear()
        self.fasta = None

    def outDir_select(self):
        self.outDirPath = QFileDialog.getExistingDirectory()
        self.outDirPath = os.path.abspath(self.outDirPath)
        if checkDir(self.outDirPath) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the outDir. Please check !',
            QtWidgets.QMessageBox.Ok)
            self.ui.outDir.clear()
            self.outDirPath = None
        else:
            self.ui.outDir.setText(self.outDirPath)

    def outDir_clear(self):
        self.ui.outDir.clear()
        self.outDirPath = None

    def getSampleName(self):
        self.sampleName = self.ui.sampleName.text().rstrip(' ')

    def getDelim(self):
        self.delim = self.ui.delim.text().rstrip(' ')

    def startSignal(self):
        self.Date_1 = datetime.datetime.today()
        self.ui.renameStatus.setStyleSheet('color: blue')
        self.ui.renameStatus.setText('Running ...')

    def endSignal(self):
        self.Date_2 = datetime.datetime.today()
        Time_duration = (self.Date_2 - self.Date_1).total_seconds()
        self.ui.renameStatus.setStyleSheet('color: green')
        self.ui.renameStatus.setText('Finished in %.2fs !' % Time_duration)

    def checkElements(self):
        self.getSampleName()
        self.getDelim()

        if self.fasta == None:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select a fasta file', QtWidgets.QMessageBox.Ok)
            return -1
        if self.outDirPath == None:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please specify an output directory', QtWidgets.QMessageBox.Ok)
            return -2
        if self.sampleName == None or self.sampleName == '':
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please specify the sample name', QtWidgets.QMessageBox.Ok)
            return -3
        if self.delim == '' or self.delim != "||":
            QtWidgets.QMessageBox.question(self, 'Warning !', "Please use '||' as a delimiter", QtWidgets.QMessageBox.Ok)
            return -4
        return 0

    def showProgress(self):
        pass

    def rename(self):
        self.startSignal()
        path = os.path.join(toolPath, 'scripts')
        script = os.path.join(path, 'renameFastaHeader.py')
        if self.delim == "||":
            cmd = f'"{sys.executable}" "{script}" -f "{self.fasta}" -n "{self.sampleName}" -o "{self.outDirPath}"'
        else:
            cmd = f'"{sys.executable}" "{script}" -f "{self.fasta}" -n "{self.sampleName}" -o "{self.outDirPath}" -d "{self.delim}"'
        
        if sys.platform == 'win32':
            subprocess.call(shlex.split(cmd), shell = True)      
        else:
            subprocess.call(shlex.split(cmd), shell = False)

    def renameWorker(self):
        if self.checkElements() == 0:
            self.worker = Worker(self.rename)
            self.threadpool.start(self.worker)
            self.worker.signals.finished.connect(self.endSignal)
            self.worker.signals.progress.connect(self.showProgress)

class NodeShapes(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.ui = Ui_nodeShapes()
        self.ui.setupUi(self)

        self.VisOptions = ['dot', 'ellipse', 'square', 
        'database', 'box', 'text', 
        'diamond', 'star', 'triangle', 
        'triangleDown']

        self.cyOptions = ["ellipse", "triangle", "round-triangle", 
        "rectangle", "round-rectangle", "bottom-round-rectangle",
        "cut-rectangle", "barrel", "rhomboid", "diamond", 
        "round-diamond", "pentagon", "round-pentagon", "hexagon", 
        "round-hexagon", "concave-hexagon", "heptagon", "round-heptagon",
        "octagon", "round-octagon", "star", "tag", "round-tag", "vee"]

        self.ui.visBB.activated[str].connect(self.backbone)
        self.ui.visSNP.activated[str].connect(self.snp)
        self.ui.visDEL.activated[str].connect(self.deletion)
        self.ui.visINS.activated[str].connect(self.insertion)
        self.ui.visINV.activated[str].connect(self.inversion)
        self.ui.visDUP.activated[str].connect(self.duplication)
        self.ui.visTRANS.activated[str].connect(self.translocation)
        self.ui.cyBB.activated[str].connect(self.backbone)
        self.ui.cySNP.activated[str].connect(self.snp)
        self.ui.cyDEL.activated[str].connect(self.deletion)
        self.ui.cyINS.activated[str].connect(self.insertion)
        self.ui.cyINV.activated[str].connect(self.inversion)
        self.ui.cyDUP.activated[str].connect(self.duplication)
        self.ui.cyTRANS.activated[str].connect(self.translocation)
        self.ui.change.clicked.connect(self.saveChanges)
        self.ui.reset.clicked.connect(self.reset)

    def backbone(self):
        self.visBB = self.ui.visBB.currentText()
        self.cyBB = self.ui.cyBB.currentText()

    def snp(self):
        self.visSNP = self.ui.visSNP.currentText()
        self.cySNP = self.ui.cySNP.currentText()

    def deletion(self):
        self.visDEL = self.ui.visDEL.currentText()
        self.cyDEL = self.ui.cyDEL.currentText()

    def insertion(self):
        self.visINS = self.ui.visINS.currentText()
        self.cyINS = self.ui.cyINS.currentText()

    def inversion(self):
        self.visINV = self.ui.visINV.currentText()
        self.cyINV = self.ui.cyINV.currentText()

    def duplication(self):
        self.visDUP = self.ui.visDUP.currentText()
        self.cyDUP = self.ui.cyDUP.currentText()

    def translocation(self):
        self.visTRANS = self.ui.visTRANS.currentText()
        self.cyTRANS = self.ui.cyTRANS.currentText()

    def saveChanges(self):
        self.backbone()
        self.snp()
        self.deletion()
        self.insertion()
        self.inversion()
        self.duplication()
        self.translocation()

        config = ConfigParser()
        config.optionxform = str
        config.read(copied_config)
        config.set("nodes", "BB_shape", self.visBB)
        config.set("nodes", "SNP_shape", self.visSNP)
        config.set("nodes", "DEL_shape", self.visDEL)
        config.set("nodes", "DUP_shape", self.visDUP)
        config.set("nodes", "INS_shape", self.visINS)
        config.set("nodes", "INV_shape", self.visINV)
        config.set("nodes", "BND_shape", self.visTRANS)
        config.set("cytoscape", "BB_shape", self.cyBB)
        config.set("cytoscape", "SNP_shape", self.cySNP)
        config.set("cytoscape", "DEL_shape", self.cyDEL)
        config.set("cytoscape", "DUP_shape", self.cyDUP)
        config.set("cytoscape", "INS_shape", self.cyINS)
        config.set("cytoscape", "INV_shape", self.cyINV)
        config.set("cytoscape", "BND_shape", self.cyTRANS)
        with open(copied_config, 'w') as configfile:
            config.write(configfile)
        QtWidgets.QMessageBox.question(self, 'Information', 'The settings of node shapes have been saved !',
        QtWidgets.QMessageBox.Ok)

    def reset(self):
        visBB = self.VisOptions
        visBB.remove('dot')
        visBB.insert(0, 'dot')
        self.ui.visBB.clear()
        self.ui.visBB.addItems(visBB)        
        cyBB = self.cyOptions
        cyBB.remove('ellipse')
        cyBB.insert(0, 'ellipse')
        self.ui.cyBB.clear()
        self.ui.cyBB.addItems(cyBB)
        
        visSNP = self.VisOptions
        visSNP.remove('dot')
        visSNP.insert(0, 'dot')
        self.ui.visSNP.clear()
        self.ui.visSNP.addItems(visSNP)
        cySNP = self.cyOptions
        cySNP.remove('ellipse')
        cySNP.insert(0, 'ellipse')
        self.ui.cySNP.clear()
        self.ui.cySNP.addItems(cySNP)

        visDEL = self.VisOptions
        visDEL.remove('triangle')
        visDEL.insert(0, 'triangle')
        self.ui.visDEL.clear()
        self.ui.visDEL.addItems(visDEL)

        cyDEL = self.cyOptions
        cyDEL.remove('concave-hexagon')
        cyDEL.insert(0, 'concave-hexagon')
        self.ui.cyDEL.clear()
        self.ui.cyDEL.addItems(cyDEL)
        
        visINS = self.VisOptions
        visINS.remove('triangleDown')
        visINS.insert(0, 'triangleDown')
        self.ui.visINS.clear()
        self.ui.visINS.addItems(visINS)

        cyINS = self.cyOptions
        cyINS.remove('triangle')
        cyINS.insert(0, 'triangle')
        self.ui.cyINS.clear()
        self.ui.cyINS.addItems(cyINS)

        visINV = self.VisOptions
        visINV.remove('text')
        visINV.insert(0, 'text')
        self.ui.visDUP.clear()
        self.ui.visDUP.addItems(visINV)
        cyINV = self.cyOptions
        cyINV.remove('vee')
        cyINV.insert(0, 'vee')
        self.ui.cyDUP.clear()
        self.ui.cyDUP.addItems(visINV)

        visDUP = self.VisOptions
        visDUP.remove('database')
        visDUP.insert(0, 'database')
        self.ui.visDUP.clear()
        self.ui.visDUP.addItems(visDUP)
        cyDUP = self.cyOptions
        cyDUP.remove('diamond')
        cyDUP.insert(0, 'diamond')
        self.ui.cyDUP.clear()
        self.ui.cyDUP.addItems(cyDUP)

        visTRANS = self.VisOptions
        visTRANS.remove('star')
        visTRANS.insert(0, 'star')
        self.ui.visTRANS.clear()
        self.ui.visTRANS.addItems(visTRANS)
        cyTRANS = self.cyOptions
        cyTRANS.remove('star')
        cyTRANS.insert(0, 'star')
        self.ui.cyTRANS.clear()
        self.ui.cyTRANS.addItems(cyTRANS)

        self.backbone()
        self.snp()
        self.deletion()
        self.insertion()
        self.inversion()
        self.duplication()
        self.translocation()

        config = ConfigParser()
        config.optionxform = str
        config.read(copied_config)
        config.set("nodes", "BB_shape", self.visBB)
        config.set("nodes", "SNP_shape", self.visSNP)
        config.set("nodes", "DEL_shape", self.visDEL)
        config.set("nodes", "DUP_shape", self.visDUP)
        config.set("nodes", "INS_shape", self.visINS)
        config.set("nodes", "INV_shape", self.visINV)
        config.set("nodes", "BND_shape", self.visTRANS)
        config.set("cytoscape", "BB_shape", self.cyBB)
        config.set("cytoscape", "SNP_shape", self.cySNP)
        config.set("cytoscape", "DEL_shape", self.cyDEL)
        config.set("cytoscape", "DUP_shape", self.cyDUP)
        config.set("cytoscape", "INS_shape", self.cyINS)
        config.set("cytoscape", "INV_shape", self.cyINV)
        config.set("cytoscape", "BND_shape", self.cyTRANS)
        with open(copied_config, 'w') as configfile:
            config.write(configfile)
        QtWidgets.QMessageBox.question(self, 'Information', 'The settings of node shapes have been initialized !',
        QtWidgets.QMessageBox.Ok)

class GraphModify(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.ui = Ui_GraphModi()
        self.ui.setupUi(self)
        self.interaction = 'No'
        self.layoutModification = 'No'
        self.addRemoveNodes = 'No'
        self.height = '980'
        self.width = '1400'
        self.maxNodesDisplay = '200'
        self.maxNodeLenDisplay = '1000000'
        self.geneNodeOverlapCntThreshold = '2'

        self.ui.interaction.activated[str].connect(self.Interaction)
        self.ui.layoutModification.activated[str].connect(self.LayoutModification)
        self.ui.addRemoveNodes.activated[str].connect(self.AddRemoveNode)
        self.ui.change.clicked.connect(self.saveChanges)
        self.ui.reset.clicked.connect(self.reset)

    def Interaction(self):
        self.interaction = self.ui.interaction.currentText()

    def LayoutModification(self):
        self.layoutModification = self.ui.layoutModification.currentText()

    def AddRemoveNode(self):
        self.addRemoveNodes = self.ui.addRemoveNodes.currentText()

    def getHeight(self):
        try:
            self.height = int(self.ui.height.text().rstrip(' '))
            if self.height <= 0:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal height input. Please check!', QtWidgets.QMessageBox.Ok)
                return -1
            self.height = str(self.height)
        except ValueError:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal height input. Please check!', QtWidgets.QMessageBox.Ok)
            return -1

    def getWidth(self):
        try:
            self.width = int(self.ui.width.text().rstrip(' '))
            if self.width <= 0: 
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal width input. Please check!', QtWidgets.QMessageBox.Ok)
                return -2                
            self.width = str(self.width)
        except ValueError:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal width input. Please check!', QtWidgets.QMessageBox.Ok)
            return -2

    def getMaxNodesDisplay(self):
        try:
            self.maxNodesDisplay = int(self.ui.maxNodesDisplay.text().rstrip(' '))
            if self.maxNodesDisplay <= 0:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal value for max num. of nodes. Please check!', QtWidgets.QMessageBox.Ok)
                return -3               
            self.maxNodesDisplay = str(self.maxNodesDisplay)
        except ValueError:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal value for max num. of nodes. Please check!', QtWidgets.QMessageBox.Ok)
            return -3

    def getMaxNodeLenDisplay(self):
        try:
            self.maxNodeLenDisplay = int(self.ui.nodeLen.text().rstrip(' '))
            if self.maxNodeLenDisplay <= 0:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal value for max length of sequence to display. Please check!', QtWidgets.QMessageBox.Ok)
                return -4                
            self.maxNodeLenDisplay = str(self.maxNodeLenDisplay)
        except ValueError:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal value for max length of sequence to display. Please check!', QtWidgets.QMessageBox.Ok)
            return -4

    def getThreshold(self):
        try:
            self.geneNodeOverlapCntThreshold = int(self.ui.threshold.text().rstrip(' '))
            if self.geneNodeOverlapCntThreshold <= 0:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal value for min num. of node overlapping with a gene. Please check!', QtWidgets.QMessageBox.Ok)
                return -5               
            self.geneNodeOverlapCntThreshold = str(self.geneNodeOverlapCntThreshold)
        except ValueError:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal value for min num. of node overlapping with a gene. Please check!', QtWidgets.QMessageBox.Ok)
            return -5        

    def saveChanges(self):
        if self.getHeight() != -1 and self.getWidth() != -2 and self.getMaxNodesDisplay() != -3 and self.getMaxNodeLenDisplay() != -4 and self.getThreshold() != -5:
            self.Interaction()
            self.AddRemoveNode()
            self.LayoutModification()

            config = ConfigParser()
            config.optionxform = str
            config.read(copied_config)
            config.set("enable", "interaction", self.interaction)
            config.set("enable", "graphLayoutModification", self.layoutModification)
            config.set("enable", "addRemoveNodes", self.addRemoveNodes)
            config.set("canvas", "height", self.height)
            config.set("canvas", "width", self.width)
            config.set("nodes", "maxNodesDisplay", self.maxNodesDisplay)
            config.set("nodes", "maxNodeLenDisplay", self.maxNodeLenDisplay)
            config.set("nodes", "geneNodeOverlapCntThreshold", self.geneNodeOverlapCntThreshold)
            with open(copied_config, 'w') as configfile:
                config.write(configfile)
            QtWidgets.QMessageBox.question(self, 'Information', 'The settings of graph modification have been saved !',
            QtWidgets.QMessageBox.Ok)

    def reset(self):
        self.ui.interaction.clear()
        self.ui.interaction.addItems(['No', 'Yes'])
        self.interaction = 'No'
        self.ui.layoutModification.clear()
        self.ui.layoutModification.addItems(['No', 'Yes'])
        self.layoutModification = 'No'
        self.ui.addRemoveNodes.clear()
        self.ui.addRemoveNodes.addItems(['No', 'Yes'])
        self.addRemoveNodes = 'No'
        self.height = '980'
        self.width = '1400'
        self.maxNodesDisplay = '200'
        self.maxNodeLenDisplay = '1000000'
        self.geneNodeOverlapCntThreshold = '2'
        self.ui.height.setText(self.height)
        self.ui.width.setText(self.width)
        self.ui.maxNodesDisplay.setText(self.maxNodesDisplay)
        self.ui.nodeLen.setText(self.maxNodeLenDisplay)
        self.ui.threshold.setText(self.geneNodeOverlapCntThreshold)

        config = ConfigParser()
        config.optionxform = str
        config.read(copied_config)
        config.set("enable", "interaction", self.interaction)
        config.set("enable", "graphLayoutModification", self.layoutModification)
        config.set("enable", "addRemoveNodes", self.addRemoveNodes)
        config.set("canvas", "height", self.height)
        config.set("canvas", "width", self.width)
        config.set("nodes", "maxNodesDisplay", self.maxNodesDisplay)
        config.set("nodes", "maxNodeLenDisplay", self.maxNodeLenDisplay)
        config.set("nodes", "geneNodeOverlapCntThreshold", self.geneNodeOverlapCntThreshold)

        with open(copied_config, 'w') as configfile:
            config.write(configfile)
        QtWidgets.QMessageBox.question(self, 'Information', 'The settings of graph modification have been initialized !',
        QtWidgets.QMessageBox.Ok)

class Main(QMainWindow):
    def __init__(self, parent = None):
        QMainWindow.__init__(self, parent)
        self.ui = Ui_panGraphViewer()
        self.ui.setupUi(self)
        self.gfa = None
        self.fasta = None
        self.fastaMissing = 0
        self.vcf = None
        self.VCFresponse = None
        self.bed = None
        self.bedError = 0
        self.backbone = None
        self.chrs = list()
        self.chr = None
        self.start = None
        self.end = None
        self.shape = 'dot'
        self.enable = 'no'
        self.node = None
        self.nodes = list()
        self.selectedNodes = list()
        self.gene = None
        self.progress = None
        self.canvas = None
        self.outDir = None
        self.outDir1 = None
        self.run = None
        self.lines = list()
        self.event_stop = threading.Event()

        self.threadpool = QThreadPool()
        #==================== menubar ====================
        #shrotcuts
        self.ui.actionNew.setShortcut('Ctrl+N')
        self.ui.actionClose.setShortcut('Ctrl+Q')
        self.ui.actionAbout.setShortcut('Ctrl+A')
        self.ui.actionManual.setShortcut('ALT+M')
        self.ui.actionFeedback.setShortcut('Ctrl+F')
        self.ui.actionContact_the_Author.setShortcut('Ctrl+E')
        self.ui.actionModify_FASTA_header.setShortcut('ALT+H')
        self.ui.actionNode_Shapes.setShortcut('ALT+N')
        self.ui.actionGraph_modification.setShortcut('Ctrl+G')
        self.ui.actionNCBI_Blast.setShortcut('Ctrl+B')
        self.ui.actionMSA.setShortcut('Ctrl+L')
        self.ui.actionVCF_to_rGFA.setShortcut('ALT+V')
        self.ui.actionPrimary_Screen_Clipping.setShortcut('ALT+P')
        self.ui.actionFull_Screen_s_Capturing.setShortcut('ALT+C')
        self.ui.actionUpdate.setShortcut('Ctrl+U')

        #actions
        self.ui.actionNew.triggered.connect(self.new)
        self.ui.actionClose.triggered.connect(self.quit)
        self.ui.actionVCF_to_rGFA.triggered.connect(self.VCF2rGFA)
        self.ui.actionModify_FASTA_header.triggered.connect(self.modifyFaHeader)
        self.ui.actionNode_Shapes.triggered.connect(self.nodeShapeModify)
        self.ui.actionGraph_modification.triggered.connect(self.graphModify)
        self.ui.actionAbout.triggered.connect(self.about)
        self.ui.actionManual.triggered.connect(self.manual)
        self.ui.actionContact_the_Author.triggered.connect(self.sendEmail)
        self.ui.actionFeedback.triggered.connect(self.feedback)
        self.ui.actionNCBI_Blast.triggered.connect(self.blast)
        self.ui.actionMSA.triggered.connect(self.msa)
        self.ui.actionPrimary_Screen_Clipping.triggered.connect(self.clipping)
        self.ui.actionFull_Screen_s_Capturing.triggered.connect(self.screenshot)
        self.ui.actionUpdate.triggered.connect(self.update)

        #=================== load gfa panel ==================
        self.ui.selectGFA.clicked.connect(self.select_gfa)
        self.ui.clearGFA.clicked.connect(self.clear_gfa)
        self.ui.selectOutPath.clicked.connect(self.outDir_select)
        self.ui.clearOutPath.clicked.connect(self.outDir_clear)

        #=================== load vcf panel ==================
        self.ui.fastaSelect.clicked.connect(self.selectFasta)
        self.ui.fastaClear.clicked.connect(self.clearFasta)
        self.ui.vcfSelect.clicked.connect(self.selectVCF)
        self.ui.vcfClear.clicked.connect(self.clearVCF)
        self.ui.outDirSelect.clicked.connect(self.outdirSelect)
        self.ui.outDirClear.clicked.connect(self.outdirClear)
        

        #=================== parse gfa panel ====================
        self.ui.startParseGFA.clicked.connect(self.beginParseGAF)

        #=================== parse vcf panel ====================
        self.ui.parseVCF.clicked.connect(self.convertParseWorker)

        #=============== parameter setting panel =================
        #basic
        self.ui.comboBoxSample.activated[str].connect(self.Backbone)
        self.ui.comboBoxChr.activated[str].connect(self.Chr)
        #sample(s) display
        self.ui.removeSample.clicked.connect(self.removeSample)
        self.ui.AddBackSample.clicked.connect(self.addSample)
        self.ui.ResetSamples.clicked.connect(self.resetSample)

        #=================== plot panel ==================
        self.ui.plotGraph.clicked.connect(self.beginShowGraph)
        self.ui.clearGraphPlot.clicked.connect(self.plotClean)

        #=================== node panel ==================
        self.ui.addNode.clicked.connect(self.addNodes)
        self.ui.checkNodesInfo.clicked.connect(self.showSeqWorker)
        self.ui.rmSelectedNodes.clicked.connect(self.removeNode)
        self.ui.nodesComboBox.activated[str].connect(self.Node)
        self.ui.clearNodeList.clicked.connect(self.clearNodesList)
        self.ui.saveNodesInfo.clicked.connect(self.saveSeqWorker)

        #================== annotation panel ===================
        self.ui.selectBedPath.clicked.connect(self.select_bed)
        self.ui.clearBedPath.clicked.connect(self.clear_bed)
        self.ui.parseGene.clicked.connect(self.parsBedWorker)
        self.ui.plotOvlap.clicked.connect(self.plotOverlapWorker)

        ##Lock frames
        self.ui.ParameterFrame.setEnabled(False)
        self.ui.plotFrame.setEnabled(False)
        self.ui.checkNodesFrame.setEnabled(False)
        self.ui.checkGeneOvlpFrame.setEnabled(False)
        self.disableParameters()

        #self.ui.VCFtab.setEnabled(False)
        #self.ui.actionVCF_to_rGFA.setEnabled(False)

        self.ui.vizCanvas.tabCloseRequested.connect(lambda index: self.removeTab(index))

        self.canvas_link = []

    def removeTab(self, index):
        self.ui.vizCanvas.removeTab(index)

        canvas = self.canvas_link[index]
        canvas.setParent(None)
        canvas.deleteLater() 
        self.canvas_link.pop(index)

    #=============================================== functions ===============================================
    def new(self): # create a new GUI in different PID
        subprocess.Popen([sys.executable, tool])

    def quit(self): # quit the program
        response = QtWidgets.QMessageBox.question(self, 'Warning !', 'Do you really want to close the appliction?',
        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

        if response == QtWidgets.QMessageBox.Yes:
             sys.exit(0)

    def about(self): # A brief instruction for the program
        self.about_window = about()
        self.about_window.show()

    def manual(self): # A guideline for the program (use a web browser to open the Manual)
        webbrowser.open(f'file://{manualPath}')

    def feedback(self): # report the issue confronted
        issueAddress = "https://github.com/TF-Chan-Lab/panGraphViewer/issues"
        webbrowser.open(issueAddress)

    def sendEmail(self): # seems macOS cannot load this
        email = "yuxuan.yuan@outlook.com"
        webbrowser.open(f'mailto:{email}', new=1)

    def blast(self): # redirect to the web-based NCBI-BLAST
        blastAddress = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome"
        webbrowser.open(blastAddress)

    def msa(self): # redirect to the web-based MAFFT
        mafftAddress = "https://www.ebi.ac.uk/Tools/msa/mafft"
        webbrowser.open(mafftAddress)

    def VCF2rGFA(self): # load vcf2rGFA convertor
        self.vcf2rGFA_window = vcf_2_rGFA()
        self.vcf2rGFA_window.show()

    def modifyFaHeader(self): # modify fasta header
        self.modifyFaHeader_window = renameFaHeader()
        self.modifyFaHeader_window.show()

    def nodeShapeModify(self): # settings to change node shapes
        self.nodeShapeModify_window = NodeShapes()
        self.nodeShapeModify_window.show()

    def graphModify(self): # settings to enable graph modification
        self.graphModify_window = GraphModify()
        self.graphModify_window.show()

    def clipping(self): # allow primary screen capture
        if sys.platform == "win32":
            if self.outDir == None or self.outDir == '':
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Please specify an output directory first !',
                QtWidgets.QMessageBox.Ok)
                return -1
            now = datetime.datetime.now()
            name = os.path.join(self.outDir, now.strftime("%Y-%m-%d-%H_%M_%S"))
            screenshotScript = os.path.join(toolPath, 'scripts', 'screenClipping.py')
            subprocess.Popen([sys.executable, screenshotScript, name])
            #self.ClippingWindow = CaptureScreen()
            #self.ClippingWindow.show()
        else:
            if sys.platform == "darwin":
                QtWidgets.QMessageBox.question(self, 'Information !', "You may press 'Command+Shift+4' to capture and the 'png' figure will be saved to your Desktop !",
                QtWidgets.QMessageBox.Ok) 
            else:
                QtWidgets.QMessageBox.question(self, 'Information !', 'Oops! Currently, we only support this in Windows platforms !',
                QtWidgets.QMessageBox.Ok)

    def screenshot(self): # allow full screen capture 
        if sys.platform != "darwin":
            screenshotScript = os.path.join(toolPath, 'scripts', 'screenshot.py')
            subprocess.Popen([sys.executable, screenshotScript])
            #self.screenshotWindow = Screenshot()
            #self.screenshotWindow.show()
        else:
            QtWidgets.QMessageBox.question(self, 'Information !', "You may press 'Command+Shift+3' to capture and the 'png' figure will be saved to your Desktop !",
            QtWidgets.QMessageBox.Ok)            

    def update(self):
        response = QtWidgets.QMessageBox.question(self, 'Warning !', 'Do you want to open the GitHub release page to check new versions ?',
        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if response == QtWidgets.QMessageBox.Yes:
            relasePage = "https://github.com/TF-Chan-Lab/panGraphViewer/releases"
            webbrowser.open(relasePage)

    def select_gfa(self): # load GFA file
        self.ui.comboBoxChr.clear()
        self.ui.comboBoxSample.clear()
        self.ui.nodesComboBox.clear()
        if sys.platform == 'win32':
            gfa = QFileDialog.getOpenFileName(self, 'Select graphical fragment assembly', '','gfa (*gfa)')
            self.gfa=codecs.decode(str(gfa)[1:-1].split(',')[0][1:-1],'unicode_escape')
        else:
            self.gfa = QFileDialog.getOpenFileName(self, 'Select graphical fragment assembly', '','gfa (*gfa)')[0]
        self.gfa = os.path.abspath(self.gfa)

        if checkFile(self.gfa) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the input (r)GFA file. Please check !',
            QtWidgets.QMessageBox.Ok)
            self.ui.gfaPath.clear()
            self.gfa = None
        else:
            self.ui.gfaPath.setText(self.gfa)

    def clear_gfa(self): # clear the infomation of the loaded rGFA file
        self.ui.comboBoxChr.clear()
        self.ui.comboBoxSample.clear()
        self.ui.nodesComboBox.clear()
        self.ui.gfaPath.clear()
        self.gfa = None

    def selectFasta(self):  # load GFA file
        self.fasta = None 
        self.ui.fastaPath.clear()
        if sys.platform == 'win32':
            fasta = QFileDialog.getOpenFileName(self, 'Select the backbone genome assembly', '',
                                                'fasta(*fa *fasta *fna)')
            self.fasta = codecs.decode(str(fasta)[1:-1].split(',')[0][1:-1], 'unicode_escape')
        else:
            self.fasta = \
            QFileDialog.getOpenFileName(self, 'Select the backbone genome assembly', '', 'fasta(*fa *fasta *fna)')[0]
        self.fasta = os.path.abspath(self.fasta)
        
        if checkFile(self.fasta) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the input FASTA file. Please check !',
            QtWidgets.QMessageBox.Ok)
            self.ui.fastaPath.clear()
            self.fasta = None
        else:
            self.ui.fastaPath.setText(self.fasta)

    def clearFasta(self):
        self.ui.fastaPath.clear()
        self.fasta = None

    def selectVCF(self):
        self.vcf = None 
        self.ui.vcfPath.clear()
        if sys.platform == 'win32':
            vcf = QFileDialog.getOpenFileName(self, 'Select the vcf file', '', 'vcf(*vcf)')
            self.vcf = codecs.decode(str(vcf)[1:-1].split(',')[0][1:-1], 'unicode_escape')
        else:
            self.vcf = QFileDialog.getOpenFileName(self, 'Select the vcf file', '', 'vcf(*vcf)')[0]
        self.vcf = os.path.abspath(self.vcf)
        if checkFile(self.vcf) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the input VCF file. Please check !',
            QtWidgets.QMessageBox.Ok) 
            self.ui.vcfPath.clear()
            self.vcf = None
        else:
            self.ui.vcfPath.setText(self.vcf)           

    def clearVCF(self):
        self.ui.vcfPath.clear()
        self.vcf = None

    def outdirSelect(self):
        self.outDir1 = None 
        self.ui.outDir.clear()
        self.outDir1 = QFileDialog.getExistingDirectory()
        self.outDir1 = os.path.abspath(self.outDir1)
        if checkDir(self.outDir1) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the outDir. Please check !',
            QtWidgets.QMessageBox.Ok)
            self.ui.outDir.clear()
            self.outDir1 = None
        else:
            self.ui.outDir.setText(self.outDir1)

    def outdirClear(self):
        self.ui.outDir.clear()
        self.outDir1 = None

    def checkVCFparseElem(self):                
        if self.vcf == None or len(self.vcf) == 0:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select a VCF file',
            QtWidgets.QMessageBox.Ok)
            return -2
        if self.outDir1 == None or len(self.outDir1) == 0:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please specify an output directory',
            QtWidgets.QMessageBox.Ok)
            return -3
        if self.fasta == None or len(self.fasta) == 0:
            self.fastaMissing = 1
            self.VCFresponse = QtWidgets.QMessageBox.question(self, 'Warning !', 'While a fasta file is missing, we can still generate a graph. However, node sequence checking will not be available. Do you want to continue ?',
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            if self.VCFresponse == QtWidgets.QMessageBox.No:
                return -1
        if self.VCFresponse == QtWidgets.QMessageBox.Yes:
            self.ui.checkNodesFrame.setEnabled(False)
        else:
            self.ui.checkNodesFrame.setEnabled(True)
        self.VCFresponse = None
        return 0

    def getBB(self):
        self.bb = self.ui.bbName.text().rstrip(' ') 
 
    def getThreads(self):
        self.threads = self.ui.threads.text().rstrip(' ')

    def checkVcfElements(self):
        if self.checkVCFparseElem() == 0:
            self.bb = None
            self.getBB()
            if self.bb == '' or self.bb == None:
                self.bb = "backbone"
            self.getThreads()
            try:
                self.threads = int(self.threads)
            except ValueError:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Threads must be an integer !',
                                               QtWidgets.QMessageBox.Ok)
                return -2
            return 0

    def convertVCF(self):
        self.Date_vcf1 = datetime.datetime.today()
        self.ui.statusVCF.setStyleSheet('color: blue')
        self.ui.statusVCF.setText('Converting ...')
        self.ui.rGFAtab.setEnabled(False)
        self.ui.parseVCF.setEnabled(False)
        self.ui.ParameterFrame.setEnabled(False)
        self.ui.plotFrame.setEnabled(False)
        self.ui.checkNodesFrame.setEnabled(False)
        self.ui.checkGeneOvlpFrame.setEnabled(False)
        path = os.path.join(toolPath, 'scripts')
        script = os.path.join(path, 'vcf2rGFA.py')
        if self.fastaMissing == 0:
            cmd = f'"{sys.executable}" "{script}" -f "{self.fasta}" -v "{self.vcf}" -n "{self.threads}" -o "{self.outDir1}" -b "{self.bb}"'
        else:
            cmd = f'"{sys.executable}" "{script}" -v "{self.vcf}" -n "{self.threads}" -o "{self.outDir1}" -b "{self.bb}"'
            self.fastaMissing = 0
        if sys.platform == 'win32':
            self.retcode = subprocess.call(shlex.split(cmd), shell = True)      
        else:
            self.retcode = subprocess.call(shlex.split(cmd), shell = False)

        #args = options(self.vcf, self.fasta, self.bb, self.outDir, self.threads) ## convert the entire vcf 
        #self.retcode = VCF2rGFA(args).run()
        #self.ui.parseVCF.setEnabled(True)

    def endParse(self):
        if self.retcode == 1:
            QtWidgets.QMessageBox.question(self, 'Error !', "The input VCF file doesn't have a sample column. Abort !",
            QtWidgets.QMessageBox.Ok)                
            self.ui.statusVCF.setStyleSheet('color: red')
            self.ui.statusVCF.setText('Aborted !')
            self.ui.parseVCF.setEnabled(True)
            self.ui.rGFAtab.setEnabled(True)
            return -1

        if self.retcode == 2:
            QtWidgets.QMessageBox.question(self, 'Error !', 'Unmatched FASTA ! Please check the log info. Abort !',
            QtWidgets.QMessageBox.Ok)
            self.ui.statusVCF.setStyleSheet('color: red')
            self.ui.statusVCF.setText('Aborted !')
            self.ui.parseVCF.setEnabled(True)
            self.ui.rGFAtab.setEnabled(True)
            return -2

        if self.retcode == 3:
            QtWidgets.QMessageBox.question(self, 'Error !', 'Chromosome length info is missing in the VCF headers. Abort !',
            QtWidgets.QMessageBox.Ok)                
            self.ui.statusVCF.setStyleSheet('color: red')
            self.ui.statusVCF.setText('Aborted !')
            self.ui.rGFAtab.setEnabled(True)
            self.ui.parseVCF.setEnabled(True)
            return -3

        self.Date_vcf2 = datetime.datetime.today()
        Time_duration = (self.Date_vcf2 - self.Date_vcf1).total_seconds()

        self.gfa1 = os.path.join(f"{self.outDir1}", f"{self.bb}_vcf2rGFA.gfa")
        self.run = PanGraph(f"{self.gfa1}", f"{self.outDir1}")
        self.allChrs = self.run.backbone['contigs']
        self.displayInfo()
        self.resetSample()
        self.ui.statusVCF.setStyleSheet('color: green')
        self.ui.statusVCF.setText('Finished in %.2fs !' % Time_duration)

        self.ui.rGFAtab.setEnabled(True)
        self.ui.parseVCF.setEnabled(True)
        self.ui.startParseGFA.setEnabled(True)
        self.ui.plotGraph.setEnabled(True)
        self.ui.ParameterFrame.setEnabled(True)
        self.ui.plotFrame.setEnabled(True)
        if self.VCFresponse == QtWidgets.QMessageBox.Yes:
            self.ui.checkNodesFrame.setEnabled(False)
        else:
            self.ui.checkNodesFrame.setEnabled(True)
        self.ui.checkGeneOvlpFrame.setEnabled(True)

    def convertParseWorker(self):
        if self.checkVcfElements() == 0:
            self.workerVCF = Worker(self.convertVCF)
            self.threadpool.start(self.workerVCF)
            self.workerVCF.signals.finished.connect(self.endParse)
            self.workerVCF.signals.progress.connect(self.showProgress)

    def enableParameters(self): # enable the parameter panel
        self.ui.comboBoxSample.setEnabled(True)
        self.ui.comboBoxChr.setEnabled(True)
        self.ui.startPosition.setEnabled(True)
        self.ui.endPosition.setEnabled(True)
        self.ui.sampleComboBox.setEnabled(True)
        self.ui.removeSample.setEnabled(True)
        self.ui.AddBackSample.setEnabled(True)
        self.ui.ResetSamples.setEnabled(True)

    def disableParameters(self): # disable the parameter panel
        self.ui.comboBoxSample.setEnabled(False)
        self.ui.comboBoxChr.setEnabled(False)
        self.ui.startPosition.setEnabled(False)
        self.ui.endPosition.setEnabled(False)
        self.ui.sampleComboBox.setEnabled(False)
        self.ui.removeSample.setEnabled(False)
        self.ui.AddBackSample.setEnabled(False)
        self.ui.ResetSamples.setEnabled(False)

    def outDir_select(self): # set the output directory
        self.outDir = None 
        self.ui.outDirPath.clear()
        self.outDir = QFileDialog.getExistingDirectory()
        self.outDir = os.path.abspath(self.outDir)
        
        if checkDir(self.outDir)== -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the outDir. Please check !',
            QtWidgets.QMessageBox.Ok)
            self.ui.outDirPath.clear()
            self.outDir = None
        else:
            self.ui.outDirPath.setText(self.outDir)

    def outDir_clear(self): # clear the information of the selected outdir
        self.ui.outDirPath.clear()
        self.outDir = None

    def checkGFAmsg(self): # check if an rGFA file is loaded and pop out a dialog box
        if self.gfa == None or len(self.gfa) == 0:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select an rGFA or a GFA file',
            QtWidgets.QMessageBox.Ok)
            return -1
        else:
            return 0

    def checkOutDirPath(self): # check if the outDir is selected and pop out a dialog box
        noGFApath, noVCFpath = False, False
        if self.outDir == None or len(self.outDir) == 0:
            noGFApath = True
        if self.outDir1 == None or len(self.outDir1) == 0:
            noVCFpath = True

        if noGFApath == True and noVCFpath == True:        
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please specify an output directory',
            QtWidgets.QMessageBox.Ok)
            return -1
        else:
            return 0

    def completeGFAparse(self): # show the status of the parsing panel
        self.ui.status.setStyleSheet('color: green')
        self.date_gfa2 = datetime.datetime.today()
        time_delta = (self.date_gfa2 - self.date_gfa1)
        total_seconds = time_delta.total_seconds()
        self.ui.status.setText('Finished in %.2fs ! ' % total_seconds)
        self.allChrs = self.run.backbone['contigs']
        self.ui.VCFtab.setEnabled(True)
        self.displayInfo()
        self.enableFrames()
        
        self.resetSample()
        if self.run.neededGFA == 1:
            QtWidgets.QMessageBox.question(self, 'Warning !', "The format of 'SN:Z' in the input rGRF file is not the one we want. Although we can display the graph, many features would not show. Please refer to the 'Manual' in Help if you want a better experience in using this tool.", QtWidgets.QMessageBox.Ok)

    def showProgress(self):
        pass

    def showCompleted(self):
        pass

    def prepareGFAparse(self):
        self.ui.nodesComboBox.clear()
        self.ui.sampleComboBox.clear()
        self.date_gfa1 = datetime.datetime.today()
        self.disableFrames()
        self.ui.status.setStyleSheet('color: blue')
        self.ui.status.setText('Parsing ...')        
        self.disableParameters()
        #self.run = panGraph(self.gfa, self.outDir)
        self.run = PanGraph(self.gfa, self.outDir) 

    def disableFrames(self):
        self.ui.startParseGFA.setEnabled(False)
        self.ui.startParseGFA.setEnabled(False)
        self.ui.ParameterFrame.setEnabled(False)
        self.ui.plotFrame.setEnabled(False)
        self.ui.checkNodesFrame.setEnabled(False)
        self.ui.checkGeneOvlpFrame.setEnabled(False)

    def enableFrames(self):
        self.ui.startParseGFA.setEnabled(True)
        self.ui.startParseGFA.setEnabled(True)
        self.ui.ParameterFrame.setEnabled(True)
        self.ui.plotFrame.setEnabled(True)
        self.ui.checkNodesFrame.setEnabled(True)
        self.ui.checkGeneOvlpFrame.setEnabled(True)                 

    def fmtCheckPrep(self):
        path = os.path.join(toolPath, 'scripts')
        script = os.path.join(path, 'gfa2rGFA.py')
        file = os.path.split(self.gfa)[1]
        converted = os.path.join(self.outDir, f'{os.path.splitext(file)[0]}_converted{os.path.splitext(file)[1]}')
        cmd = f'"{sys.executable}" "{script}" -in_gfa "{self.gfa}" -out_rgfa "{converted}"'
        if sys.platform == 'win32':
            self.cvtcode = subprocess.call(shlex.split(cmd), shell = True)
        else:
            self.cvtcode = subprocess.call(shlex.split(cmd), shell = False)
        
        #print("The cvtcode is: '%s'" % self.cvtcode)

        if self.cvtcode not in [4, 5, 6, 7]:
            self.prepareGFAparse()

    def fmtcvtDone(self):
        try:
            if self.cvtcode == 4: 
                self.ui.status.setText('Crashed')
                self.enableFrames()
                QtWidgets.QMessageBox.question(self, 'Error!', 'Missing GFA field in the input file. Abort!',
                QtWidgets.QMessageBox.Ok)    
                return -2
            if self.cvtcode == 5:  
                self.ui.status.setText('Crashed')
                self.enableFrames()
                QtWidgets.QMessageBox.question(self, 'Error!', 'both seq and LN tag are missing in the input file. Abort!',
                QtWidgets.QMessageBox.Ok) 
                return -3 
            if self.cvtcode == 6: 
                self.ui.status.setText('Crashed')
                self.enableFrames()
                QtWidgets.QMessageBox.question(self, 'Error!', 'The file is in GFA v1, but error occurs during conversion. Abort!',
                QtWidgets.QMessageBox.Ok)  
                return -4
            if self.cvtcode == 7:
                self.ui.status.setText('Crashed')
                self.enableFrames()
                QtWidgets.QMessageBox.question(self, 'Error!', 'Path info is missing from the input GFA1 file, which is necessary for conversion. Abort!',
                QtWidgets.QMessageBox.Ok) 
                return -5
            else:
                self.ui.status.setText('Converted')  
                file = os.path.split(self.gfa)[1]   
                self.gfa = os.path.join(self.outDir, f'{os.path.splitext(file)[0]}_converted{os.path.splitext(file)[1]}') 
                self.ui.gfaPath.setText(self.gfa)                          
        except:
            pass 
               
        self.run = PanGraph(f"{self.gfa}", f"{self.outDir}")
        self.allChrs = self.run.backbone['contigs']
        self.displayInfo()
        self.resetSample()
        self.completeGFAparse()


    def BeginParseGAF(self): # start to parse the rGFA file
        if self.checkGFAmsg() == 0 and self.checkOutDirPath() == 0:
            self.ui.VCFtab.setEnabled(False)
            self.workerParseGFA = Worker(self.prepareGFAparse)
            self.threadpool.start(self.workerParseGFA)
            self.workerParseGFA.signals.finished.connect(self.completeGFAparse)
            self.workerParseGFA.signals.progress.connect(self.showProgress)



    def beginParseGAF(self): # start to parse the rGFA file
        if self.checkGFAmsg() == 0 and self.checkOutDirPath() == 0:
            #if self.formatCheck() != -1: 
            run = GFA2rGFA(self.gfa, self.outDir)
            if run.checkGfaFormat() == 1:
                self.BeginParseGAF()
            if run.checkGfaFormat() == 2:
                rect = QtWidgets.QMessageBox.question(self, 'Warring!', 'The input file is in GFA1 format and we need to do some checking and conversion. This will take a while, do you want to continue?',
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                if rect == QtWidgets.QMessageBox.Yes:
                    self.ui.status.setStyleSheet('color: blue')
                    self.ui.status.setText('Checking ...')
                    self.ui.VCFtab.setEnabled(False)
                    self.workerParseGFA = Worker(self.fmtCheckPrep)
                    self.threadpool.start(self.workerParseGFA)
                    self.workerParseGFA.signals.finished.connect(self.fmtcvtDone)
                    self.workerParseGFA.signals.progress.connect(self.showProgress)
            if run.checkGfaFormat() == 3:
                self.ui.status.setText('Crashed')
                QtWidgets.QMessageBox.question(self, 'Error!', 'The input file is in an unknown format. Abort!',
                QtWidgets.QMessageBox.Ok)                 

    def displayInfo(self): # enable and show nodes
        self.enableParameters()
        self.ui.comboBoxSample.clear()
        self.ui.sampleComboBox.clear()
        self.lines = list()
        try:
            self.lines = [i for i in self.run.nameCols.keys()]
            self.backboneList = [self.run.backbone['name']]
            self.leftSamples = self.lines
            if '' not in self.backboneList:
                self.backboneList.insert(0, '')
            self.ui.comboBoxSample.addItems(self.backboneList) 
            if '' not in self.leftSamples:
                self.leftSamples.insert(0, '')                      
            self.ui.sampleComboBox.addItems([i for i in self.leftSamples if i != self.run.backbone['name']])
            self.ui.nodesComboBox.clear()
            self.nodes = [i for i in self.run.inf['NodeID']]
            if '' not in self.nodes:
                self.nodes.insert(0, '')
            self.ui.nodesComboBox.addItems(self.nodes)
        except TypeError:
            pass

    def Backbone(self): # indicate the backbone sample and show the chrs
        self.backbone = self.ui.comboBoxSample.currentText()
        if self.backbone == None or self.backbone == '':
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please specify the name of the backbone sample!',
            QtWidgets.QMessageBox.Ok)
            return -1
        self.ui.comboBoxChr.clear()
        self.clearNodesList()
        self.chrs = list()
        
        self.chrs = natsorted(self.allChrs)
        if '' not in self.chrs:
            self.chrs.insert(0, '')
        self.ui.comboBoxChr.addItems(self.chrs)

    def Chr(self): # select chr
        self.chr = self.ui.comboBoxChr.currentText()
    
    def checkChr(self):
        if self.chr == list() or self.chr == '':
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select a chromosome name to plot!',
            QtWidgets.QMessageBox.Ok)
            return -1

    def checkStartEnd(self): # check the input coordinates
        self.start = self.ui.startPosition.text()
        try:
            if self.start not in [None, ''] and int(self.start) < 0:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal Start postion. Must >= 1. Please Check!',
                QtWidgets.QMessageBox.Ok)
                return -1
        except ValueError:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal start position. Must be an integer. Please Check!',
            QtWidgets.QMessageBox.Ok)
            return -1        
        self.end = self.ui.endPosition.text()
        if self.end not in [None, '']: 
            try:
                if self.start != None and self.start != '':
                    if int(self.start) > int(self.end):
                        QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal End postion. Must be bigger than the Start Position. Please Check!',
                        QtWidgets.QMessageBox.Ok)
                        return -1
            except ValueError:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Illegal End position. Must be an integer. Please Check!',
                QtWidgets.QMessageBox.Ok)
                return -1

    def beginShowGraph(self): # compute and show the graph
        self.backbone = self.ui.comboBoxSample.currentText()
        if self.backbone == '' or self.backbone == None:
            self.ui.plotGraph.setEnabled(True)
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please specify a backbone sample',
            QtWidgets.QMessageBox.Ok)
        else:
            self.Chr()
            if self.checkOutDirPath() != -1 and self.checkChr() != -1:
                if self.checkStartEnd() != -1:
                    self.ui.plotGraph.setEnabled(False)
                    self.ui.checkNodesInfo.setEnabled(False)
                    self.start = None if self.ui.startPosition.text() == '' else int(self.ui.startPosition.text())
                    self.end = None if self.ui.endPosition.text() == '' else int(self.ui.endPosition.text())
                    self.ui.plotStatusLabel.setStyleSheet('color: blue')
                    self.ui.plotStatusLabel.setText('Plotting ...')
                    self.date_plot1 = datetime.datetime.today()
                    self.leftSamples = [i for i in self.leftSamples if i != '']
                    self.workerShowGraph = Worker(self.run.drawGraph, self.leftSamples, self.chr, self.start, self.end)
                    self.threadpool.start(self.workerShowGraph)
                    self.workerShowGraph.signals.finished.connect(self.completePlot)
                    self.workerShowGraph.signals.progress.connect(self.showProgress)                  

    def completePlot(self): # show the status of plotting
        self.date_plot2 = datetime.datetime.today()
        time_delta = (self.date_plot2 - self.date_plot1)
        total_seconds = time_delta.total_seconds()

        if self.run.plotErrorSignal == 0:
            if self.run.emptyGraphSignal == 1:
                QtWidgets.QMessageBox.question(self, 'Error !', "There is no node left in the selected region. Abort !",
                QtWidgets.QMessageBox.Ok)
                self.ui.plotStatusLabel.setStyleSheet('color: red')
                self.ui.plotStatusLabel.setText('Abort !')
                self.ui.plotGraph.setEnabled(True)
                self.run.emptyGraphSignal = 0
                return -1
            if self.run.overNodeLimitSignal == 1:
                subNodesCount = len(self.run.subNodes)
                if subNodesCount > self.run.maxNodesLimit:
                    response = QtWidgets.QMessageBox.question(self, 'Warning !', f"The num. of nodes in the selected region is '{subNodesCount}', which is bigger than the limit '{self.run.maxNodesLimit}' in settings. We are afraid that we may not be able to render such a large graph at the moment. Please modify Start/End Position to limit the number of nodes", QtWidgets.QMessageBox.Ok)
                else:
                    response = QtWidgets.QMessageBox.question(self, 'Warning !', f"The num. of nodes in the selected region is '{subNodesCount}', which is bigger than the max. num. of node setting '{self.run.maxNodesDisplay}' in a vis.js plot. We will plot a Cytoscape graph instead of a Vis graph to ensure a better visualization. Do you want to continue?",
                    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                self.run.overNodeLimitSignal = 0
                if response == QtWidgets.QMessageBox.No or response == QtWidgets.QMessageBox.Ok:
                    self.ui.plotStatusLabel.setStyleSheet('color: red')
                    self.ui.plotStatusLabel.setText('Cancelled')
                    self.ui.plotGraph.setEnabled(True)
                    self.ui.checkNodesInfo.setEnabled(True)
                    return -2
            canvas = QtWebEngineWidgets.QWebEngineView()
            if self.start == None and self.end == None:
                i = self.ui.vizCanvas.addTab(canvas, f"{self.chr}")
            if self.start != None and self.end != None:
                i = self.ui.vizCanvas.addTab(canvas, f"{self.chr}: {self.start} - {self.end}")
            if self.start != None and self.end == None:
                i = self.ui.vizCanvas.addTab(canvas, f"{self.chr}: {self.start} -")
            if self.start == None and self.end != None:
                i = self.ui.vizCanvas.addTab(canvas, f"{self.chr}: 1 - {self.end}")
            
            self.canvas_link.append(canvas)
            self.ui.vizCanvas.setCurrentIndex(i)

            html = self.run.drawGraphResult['outHtml']

            canvas.load(QtCore.QUrl().fromLocalFile(html))
            self.ui.plotStatusLabel.setStyleSheet('color: green')
            self.ui.plotStatusLabel.setText('Finished in %.2fs ! ' % total_seconds)
            sleep(1)
            canvas.show()
            sleep(1)
        else:
            self.ui.plotStatusLabel.setStyleSheet('color: red')
            self.ui.plotStatusLabel.setText('Crashed !')
            sleep(1.5)
            if self.run.plotErrorSignal != -4:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'No nodes found in the selected chromosomal region. Please Check!',
                QtWidgets.QMessageBox.Ok)
            else:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Maximum recurison depth exceeded in the selected chromosomal region!',
                QtWidgets.QMessageBox.Ok)
        self.ui.plotGraph.setEnabled(True)
        self.ui.checkNodesInfo.setEnabled(True)

    def plotClean(self): # Clean the display canvas
        self.ui.plotStatusLabel.setStyleSheet('color: black')
        self.ui.plotStatusLabel.setText('Idle')
        while self.ui.vizCanvas.count() > 0:
            self.ui.vizCanvas.removeTab(0)

    def Node(self): # select a node
        self.node = self.ui.nodesComboBox.currentText()

    def manualNodes(self): # manually input nodes
        tmpList=[]
        manualInput = self.ui.textEdit.toPlainText().split('\n')
        for i in manualInput:
            i = ''.join(i.split())
            if i !='' and not i.startswith("List") :
                tmpList.append(i)
                if i not in self.nodes:
                    QtWidgets.QMessageBox.question(self, 'Warning !', f"Illegal node ID: '{i}' in the list. Please check!",
                    QtWidgets.QMessageBox.Ok)
                    self.selectedNodes = []
                    return -1
        self.selectedNodes = tmpList

    def addNodes(self): # add a node to the NodeList
        #if self.manualNodes() != -1:
        self.Node()
        if self.node == None or self.node == '':
            QtWidgets.QMessageBox.question(self, 'Warning !', 'No node is selected. Please Check!',
            QtWidgets.QMessageBox.Ok)
            return 
        else:
            if self.node in self.selectedNodes:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'The node selected is in the list. Please Check!',
                QtWidgets.QMessageBox.Ok)
                return 
            else:
                info = ''
                self.selectedNodes.append(self.node)
                #if self.manualNodes() != -1:
                for i in self.selectedNodes:
                    if info != '':
                        info = info + '\n'+ i
                    else:
                        info = i
                self.ui.textEdit.setPlainText("%s" % info)
                #else:  
                #    return 

    def removeNode(self): # remove a node in the NodeList
        if self.selectedNodes == list():
            QtWidgets.QMessageBox.question(self, 'Warning !', 'There is no node selected/input in this list!',
            QtWidgets.QMessageBox.Ok)
            return 
        else:
            self.Node()
            if self.node not in self.selectedNodes:
                QtWidgets.QMessageBox.question(self, 'Warning !', f"No node: '{self.node}' in the list to be removed. Please check!",
                QtWidgets.QMessageBox.Ok)
                return 
            else:
                #if self.manualNodes() != -1:
                info = ''
                tmp = [i for i in self.selectedNodes if i != self.node]
                self.selectedNodes = tmp
                for i in self.selectedNodes:
                    if info != '':
                        info = info + '\n'+ i
                    else:
                        info = i
                self.ui.textEdit.setPlainText(f'{info}')
                #else:
                #    return 

    def clearNodesList(self): # empty the NodeList
        self.selectedNodes = list()
        self.ui.textEdit.setPlainText("List of the selected node(s)")

    def getNodesInfo(self): # dialog box to show if a node in the NodeList
        if self.manualNodes() != -1:
            if self.backbone == None or self.backbone == '':
                QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select a backbone sample in the basic settings first!',
                QtWidgets.QMessageBox.Ok)
                return -1   
        else:
            return -2
        if self.selectedNodes == list() or self.selectedNodes == '':
            QtWidgets.QMessageBox.question(self, 'Warning !', 'There is no node selected/input in this list!',
            QtWidgets.QMessageBox.Ok)
            return -3
        self.ui.nodeStatus.setStyleSheet('color: blue')
        self.ui.nodeStatus.setText('Running ...')
        self.date_node1 = datetime.datetime.today()               
        return 1

    def getSequences(self): # show the sequence information of nodes
        self.manualNodes()
        self.run.searchNodes(self.selectedNodes)
    
    def showSequences(self):
        self.showInfo = NodesInfoBrowser()
        if self.run.bigNodeList == []:
            self.showInfo.ui.NodeInfotextBrowser.setPlainText(f'{self.run.nodesInfo}')
            self.showInfo.show()
        else:
            response = QtWidgets.QMessageBox.question(self, "Warning !", f"The length of node {self.run.bigNodeList} is greater than '{self.run.bigNodeSeqLen} bp' in settings, it may frozen the program when displaying. Do you want to continue ?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            if response == QtWidgets.QMessageBox.Yes:
                self.showInfo.ui.NodeInfotextBrowser.setPlainText(f'{self.run.nodesInfo}')
                self.showInfo.show()
            self.run.bigNodeList = []
        self.ui.nodeStatus.setStyleSheet('color: green')
        self.date_node2 = datetime.datetime.today()
        time_delta = (self.date_node2 - self.date_node1).total_seconds()
        self.ui.nodeStatus.setText('Finished in %.2fs ! ' % time_delta)
        self.run.nodesInfo == None

    def showSeqWorker(self):
        if self.getNodesInfo() == 1:
            self.workerShowSeq = Worker(self.getSequences)
            self.threadpool.start(self.workerShowSeq)
            self.workerShowSeq.signals.finished.connect(self.showSequences)
            self.workerShowSeq.signals.progress.connect(self.showProgress)

    def saveSequences(self): # save the selected node information to outDir
        if self.outDir == None:
            QtWidgets.QMessageBox.question(self, 'Error !', 'Please select an output directory!',
            QtWidgets.QMessageBox.Ok)
            self.ui.nodeStatus.setStyleSheet('color: red')
            self.ui.nodeStatus.setText('Aborted ! ')

        else:
            now = datetime.datetime.now()
            self.run.saveNodesInfo(self.selectedNodes, now.strftime("%Y-%m-%d-%H_%M_%S"))
            self.ui.nodeStatus.setStyleSheet('color: green')
            self.date_node4 = datetime.datetime.today()
            time_delta = (self.date_node4 - self.date_node1).total_seconds()
            self.ui.nodeStatus.setText('Finished in %.2fs ! ' % time_delta)
            self.run.nodesInfo = None
            QtWidgets.QMessageBox.question(self, 'Information', 'The sequences of seleted nodes have been saved !',
            QtWidgets.QMessageBox.Ok)

    def saveSeqWorker(self):
        if self.getNodesInfo() == 1:
            self.workerSaveSeq = Worker(self.getSequences)
            self.threadpool.start(self.workerSaveSeq)
            self.workerSaveSeq.signals.finished.connect(self.saveSequences)
            self.workerSaveSeq.signals.progress.connect(self.showProgress)        

    def removeSample(self): # remove a sample from the sampleList
        if self.backbone == None or self.backbone == '':
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select a backbone sample in the basic settings first!',
            QtWidgets.QMessageBox.Ok)
        else:
            sample = self.ui.sampleComboBox.currentText()
            if sample != '' and sample == self.backbone:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'The backbone sample cannot be removed!',
                QtWidgets.QMessageBox.Ok)
            else:
                if sample == '':
                    QtWidgets.QMessageBox.question(self, 'Warning !', 'There is no sample selected to be removed!',
                    QtWidgets.QMessageBox.Ok)
                else:

                    if sample not in self.leftSamples:
                        QtWidgets.QMessageBox.question(self, 'Warning !', 'There is no such a sample selected to be removed!',
                        QtWidgets.QMessageBox.Ok)
                    else:
                        self.leftSamples = [i for i in self.leftSamples if i != sample]
                        samples = 'The remaining samples are:\n'
                        for i in self.leftSamples:
                            if i != '':
                                samples = samples + i + '\n'
                        self.ui.textEdit_2.setPlainText(f'{samples}')

    def addSample(self): # add a sample to the sampleList
        sample = self.ui.sampleComboBox.currentText()
        if sample == '':
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select one sample to be added!',
            QtWidgets.QMessageBox.Ok)
        else:
            if sample in self.leftSamples:
                QtWidgets.QMessageBox.question(self, 'Warning !', 'The sample is in the list. Please check!',
                QtWidgets.QMessageBox.Ok)
            else:
                self.leftSamples.append(sample)
                samples = 'The remaining samples are:\n'
                for i in self.leftSamples:
                    if i != '':
                        samples = samples + i + '\n'
                self.ui.textEdit_2.setPlainText(f'{samples}')

    def resetSample(self): # reset the sampleList
        samples = 'The remaining samples are:\n'
        for i in self.lines:
            if i != '':
                samples = samples + i + '\n'
        self.ui.textEdit_2.setPlainText(f'{samples}')
        self.leftSamples = self.lines

    def select_bed(self): # select a bed/gff/gtf file
        self.ui.geneComboBox.clear()
        if sys.platform == 'win32':
            bed = QFileDialog.getOpenFileName(self, 'Select an annotation file', '','annotation (*bed *gtf *gff3)')
            self.bed = codecs.decode(str(bed)[1:-1].split(',')[0][1:-1],'unicode_escape')
        else:
            self.bed = QFileDialog.getOpenFileName(self, 'Select an annoation file', '','annoation (*bed *gtf *gff3)')[0]
        self.bed = os.path.abspath(self.bed)
    
        if checkFile(self.bed) == -1:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Something is wrong with the input annotation file. Please check !',
            QtWidgets.QMessageBox.Ok)
            self.ui.bedPath.clear()
            self.bed = None
        else:
            self.ui.bedPath.setText(self.bed)

    def clear_bed(self): # clear the information of the selected bed
        self.ui.bedPath.clear()
        self.ui.geneComboBox.clear()
        self.bed = None
        self.ui.bedParseStatus.setStyleSheet('color: black')
        self.ui.bedParseStatus.setText('Idle')

    def checkBed(self):
        if self.bed == None or len(self.bed) == 0:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select an annotation file',
            QtWidgets.QMessageBox.Ok)
            return -1
        return 0

    def parseBed(self):
        self.ui.bedParseStatus.setStyleSheet('color: blue')
        self.ui.bedParseStatus.setText('Running ...')
        self.date_bed1 = datetime.datetime.today()
        self.run.loadBedGff(self.bed)
        self.threshold = self.run.getGeneNodeOverlapThreadhold()
        self.run.loadRGFA()
        self.run.genGraph()
        self.run.updateNodes()
        self.genes = list()
        genes = self.run.nodeGeneOverlap(self.bed, self.threshold)
        if len(genes) == 0:
            self.bedErr=1
            logging.error(f"No gene found overlapping with at least '{self.threshold}' nodes !")
            #QtWidgets.QMessageBox.question(self, 'Warning !', f"No gene found overlapping with at least '{threshold}' nodes !",
            #QtWidgets.QMessageBox.Ok)
        else:
            now = datetime.datetime.now()
            geneFile = os.path.join(self.outDir, f'{now.strftime("%Y-%m-%d-%H_%M_%S")}.geneList.txt')
            with open(geneFile, 'w') as f:
                for i in sorted(genes.keys()):
                    self.genes.append(i)
                    f.write('%s\n' % i)
            if '' not in self.genes:
                self.genes.insert(0, '')
            self.ui.geneComboBox.addItems(self.genes)
    
    def completeBedParsing(self):
        self.date_bed2 = datetime.datetime.today()
        time_delta = (self.date_bed2 - self.date_bed1)
        total_seconds = time_delta.total_seconds()
        self.ui.bedParseStatus.setStyleSheet('color: green')
        self.ui.bedParseStatus.setText('Finished in %.2fs !' % total_seconds)
        try: 
            if self.bedErr == 1:
                self.bedErr = 0
                self.ui.geneComboBox.clear()
                QtWidgets.QMessageBox.question(self, 'Warning !', f"No gene found overlapping with at least '{self.threshold}' nodes !",
                QtWidgets.QMessageBox.Ok)
        except AttributeError:
            pass

    def parsBedWorker(self):
        if self.checkBed() == 0:
            self.workerParseBed = Worker(self.parseBed)
            self.threadpool.start(self.workerParseBed)
            self.workerParseBed.signals.finished.connect(self.completeBedParsing)
            self.workerParseBed.signals.progress.connect(self.showProgress)

    def checkplotOverlapCondition(self):
        self.gene = self.ui.geneComboBox.currentText()
        if self.gene == None or len(self.gene) == 0:
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please select a gene to check overlap with nodes !',
            QtWidgets.QMessageBox.Ok)
            return -1           
        self.backbone = self.ui.comboBoxSample.currentText()
        if self.backbone == '' or self.backbone == None:
            self.ui.plotGraph.setEnabled(True)
            QtWidgets.QMessageBox.question(self, 'Warning !', 'Please specify a backbone sample',
            QtWidgets.QMessageBox.Ok)
            return -2
        return 0
     
    def plotOverlap(self):
        self.date_overlap1 = datetime.datetime.today()
        self.ui.bedParseStatus.setStyleSheet('color: blue')
        self.ui.bedParseStatus.setText('Plotting ...')            
        try:
            self.run.overlapGenes(self.gene)
        except KeyError:
            self.bedError = 1

    def completePlottingOverlap(self):
        if self.bedError == 1:
            QtWidgets.QMessageBox.question(self, 'Error !', 'Please check if the annotation file is the correct one.',
            QtWidgets.QMessageBox.Ok)
            self.ui.bedParseStatus.setStyleSheet('color: red')
            self.ui.bedParseStatus.setText('Plotting crashed !') 
            self.ui.geneComboBox.clear()
            self.bedError = 0
            return -2             
        if self.run.noOverlap == True:
            QtWidgets.QMessageBox.question(self, 'Warning !', f'No node overlap with gene: {self.gene}',
            QtWidgets.QMessageBox.Ok)
            self.run.noOverlap = False
            self.ui.bedParseStatus.setStyleSheet('color: red')
            self.ui.bedParseStatus.setText('Plotting crashed !') 
            return -3  

        self.run.drawOverlapGenes(self.gene)
        sleep(1)
        self.date_overlap2 = datetime.datetime.today()
        time_delta = (self.date_overlap2 - self.date_overlap1)
        total_seconds = time_delta.total_seconds()
        self.ui.bedParseStatus.setStyleSheet('color: green')
        self.ui.bedParseStatus.setText('Finished in %.2fs !' % total_seconds)
        canvas = QtWebEngineWidgets.QWebEngineView()
        i = self.ui.vizCanvas.addTab(canvas,f"Overlap with Gene: {self.gene}")
        self.ui.vizCanvas.setCurrentIndex(i)
        html = os.path.join(self.outDir, f'drawOverlap_with_Gene-{self.gene}.html')
        canvas.load(QtCore.QUrl().fromLocalFile(html))
        canvas.show()
        sleep(1)

        canvas = QtWebEngineWidgets.QWebEngineView()
        chr = self.run.bed[self.gene]['Chr']
        start = self.run.bed[self.gene]['Start']
        end = self.run.bed[self.gene]['End']
        i = self.ui.vizCanvas.addTab(canvas,f"Subgraph within the region of Gene: {self.gene}")
        self.ui.vizCanvas.setCurrentIndex(i)
        self.run.drawGraph(self.leftSamples, chr, start, end)
        html = self.run.drawGraphResult['outHtml']
        canvas.load(QtCore.QUrl().fromLocalFile(html))
        canvas.show()

    def plotOverlapWorker(self):
        if self.checkplotOverlapCondition() == 0:
            self.workerPlotOverlap = Worker(self.plotOverlap)
            self.threadpool.start(self.workerPlotOverlap)
            self.workerPlotOverlap.signals.finished.connect(self.completePlottingOverlap)
            self.workerPlotOverlap.signals.progress.connect(self.showProgress)


if __name__=="__main__":
    app = QApplication([])
    window = Main()
    window.show()
    sys.exit(app.exec_())
