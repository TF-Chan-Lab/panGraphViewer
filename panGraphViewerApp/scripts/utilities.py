#!/usr/bin/env python3

import os
import sys
import logging
from subprocess import Popen, PIPE

#============================= Function =================================
##logging info
DEBUG="" #change it when debugging
logFormat = "%(asctime)s [%(levelname)s] %(message)s"
level = "DEBUG" if DEBUG != "" else "INFO"
logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )
#========================================================================


def checkDir(Dirname):
    logging.info("Checking folder: '%s'" % Dirname)
    dirname = os.path.abspath(Dirname)
    if not os.path.isdir(dirname):
        logging.error("Oops! Folder: '%s' does not exit. Please check!" % Dirname)
        return -1
    if not os.access(dirname, os.W_OK):
        logging.error("Oops! Folder: '%s' is not writable. Please check!" % Dirname)
        return -1

def checkFile(Filename):
    logging.info("Checking file: '%s'" % Filename)
    filename = os.path.abspath(Filename)
    if not os.path.isfile(filename):
        logging.error("Oops! File: '%s' does not exit. Please check!" % Filename)
        return -1
    if not os.access(filename, os.R_OK):
        logging.error("Oops! File: '%s' is not readable. Please check!" % Filename)
        return -1

def checkInt(value):
    try:
        val = int(value)
    except ValueError:
        logging.error("Oops! The input value '%s' is not an integer. Please check!" % value)
        return -1

class Utilities:
    @staticmethod
    def runCmdLine(cmd, input=None, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE):
        p = Popen(cmd, shell=shell, stdin=stdin, stdout=stdout, stderr=stderr, universal_newlines=True)
        out, _err = p.communicate(input=input)

        output = out.strip().split('\n')

        return output

    @staticmethod
    def get_var(cfg, group, var, must_have=False):
        global config

        if not config:
            config = ConfigParser()
            config.read(cfg)

        try:
            return config.get(group, var)
        except:
            if must_have:
                raise
            else:
                logging.warning(f'Config value [{group}][{var}] not found')

        return None


if __name__=="__main__":
    #print(Utilities.runCmdLine('pwd'))
    print(Utilities.get_var('nodes','SN_delim'))
