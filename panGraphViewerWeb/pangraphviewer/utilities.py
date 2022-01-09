import subprocess
import os
import sys
import shlex
import logging

from configparser import ConfigParser

#============================= Function =================================
##logging info
DEBUG="" #change it when debugging
logFormat = "%(asctime)s [%(levelname)s] %(message)s"
level = "DEBUG" if DEBUG != "" else "INFO"
logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )
#========================================================================

#config = None
script_directory = os.path.dirname(os.path.realpath(__file__))
copied = os.path.join(script_directory, '..', "config.ini")

def getVar(cfg, group, var, mustHave=False, forceRead=True):
    #global config
    config = None

    if not config or forceRead:
        config = ConfigParser()
        config.read(cfg)

    try:
        return config.get(group, var)
    except:
        if mustHave:
            raise
        else:
            logging.info(f'Config value [{group}][{var}] not found')
        return None

def rev_comp(seq):
    trans = str.maketrans('ACGTN*', 'TGCAN*')
    return seq.translate(trans)[::-1]

def run_shell_cmd(cmd, stdinInput=None, logError=True):
    cmd_list = [str.strip('"').strip("'") for str in shlex.split(cmd, posix=False)]

    ret = subprocess.run(cmd_list, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, text=True)

    if ret.returncode and logError:
        logging.error(f'cmd: {cmd}')
        logging.error(f'stderr: {ret.stderr.strip()}')

    return {'returncode':ret.returncode, 'out':ret.stdout.strip().split('\n'), 'err':ret.stderr.strip().split('\n')}

def convert_vcf_to_gfa(vcf, chr, backbone, outdir, prefix, fasta=None):
    pyscript = 'pangraphviewer/vcf2rGFA.py'
    cmd = f"'{sys.executable}' '{pyscript}' -b '{backbone}' -v '{vcf}' -o '{outdir}' -p '{prefix}'"
    if chr: cmd += f" -c '{chr}'"
    if fasta: cmd += f" -f '{fasta}'"

    return run_shell_cmd(cmd)

def convert_gfa_to_rgfa(gfa, rgfa):
    pyscript = 'pangraphviewer/gfa2rGFA.py'
    cmd = f"'{sys.executable}' '{pyscript}' -in_gfa '{gfa}' -out_rgfa '{rgfa}'"

    return run_shell_cmd(cmd, logError=False)

    #return subprocess.call(shlex.split(cmd), shell = False)

def makedirs(dir):
    try:
        umask_orig = os.umask(0)
        os.makedirs(dir, mode=0o777, exist_ok=True)
    except:
        raise Exception(f'Cannot create directory {dir}')
    finally:
        os.umask(umask_orig)
