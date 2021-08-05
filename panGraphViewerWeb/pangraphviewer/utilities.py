import subprocess
import os
import sys
import shlex
import logging

#============================= Function =================================
##logging info
DEBUG="" #change it when debugging
logFormat = "%(asctime)s [%(levelname)s] %(message)s"
level = "DEBUG" if DEBUG != "" else "INFO"
logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )
#========================================================================

def run_shell_cmd(cmd, stdinInput=None):
    cmd_list = [str.strip('"').strip("'") for str in shlex.split(cmd, posix=False)]

    ret = subprocess.run(cmd_list, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, text=True)

    if ret.returncode:
        logging.error(f'cmd: {cmd}')
        logging.error(f'stderr: {ret.stderr.strip()}')

    return {'returncode':ret.returncode, 'out':ret.stdout.strip().split('\n'), 'err':ret.stderr.strip().split('\n')}

def convert_vcf_to_gfa(vcf, chr, backbone, outdir, prefix, fasta=None):
    pyscript = 'pangraphviewer/vcf2rGFA.py'
    cmd = f"'{sys.executable}' '{pyscript}' -b '{backbone}' -v '{vcf}' -o '{outdir}' -p '{prefix}'"
    if chr: cmd += f" -c '{chr}'"
    if fasta: cmd += f" -f '{fasta}'"

    result = run_shell_cmd(cmd)

    return result

def makedirs(dir):
    try:
        umask_orig = os.umask(0)
        os.makedirs(dir, mode=0o777, exist_ok=True)
    except:
        raise Exception(f'Cannot create directory {dir}')
    finally:
        os.umask(umask_orig)
