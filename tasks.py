#!/usr/bin/env python

import celery
import os
import shutil
import subprocess
import sys
from time import sleep

from Bio import SeqIO

def run_cmd(cmd, ignore_error=False):
    """
    Run a command line command
    Returns True or False based on the exit code
    """
    proc = subprocess.Popen(cmd,shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
    out = proc.communicate()
    return_code = proc.returncode

    t = open('log.txt', 'w')
    e = open('log.err', 'w')

    t.write('%s\n'%str(out[0]))
    if return_code != 0 and not ignore_error:
        e.write('Command (%s) failed w/ error %d\n'
                        %(cmd, return_code))
        e.write('%s\n'%str(out[1]))
        e.write('\n')

    return bool(not return_code)

@celery.task(name='tasks.run_contiguator')
def run_contiguator(wdir, contigs, refs, evalue=1e-20, contiglength=1000,
        contigcoverage=20, hitlength=1100, multitreshold=1.5,
        non=False, numN=100, pcr=False, inner=False, blastn=False,
        threads=1,
        optsize=20,
        minsize=18,
        maxsize=27,
        opttemp=60,
        mintemp=57,
        maxtemp=63,
        flanksize=1000,
        minprod=1000,
        maxprod=7000,
        optgc=50,
        mingc=20,
        maxgc=80,
        gcclamp=1,
        exclude=100,
        jobname=''):
    sdir = os.getcwd()

    # Move all the contiguator files
    shutil.copy(os.path.join(sdir, 'contiguator-app', 'CONTIGuator.py'), wdir)
    shutil.copy(os.path.join(sdir, 'contiguator-app', 'abacas'), wdir)
 
    # Move to working directory
    os.chdir(wdir)

    # Program version
    try:
        vcmd = 'python CONTIGuator.py --version'
        proc = subprocess.Popen(vcmd,shell=(sys.platform!="win32"),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                cwd=wdir)
        out = proc.communicate()
        sleep(1)
        ver = out[0].strip() + out[1].strip()
    except:
        ver = 'Unknown'

    #Run CONTIGuator
    cmd = 'python CONTIGuator.py -c "%s"'%(os.path.basename(contigs))
    for ref in refs:
        cmd += ' -r "%s"'%(os.path.basename(ref))
    cmd += ' -e %s'%str(evalue)
    cmd += ' -L %d'%int(contiglength)
    cmd += ' -C %d'%int(contigcoverage)
    cmd += ' -B %d'%int(hitlength)
    cmd += ' -R %s'%str(multitreshold)
    cmd += ' -n %d'%int(numN)
    # Check threads
    if int(threads) > 3:
        threads = 3
    cmd += ' -t %d'%int(threads)
    if bool(blastn):
        cmd += ' -b'
    if bool(non):
        cmd += ' -N'
    if bool(pcr):
        cmd += ' -P -A'
        if bool(inner) == True:
            cmd += ' -I'
        # Other PCR commands written to file
        f=open(os.path.join(wdir,'newlines'),'w')
        f.write('%d\n%d\n%d\n%f\n%f\n%f\n%d\n%d\n%d\n%f\n%f\n%f\n%d\n%d\n'%(
                            int(optsize),
                            int(minsize),
                            int(maxsize),
                            float(opttemp),
                            float(mintemp),
                            float(maxtemp),
                            int(flanksize),
                            int(minprod),
                            int(maxprod),
                            float(optgc),
                            float(mingc),
                            float(maxgc),
                            int(gcclamp),
                            int(exclude)))
        f.close()
    cmd += ' -M -V -D -G'


    # Job details files
    #try:
    lgenparams = ['Job name', 'CONTIGuator command', 'Version']
    genparams = {'Job name':jobname,
                'CONTIGuator command':cmd,
                'Version':ver}
    fout = open(os.path.join(wdir,'genparams.tsv'), 'w')
    for k in lgenparams:
        fout.write('%s\t%s\n'%(k,genparams[k]))
    fout.close()
                
    lrunparams = ['Blast e-value', 'Use blastn', 'Blast threads',
                  'Contig length threshold',
                  'Contig coverage threshold (%)',
                  'Hit length threshold',
                  'Multiple replicon threshold',
                  'Gaps size on overlapping contigs',
                  'Do not use N to separate the contigs']
    runparams = {'Blast e-value':evalue,
                 'Use blastn':blastn,
                 'Contig length threshold':contiglength,
                 'Contig coverage threshold (%)':contigcoverage,
                 'Hit length threshold':hitlength,
                 'Multiple replicon threshold':multitreshold,
                 'Do not use N to separate the contigs':non,
                 'Gaps size on overlapping contigs':numN,
                 'Blast threads':threads}
    fout = open(os.path.join(wdir,'runparams.tsv'), 'w')
    for k in lrunparams:
        fout.write('%s\t%s\n'%(k,runparams[k]))
    fout.close()
    
    if pcr:
        lpcrparams = ['Compute also the inner primers',
                    'Optimum primer size',
                    'Minimum primer size',
                    'Maximum primer size',
                    'Optimum melting temperature',
                    'Minimum melting temperature',
                    'Maximum melting temperature',
                    'Flanking region size',
                    'Minimum product size',
                    'Maximum product size',
                    'Optimum primer GC content (%)',
                    'Minimum primer GC content (%)',
                    'Maximum primer GC content (%)',
                    'GC clamp',
                    'Bases excluded from the end of the contig']
        pcrparams = {'Compute also the inner primers':inner,
                    'Optimum primer size':optsize,
                    'Minimum primer size':minsize,
                    'Maximum primer size':maxsize,
                    'Optimum melting temperature':opttemp,
                    'Minimum melting temperature':mintemp,
                    'Maximum melting temperature':maxtemp,
                    'Flanking region size':flanksize,
                    'Minimum product size':minprod,
                    'Maximum product size':maxprod,
                    'Optimum primer GC content (%)':optgc,
                    'Minimum primer GC content (%)':mingc,
                    'Maximum primer GC content (%)':maxgc,
                    'GC clamp':gcclamp,
                    'Bases excluded from the end of the contig':exclude}
        fout = open(os.path.join(wdir,'pcrparams.tsv'), 'w')
        for k in lpcrparams:
            fout.write('%s\t%s\n'%(k,pcrparams[k]))
        fout.close()
    #except:pass
    
    if not run_cmd(cmd):
        os.chdir(sdir)
        return False

    # Prepare some data (the long ones to save time on results requests)
    # Convert the pdf files in png
    # Cycle through the directories
    for mapdir in os.listdir(wdir):
        mdir = os.path.join(wdir,mapdir)
        if os.path.isdir(mdir):
            for pdf in os.listdir(mdir):
                if pdf.endswith('.pdf'):
                    # Convert
                    convcmd = 'convert -density 300x300 "%s" -resize x1000 -density 150x150 -trim "%s.png"'                    
                    convcmd = convcmd%(pdf, pdf)
                    proc = subprocess.Popen(convcmd,
                                shell=(sys.platform!="win32"),
                                stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, cwd=mdir)
                    out = proc.communicate()
                    
                    convcmd1 = 'convert -density 300x300 "%s" -resize x400 -density 150x150 -trim "%s_small.png"'
                    convcmd1 = convcmd1%(pdf, pdf)
                    proc = subprocess.Popen(convcmd1,
                                shell=(sys.platform!="win32"),
                                stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, cwd=mdir)
                    out = proc.communicate()
    
    try:
        # Be kind, remove the original files...
        os.remove(contigs)
        for ref in refs:
            os.remove(ref)

        # ...and the contiguator bundle
        os.remove('CONTIGuator.py')
        os.remove('abacas')
    except:pass
    
    # Prepare an archive
    import tarfile
    tarname = os.path.join(wdir,'CONTIGuator_results.tar.gz')
    tar = tarfile.open(tarname,'w:gz')
    for fname in os.listdir(wdir):
        if fname in ['CONTIGuator.py', 'CONTIGuator_results.tar.gz',
                    'abacas', 'newlines',
                    'summary.tsv', 'genparams.tsv', 'runparams.tsv',
                    'pcrparams.tsv']:
            continue
        if os.path.isdir(os.path.join(wdir, fname)):
            idir = os.path.join(wdir, fname)
            for ifname in os.listdir(idir):
                tar.add(os.path.join(idir, ifname),
                        arcname=os.path.join(fname, ifname))
        else:
            tar.add(os.path.join(wdir,fname),
                        arcname=fname)
    tar.close()
    
    # Return back to the original directory
    os.chdir(sdir)

    return True
