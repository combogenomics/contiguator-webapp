#!/usr/bin/env python

import os
import sys
import subprocess
import json
import time
from flask import Flask, request, session, g, redirect, url_for, abort, \
     render_template, flash, escape, Response, send_from_directory, send_file
from werkzeug.utils import secure_filename

from Bio import SeqIO

from utils import generate_hash
from utils import generate_time_hash

from store import add_job
from store import retrieve_job
from store import cumulative_jobs
from store import unique_ips
from store import unique_emails

import settings

app = Flask(__name__)

# App config from settings.py
app.config.from_object(settings)

# Production settings that override the testing ones
try:
    import production
    app.config.from_object(production)
except ImportError:
    pass

# Mail log setup
try:
    import mail_log as ml
    if not app.debug:
        import logging
        mail_handler = ml.TlsSMTPHandler(ml.MAIL_HOST,
                               ml.MAIL_FROM,
                               ml.ADMINS, 'CONTIGuator-webapp Failed!',
                               credentials=(ml.MAIL_USER,
                                            ml.MAIL_PWD))
        mail_handler.setLevel(logging.ERROR)
        app.logger.addHandler(mail_handler)
except ImportError:
    pass

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/run', methods=['GET', 'POST'])
def run():
    # Here handle submissions and run the analysis
    # Send emails on failures, success
    # Use redis to store user stats (hashed for privacy)
    if request.method == 'POST':
        # First things first, compute user hash
        req_id = generate_time_hash(request.remote_addr)

        # To avoid slow-downs in the running directory
        # create subdirs w/ the first 2 chars of the hash
        h2c = req_id[:2]
        try:
            os.mkdir(os.path.join(app.config['UPLOAD_FOLDER'],
                                  h2c))
        except:
            pass

        # Prepare the working directory
        # Our hash scheme ensures that it should be unique
        wdir = os.path.join(app.config['UPLOAD_FOLDER'],
                            h2c, req_id)
        wdir = os.path.abspath(wdir)
        os.mkdir(wdir)

        # Sanity check for form entries
        float_entries = ['maxtemp', 'mingc', 'mintemp', 
                   'maxgc', 'opttemp', 'multitreshold',
                   'evalue', 'optgc']
        int_entries = ['contigcoverage',
                   'minsize', 'threads',
                   'exclude', 'contiglength',
                   'optsize',
                   'flanksize', 'maxsize', 'minprod', 'maxprod',
                   'numN', 'hitlength', 'gcclamp']
        bool_entries = ['pcr', 'blastn', 'non']
        str_entries = ['email', 'jobname']
        for f, entries in ((float, float_entries),
                             (str, str_entries),
                             (bool, bool_entries),
                             (int, int_entries)):
            for e in entries:
                try:
                    if e in request.form:
                        f(request.form[e])
                except Exception as e:
                    flash(u'Something went wrong while processing your options %s'%e,
                          'danger')
                    return redirect(url_for('index'))

        # Save input files
        draft = request.files['contigs']
        if draft:
            filename = secure_filename(draft.filename)
            draft.save(os.path.join(wdir, filename))
            dname = filename
        else:
            flash(u'Forgot to provide your draft genome?',
                  'danger')
            return redirect(url_for('index'))
       
        # Save the genomes files
        genomes = set()
        try:
            for genome in request.files.getlist('reference'):
                filename = secure_filename(genome.filename)
                genome.save(os.path.join(wdir, filename))
                genomes.add(filename)
        except:
            flash(u'Forgot to provide your reference genome?',
                 'danger')
            return redirect(url_for('index'))
        
        # Save the genomes files
        ptts = request.files['pttfile']
        if ptts:
            try:
                for ptt in request.files.getlist('pttfile'):
                    filename = secure_filename(ptt.filename)
                    ptt.save(os.path.join(wdir, filename))
            except:
                flash(u'Something went wrong with your ptt files',
                     'danger')
                return redirect(url_for('index'))
        
        # Check email, hash it
        email = request.form['email']
        if email:
            hemail = generate_hash(email)
        else:
            flash(u'Something went wrong with your email', 'danger')
            return redirect(url_for('index'))

        # Submit the job
        # Then redirect to the waiting page
        try:
            # Handle the bool flags
            if 'non' in request.form:non = True
            else:non = False
            if 'pcr' in request.form:pcr = True
            else:pcr = False
            if 'inner' in request.form:inner = True
            else:inner = False
            if 'blastn' in request.form:blastn = True
            else:blastn = False

            evalue=request.form.get('evalue', 1e-20)
            contiglength=request.form.get('contiglength', 1000)
            contigcoverage=request.form.get('contigcoverage', 20)
            hitlength=request.form.get('hitlength', 1100)
            multitreshold=request.form.get('multitreshold', 1.5)
            non=non
            numN=request.form.get('numN', 100)
            pcr=pcr
            inner=inner
            blastn=blastn
            threads=request.form.get('threads', 1)
            optsize=request.form.get('optsize', 20)
            minsize=request.form.get('minsize', 18)
            maxsize=request.form.get('maxsize', 27)
            opttemp=request.form.get('opttemp', 60)
            mintemp=request.form.get('mintemp', 57)
            maxtemp=request.form.get('maxtemp', 63)
            flanksize=request.form.get('flanksize', 1000)
            minprod=request.form.get('minprod', 1000)
            maxprod=request.form.get('maxprod', 7000)
            optgc=request.form.get('optgc', 50)
            mingc=request.form.get('mingc', 20)
            maxgc=request.form.get('maxgc', 80)
            gcclamp=request.form.get('gcclamp', 1)
            exclude=request.form.get('exclude', 100)
            jobname=request.form.get('jobname', 'CONTIGuator job')

            cmd = 'python tasks.py %s %s %s ' % (req_id,
                                                 wdir,
                                                 dname)
            cmd += '%s %s %s %s ' % (evalue, contiglength, contigcoverage, hitlength)
            cmd += '%s %s %s %s %s %s ' % (multitreshold, non, numN, pcr, inner, blastn)
            cmd += '%s %s %s %s %s %s ' % (threads, optsize, minsize, maxsize, opttemp, mintemp)
            cmd += '%s %s %s %s %s %s %s ' % (maxtemp, flanksize, minprod, maxprod, optgc, mingc, maxgc)
            cmd += '%s %s "%s" ' % (gcclamp, exclude, jobname)
            cmd += ' '.join(genomes)
            f = open(os.path.join(wdir, 'cmd.sh'), 'w')
            f.write(cmd + '\n')
            f.close()
            cmd = 'at -q b -M now -f %s' % os.path.join(wdir, 'cmd.sh')
            proc = subprocess.Popen(cmd,
                                    shell=(sys.platform!="win32"),
				    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
				    stderr=subprocess.PIPE)
	    out = proc.communicate()
	    
	    return_code = proc.returncode
	    if return_code != 0:
		raise Exception('%s'%str(out[1]))
        except Exception as e:
            flash(u'Could not submit your job "%s"' % e,
                  'danger')
            return redirect(url_for('index'))
            
        try:
            # Send details to redis
            add_job(req_id, request.remote_addr, hemail)
        except Exception as e:
            flash(u'Could not save your job details (%s)' % e, 'danger')
            return redirect(url_for('index')) 

        return redirect(url_for('results',
                        req_id=req_id))

    # No POST, return to start
    flash(u'No job details given, would you like to start a new one?',
          'warning')
    return redirect(url_for('index'))

@app.route('/log/<req_id>')
def log(req_id):
    # Get details from redis
    j = retrieve_job(req_id)

    # Return the log, if present
    h2c = req_id[:2]
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id)):
        flash('Could not retrieve the log: is your job older than one week?', 'danger')
        return render_template('index.html')
    if 'CONTIGuator.log' not in os.listdir(os.path.join(
                                      app.config['UPLOAD_FOLDER'],
                                      h2c, req_id)):
        flash('Could not retrieve the log.txt file', 'danger')
        return render_template('error.html', req_id=req_id)

    path = os.path.join(app.config['UPLOAD_FOLDER'],
                        h2c, req_id, 'CONTIGuator.log')
    return Response(''.join(open(path).readlines()),
                    mimetype='text/plain')

@app.route('/err/<req_id>')
def err(req_id):
    # Get details from redis
    j = retrieve_job(req_id)

    # Return the log, if present
    h2c = req_id[:2]
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id)):
        flash('Could not retrieve the log: is your job older than one week?', 'danger')
        return render_template('index.html')
    if 'log.err' not in os.listdir(os.path.join(
                                      app.config['UPLOAD_FOLDER'],
                                      h2c, req_id)):
        flash('Could not retrieve the log.err file', 'danger')
        return render_template('error.html', req_id=req_id)

    path = os.path.join(app.config['UPLOAD_FOLDER'],
                        h2c, req_id, 'log.err')
    return Response(''.join(open(path).readlines()),
                    mimetype='text/plain')

@app.route('/archive/<req_id>')
def archive(req_id):
    # Get details from redis
    j = retrieve_job(req_id)

    # Return the archive, if present
    h2c = req_id[:2]
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id)):
        flash('Could not retrieve the archive: is your job older than one week?', 'danger')
        return render_template('index.html')

    if 'CONTIGuator_results.tar.gz' not in os.listdir(os.path.join(
                                      app.config['UPLOAD_FOLDER'],
                                      h2c, req_id)):
        flash('Could not retrieve the archive', 'danger')
        return render_template('error.html', req_id=req_id)

    path = os.path.join(app.config['UPLOAD_FOLDER'],
                        h2c, req_id)
    return send_from_directory(path,
                               'CONTIGuator_results.tar.gz',
                               as_attachment=True)

@app.route('/unmapped/<req_id>')
def unmapped(req_id):
    # Get details from redis
    j = retrieve_job(req_id)

    # Return the file, if present
    h2c = req_id[:2]
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id)):
        flash('Could not retrieve the file: is your job older than one week?', 'danger')
        return render_template('index.html')
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id, 'UnMappedContigs')):
        flash('No contig has been mapped to yor genome (or error occurred)); you'+
               'may want to relaunch a new job with different parameters', 'danger')
        return render_template('index.html')
    if 'UnMappedContigs.txt' not in os.listdir(os.path.join(
                                      app.config['UPLOAD_FOLDER'],
                                      h2c, req_id, 'UnMappedContigs')):
        flash('Could not retrieve the UnMappedContigs.txt file',
                'danger')
        return render_template('error.html', req_id=req_id)

    path = os.path.join(app.config['UPLOAD_FOLDER'],
                        h2c, req_id, 'UnMappedContigs',
                        'UnMappedContigs.txt')
    return Response(''.join(open(path).readlines()),
                    mimetype='text/plain')

@app.route('/mapped/<req_id>/<mdir>')
def mapped(req_id, mdir):
    # Get details from redis
    j = retrieve_job(req_id)

    # Return the log, if present
    h2c = req_id[:2]
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id)):
        flash('Could not retrieve the file: is your job older than one week?', 'danger')
        return render_template('index.html')
    if 'MappedContigs.txt' not in os.listdir(os.path.join(
                                      app.config['UPLOAD_FOLDER'],
                                      h2c, req_id, mdir)):
        flash('Could not retrieve the MappedContigs.txt file',
                'danger')
        return render_template('error.html', req_id=req_id)

    path = os.path.join(app.config['UPLOAD_FOLDER'],
                        h2c, req_id, mdir, 'MappedContigs.txt')
    return Response(''.join(open(path).readlines()),
                    mimetype='text/plain')

@app.route('/pdf/<req_id>/<mdir>/<fname>')
def pdf(req_id, mdir, fname):
    # Get details from redis
    j = retrieve_job(req_id)

    # Return the log, if present
    h2c = req_id[:2]
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id)):
        flash('Could not retrieve the document: is your job older than one week?', 'danger')
        return render_template('index.html')
    if fname not in os.listdir(os.path.join(
                                      app.config['UPLOAD_FOLDER'],
                                      h2c, req_id, mdir)):
        flash('Could not retrieve the %s file'%fname,
                'danger')
        return render_template('error.html', req_id=req_id)

    path = os.path.join(app.config['UPLOAD_FOLDER'],
                        h2c, req_id, mdir, fname)
    return send_file(path, mimetype='application/pdf') 

@app.route('/png/<req_id>/<mdir>/<fname>')
def png(req_id, mdir, fname):
    # Get details from redis
    j = retrieve_job(req_id)

    # Return the log, if present
    h2c = req_id[:2]
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id)):
        flash('Could not retrieve the picture: is your job older than one week?', 'danger')
        return render_template('index.html')
    if fname not in os.listdir(os.path.join(
                                      app.config['UPLOAD_FOLDER'],
                                      h2c, req_id, mdir)):
        flash('Could not retrieve the %s file'%fname,
                'danger')
        return render_template('error.html', req_id=req_id)

    path = os.path.join(app.config['UPLOAD_FOLDER'],
                        h2c, req_id, mdir, fname)
    return send_file(path, mimetype='image/png') 

@app.route('/scaffold/<req_id>/<mdir>')
def scaffold(req_id, mdir):
    # Get details from redis
    j = retrieve_job(req_id)

    # Return the log, if present
    h2c = req_id[:2]
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id)):
        flash('Could not retrieve the scaffold: is your job older than one week?', 'danger')
        return render_template('index.html')
    if 'PseudoContig.fsa' not in os.listdir(os.path.join(
                                      app.config['UPLOAD_FOLDER'],
                                      h2c, req_id, mdir)):
        flash('Could not retrieve the scaffold file',
                'danger')
        return render_template('error.html', req_id=req_id)

    path = os.path.join(app.config['UPLOAD_FOLDER'],
                        h2c, req_id, mdir, 'PseudoContig.fsa')
    return Response(''.join(open(path).readlines()),
                    mimetype='text/plain')

@app.route('/pcr/<req_id>/<mdir>')
def pcr(req_id, mdir):
    # Get details from redis
    j = retrieve_job(req_id)

    # Return the log, if present
    h2c = req_id[:2]
    if not os.path.exists(os.path.join(
                              app.config['UPLOAD_FOLDER'],
                              h2c, req_id)):
        flash('Could not retrieve the file: is your job older than one week?', 'danger')
        return render_template('index.html')
    if 'PCRPrimers.tsv' not in os.listdir(os.path.join(
                                      app.config['UPLOAD_FOLDER'],
                                      h2c, req_id, mdir)):
        flash('Could not retrieve the PCR summary file',
                'danger')
        return render_template('error.html', req_id=req_id)

    path = os.path.join(app.config['UPLOAD_FOLDER'],
                        h2c, req_id, mdir, 'PCRPrimers.tsv')
    return Response(''.join(open(path).readlines()),
                    mimetype='text/plain')

@app.route('/legend')
def legend():
    return render_template('legend.html')

@app.route('/results/<req_id>')
def results(req_id):
    # Here show the results or the wait page
    # Get the right job using the session or the hash key
    
    # Get details from redis
    j = retrieve_job(req_id)

    h2c = req_id[:2]
    status = j['status']
    if status == 'Job done':
        # run results logics
        try:
            result = json.load(open(os.path.join(app.config['UPLOAD_FOLDER'],
                                                 h2c, req_id, 'result.json')))
        except Exception as e:
            app.logger.error('Internal server error: %s\nRequest ID: %s' % (e, req_id))
            flash(u'Internal server error: %s'%e, 'danger')
            return render_template('error.html', req_id=req_id)

        # Handle the results and prepare the data for page rendering
        wdir = os.path.join(app.config['UPLOAD_FOLDER'],
                            req_id[:2], req_id)
        wdir = os.path.abspath(wdir)
        
        # Summary data
        genparams = []
        if 'genparams.tsv' in os.listdir(wdir):
            for l in open(os.path.join(wdir, 'genparams.tsv')):
                s = l.strip().split('\t')
                if len(s) == 1:
                    genparams.append((s[0], ''))
                else:
                    genparams.append((s[0], s[1]))
        # Run parameters
        runparams = []
        if 'runparams.tsv' in os.listdir(wdir):
            for l in open(os.path.join(wdir, 'runparams.tsv')):
                s = l.strip().split('\t')
                if len(s) == 1:
                    runparams.append((s[0], ''))
                else:
                    runparams.append((s[0], s[1]))
        # PCR parameters
        pcrparams = []
        if 'pcrparams.tsv' in os.listdir(wdir):
            for l in open(os.path.join(wdir, 'pcrparams.tsv')):
                s = l.strip().split('\t')
                if len(s) == 1:
                    pcrparams.append((s[0], ''))
                else:
                    pcrparams.append((s[0], s[1]))
        
        # Summary stats
        summary = []
        if 'summary.tsv' in os.listdir(wdir):
            b = True
            for l in open(os.path.join(wdir, 'summary.tsv')):
                if b:
                    b = False
                    continue
                summary.append(l.strip().split('\t'))
        # No summary stats = errors? (or unmapped)     
        if len(summary) == 0: 
            return render_template('nomapped.html', req_id=req_id)
 
        # Are there any Map_ directories there?
        maps = filter(lambda x: os.path.isdir(
                                    os.path.join(wdir, x)) and
                                x.startswith('Map_'),
                                os.listdir(wdir))
        # Collect data for each mapped reference replicon
        mapdata = []
        for m in maps:
            d = {}
            d['dir'] = m
            # Reference length
            d['reflen'] = len(SeqIO.read(open(os.path.join(wdir,
                                m, 'Reference.embl')),'embl'))
            # Mapped contigs
            mcontigs = [int(x.strip().split('\t')[1]) for x in open(os.path.join(wdir,
                                    m, 'MappedContigs.txt'))]
            d['mapped'] = len(mcontigs)
            d['mappedbp'] = sum(mcontigs)
            # PDF/PNG map
            d['png'] = list(filter(lambda x: x.endswith('.png') and
                                    not x.endswith('_small.png'),
                            os.listdir(os.path.join(wdir, m))))[0]
            d['pdf'] = d['png'].rstrip('.png')
            # PCR?
            if 'PCRPrimers.tsv' in os.listdir(os.path.join(wdir, m)):
                d['pcr'] = len(filter(lambda x: len(x) >= 5 and
                                'Left Contig' not in x[0],
                                [x.rstrip().split('\t')
                                    for x in open(os.path.join(wdir, m,
                                        'PCRPrimers.tsv'))]))
            # Wrap it all together, keep the same order
            mapdata.append(d)

        return render_template('result.html', req_id=req_id,
                                            mapdata=mapdata,
                                            genparams=genparams,
                                            runparams=runparams,
                                            pcrparams=pcrparams,
                                            summary=summary,
                                            maps=maps)
    elif status == 'Job failed':
        error_msg = j.get('error', '')
        app.logger.error('Internal server error: %s\nRequest ID: %s' % (error_msg,
                         req_id))
        flash(u'Internal server error: %s' % error_msg,
               'danger')
        return render_template('error.html', req_id=req_id)
    else:
        # If too much time has passed, it means that the job has either failed
        # or anything like that
        cur_time = time.time()
        start_time = float(j['time'])
        delta_time = cur_time - start_time
        if (60 * 15) < delta_time < (60 * 30):
            flash(u'Your job exceeded 15 minutes, something might have gone wrong!' +
                  u'Will try 15 more minutes before giving up',
                  'danger')
        elif delta_time > (60 * 30):
            flash(u'Your job exceeded 30 minutes, something must have gone wrong!' +
                  u'If your genomes are many (and big) you might want to run Medusa locally',
                  'danger')
            return render_template('error.html', req_id=req_id)
        return render_template('waiting.html', status=status)

@app.route('/stats')
def stats():
    return render_template('stats.html')

@app.route('/stats/jobs')
def jobs():
    return Response(json.dumps([{'date':x[1],
                                 'jobs':x[0]} for x in cumulative_jobs()]),
                    mimetype='text/plain')

@app.route('/stats/ips')
def ips():
    return Response(json.dumps([{'date':x[1],
                                 'ips':x[0]} for x in unique_ips()]),
                    mimetype='text/plain')

@app.route('/stats/emails')
def emails():
    return Response(json.dumps([{'date':x[1],
                                 'emails':x[0]} for x in unique_emails()]),
                    mimetype='text/plain')

@app.route('/admin')
def admin():
    # Here admin section: upload a new medusa
    # Clean manually the jobs
    flash('Not implemented yet', 'warning')
    return render_template('index.html')

if __name__ == '__main__':
    app.run()
