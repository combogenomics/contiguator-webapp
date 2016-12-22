#!/usr/bin/env python

import redis
import time

def add_job(req_id, ip, email, passphrase=None):
    r = redis.Redis()
    
    jid = 'contiguator_%s'%req_id

    r.zadd('contiguatorjobs', jid, time.time())

    r.hset(jid, 'ip', ip)
    r.hset(jid, 'email', email)
    r.hset(jid, 'date', time.asctime())
    r.hset(jid, 'time', time.time())
    r.hset(jid, 'status', 'Job not started')

    if passphrase is not None:
        r.hset(jid, 'passphrase', passphrase)

def update_job(req_id, key, value):
    r = redis.Redis()
    
    jid = 'contiguator_%s'%req_id

    r.hset(jid, key, value)

def retrieve_job(req_id):
    r = redis.Redis()

    return r.hgetall('contiguator_%s'%req_id)

def slice_it(li, cols=10):
    start = 0
    for i in xrange(cols):
        stop = start + len(li[i::cols])
        yield li[start:stop]
        start = stop

def cumulative_jobs(size=100):
    r = redis.Redis()

    i = 0

    res = []

    for j in r.zrange('contiguatorjobs', 0, -1):
        i += 1
        date = r.hget(j, 'date')
        res.append( (i, date) )

    for s in slice_it(res, cols=size):
        if len(s) > 0:
            yield s[-1]

def unique_emails(size=100):
    r = redis.Redis()

    ip = set()

    res = []

    for j in r.zrange('contiguatorjobs', 0, -1):
        i = r.hget(j, 'email')
        date = r.hget(j, 'date')
        if i not in ip:
            ip.add(i)
            res.append( (len(ip), date) )

    for s in slice_it(res, cols=size):
        if len(s) > 0:
            yield s[-1] 

def unique_ips(size=100):
    r = redis.Redis()

    ip = set()

    res = []

    for j in r.zrange('contiguatorjobs', 0, -1):
        i = r.hget(j, 'ip')
        date = r.hget(j, 'date')
        if i not in ip:
            ip.add(i)
            res.append( (len(ip), date) )

    for s in slice_it(res, cols=size):
        if len(s) > 0:
            yield s[-1]
