#!/usr/bin/env python

import redis
import time

def add_job(req_id, ip, email, task_id, passphrase=None):
    r = redis.Redis()
    
    jid = 'contiguator_%s'%req_id

    r.zadd('contiguatorjobs', jid, time.time())

    r.hset(jid, 'ip', ip)
    r.hset(jid, 'email', email)
    r.hset(jid, 'task_id', task_id)
    r.hset(jid, 'date', time.asctime())

    if passphrase is not None:
        r.hset(jid, 'passphrase', passphrase)

def retrieve_job(req_id):
    r = redis.Redis()

    return r.hgetall('contiguator_%s'%req_id)
