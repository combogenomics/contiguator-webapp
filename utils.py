#!/usr/bin/env python

def generate_time_hash(data):
    '''Return a sha256 hash
    
    Salted hash including the client data
    and date/time
    '''
    import hashlib
    from datetime import datetime
    from time import strftime
    
    return hashlib.sha256("salegrosso" +
             data +
             datetime.now().strftime("%Y%m%d%H%M%S%f")).hexdigest()

def generate_hash(data):
    '''Return a sha256 hash
    
    Salted hash including only the client data
    '''
    import hashlib
    
    return hashlib.sha256("salegrosso" +
             data).hexdigest()
