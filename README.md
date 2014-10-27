contiguator-webapp
=============

Scaffold the contigs!

Test the server
---------------

Please make sure that redis is up and running.
Also, put in the "static" directory a 3.x version of bootstrap and jquery.min.js in the static/js directory.
Create a directory called contiguator-app, with the CONTIGuator.py file (web branch) and the abacas script.

    virtualenv venv
    source venv/bin/activate
    pip install -r requirements.txt
    python contiguator.py

In another terminal:
    
    source venv/bin/activate
    celery -A contiguator.celery worker

That's it! Open a browser on the same machine and got to 127.0.0.1:5000

Production
----------

First, install all the dependencies

    sudo pip install -r requirements.txt

Then install in apache the contiguator.conf file (changing the paths).
Move the celeryd.conf file in /etc/default/celeryd, and install the celeryd script in /etc/init.d, using your OS command to install the service.

Create a production.py file which can then be used to override the settings.py debug options.

Restart apache and start celery and redis.

To update the server once the upstream repository has been updated, just run git pull and the restart apache and celery.

Edit the mail_log.py file to setup the error logging through email.

You may also want to set up a cron job to wipe out the uploads directly every now and then
