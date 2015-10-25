#!/usr/bin/env python

from celery import Celery

def make_celery(app):
    celery = Celery(app.import_name, broker=app.config['CELERY_BROKER_URL'],
                    backend=app.config['CELERY_RESULT_BACKEND'])
    celery.conf.update(app.config)
    TaskBase = celery.Task
    class ContextTask(TaskBase):
        abstract = True
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)
    celery.Task = ContextTask
    return celery

def is_task_failed(task, task_id):
    res = task.AsyncResult(task_id)
    return res.failed()

def is_task_ready(task, task_id):
    res = task.AsyncResult(task_id)
    return res.ready()
