#!/usr/bin/env python

DEBUG=True
SECRET_KEY='development key'
CELERY_BROKER_URL='redis://localhost:6379'
CELERY_RESULT_BACKEND='redis://localhost:6379'
UPLOAD_FOLDER = 'uploads'
MAX_CONTENT_LENGTH = 1024 * 1024 * 1024
