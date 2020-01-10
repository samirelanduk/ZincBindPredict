import os
try:
    from .secrets import SECRET_KEY
except:
    import binascii
    SECRET_KEY = binascii.hexlify(os.urandom(24)).decode()

VERSION = "0.1.0"

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ALLOWED_HOSTS = []

DEBUG = True

ROOT_URLCONF = "server.urls"

INSTALLED_APPS = [
 "django.contrib.staticfiles",
 "graphene_django",
 "server"
]

DATE_FORMAT = "D j M, Y"
USE_TZ = True
TIME_ZONE = "UTC"

MIDDLEWARE = [
 "django.middleware.common.CommonMiddleware",
]

STATIC_URL = "/static/"
STATIC_ROOT = os.path.abspath(f"{BASE_DIR}/../static")

TEMPLATES = [{
 "BACKEND": "django.template.backends.django.DjangoTemplates",
 "APP_DIRS": True
}]

JOB_EXPIRATION = 10

GRAPHENE = {
 "SCHEMA": "server.schema.schema"
}
