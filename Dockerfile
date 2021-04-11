FROM ubuntu

RUN mkdir -p /home/app

WORKDIR /home/app

RUN apt update
RUN apt -y install python3-pip

COPY ./requirements.txt ./requirements.txt
RUN pip3 install gunicorn
RUN pip3 install Cython
RUN pip3 install -r requirements.txt


COPY ./data ./data
COPY ./predict ./predict
COPY ./server ./server
COPY ./manage.py ./manage.py

RUN sed -i s/'DEBUG = True'/'DEBUG = False'/g ./server/settings.py
RUN sed -i s/'ALLOWED_HOSTS = \[\]'/'ALLOWED_HOSTS = \["*"]'/g ./server/settings.py

RUN python3 manage.py collectstatic --noinput

EXPOSE 80

CMD ["gunicorn", "--bind", ":80", "--workers", "3", "server.wsgi:application"]