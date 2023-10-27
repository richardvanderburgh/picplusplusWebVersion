FROM picbuild AS source-stage

FROM ubuntu:20.04

COPY --from=source-stage /PIC++Main /PIC++Main
COPY start.sh /
COPY djangodocker /djangodocker/
COPY templates/index.html /templates/

RUN apt-get update && \
apt-get install -y python3 python3-pip 

RUN pip3 install Django \
&& pip3 install django-cors-headers \
&& pip3 install gunicorn==20.1.0 \
&& pip install whitenoise \
&& pip install tzdata

EXPOSE 8000

# Start Gunicorn to serve the Django app
CMD ["/start.sh"]
