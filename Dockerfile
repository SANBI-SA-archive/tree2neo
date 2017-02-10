FROM python:2.7
MAINTAINER Thoba Lose "thoba@sanbi.ac.za"

RUN apt-get update -y --fix-missing \
    && apt-get upgrade -y \
    && mkdir -p /code/data \
    && pip install -U pip

RUN wget "https://drive.google.com/uc?export=download&id=0B5cdXx4b_kIyaDZGWEtBZF85TEU" \
    -O vcf.tar.bz2
RUN tar xvfj vcf.tar.bz2

#RUN wget "https://zenodo.org/record/252101/files/data.tar.bz2" \
#    -O data.tar.bz2
#RUN tar xvfj data.tar.bz2

COPY requirements.txt /code
RUN pip install -r /code/requirements.txt

COPY . /code
WORKDIR /code
RUN pip install --editable .
CMD ["vcf2neo" ,"init", "-D", "/vcf"]
