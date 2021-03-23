FROM python:3.8.5-slim-buster

RUN pip install -U \
    pip \
    setuptools \
    wheel

WORKDIR /app
COPY requirements.txt ./
COPY setup.R ./
RUN pip install -r requirements.txt
RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-bioc-deseq2
#RUN Rscript setup.R

ARG GIT_HASH
ENV GIT_HASH=${GIT_HASH:-dev}

#USER user
COPY . .

CMD streamlit run app.py