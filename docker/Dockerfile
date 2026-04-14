FROM python:3.11-slim-bullseye

LABEL maintainer="Nadeem Khan <nkhan119@uottawa.ca>"
LABEL description="Multivariate GWAS Statistical Analysis — CDC 1.0.0"
LABEL version="1.0.0"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        libopenblas-dev \
        liblapack-dev \
        pkg-config \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        git \
        curl \
        && rm -rf /var/lib/apt/lists/*

# Core_stack
RUN pip install --upgrade pip setuptools wheel

RUN pip install \
    numpy==1.26.4 \
    scipy==1.13.1 \
    pandas==2.2.2 \
    statsmodels==0.14.2 \
    scikit-learn==1.5.0

# GWAS_specific
RUN pip install \
    pyarrow==16.1.0

# Causal_inference
RUN pip install \
    dowhy==0.11.1 \
    econml==0.15.1 \
    networkx==3.3

# Visualization
RUN pip install \
    matplotlib==3.9.0 \
    plotly==5.22.0 \
    kaleido==0.2.1

# Excel_table
RUN pip install \
    openpyxl==3.1.3 \
    xlsxwriter==3.2.0

# Stats_extras
RUN pip install \
    pingouin==0.5.4 \
    lifelines==0.29.0 \
    linearmodels==6.0 \
    jinja2==3.1.4 \
    tqdm==4.66.4


WORKDIR /workspace
CMD ["/bin/bash"]
