FROM python:3.7

# Install bedtools
RUN apt-get update && \
    apt-get install bedtools build-essential default-jre default-jdk pigz

# Install fastqc
RUN cd /opt && \
    wget -q http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod -R 777 /opt/FastQC
ENV PATH="/opt/FastQC:${PATH}"
RUN rm -rf /opt/fastqc_v0.11.9.zip

RUN pip install --upgrade pip && \
    pip install --no-cache-dir biopython scipy pybedtools pandas smart-open shortuuid

RUN ulimit -s unlimited

#COPY src/*.py /src/

#ENTRYPOINT ["/src/methylFASTQ.py"]

# python src/methylFASTQ.py -i inputs/hg38.unmasked.fa -o outputs/ --regions inputs/massie-lvl2-med-str-probes_hg38_chr1.bed  --coverage 250