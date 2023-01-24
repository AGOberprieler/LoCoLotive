FROM debian:stable-slim

WORKDIR /wd

RUN apt-get update && apt install -y --no-install-recommends bedops bedtools bzip2 ca-certificates curl dirmngr gawk genometools gpg gpg-agent mafft ncbi-blast+ python3 python3-pip r-base r-base-dev wget

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir numpy scipy

RUN Rscript -e 'install.packages("ape", dependencies=TRUE, repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("seqinr", dependencies=TRUE, repos="https://cloud.r-project.org")'

RUN wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 \
    && tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 \
    && mv bin/* /usr/bin/

RUN for server in $(shuf -e ha.pool.sks-keyservers.net \
                            hkp://p80.pool.sks-keyservers.net:80 \
                            keyserver.ubuntu.com \
                            hkp://keyserver.ubuntu.com:80 \
                            pgp.mit.edu) ; do \
        gpg --keyserver "$server" --recv-keys B42F6819007F00F88E364FD4036A9C25BF357DD4 && break || : ; \
    done

RUN curl -o /usr/local/bin/gosu -SL "https://github.com/tianon/gosu/releases/download/1.4/gosu-$(dpkg --print-architecture)" \
    && curl -o /usr/local/bin/gosu.asc -SL "https://github.com/tianon/gosu/releases/download/1.4/gosu-$(dpkg --print-architecture).asc" \
    && gpg --batch --verify /usr/local/bin/gosu.asc /usr/local/bin/gosu \
    && rm /usr/local/bin/gosu.asc \
    && chmod +x /usr/local/bin/gosu

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENV PYTHONUNBUFFERED=1
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
