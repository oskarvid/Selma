FROM oskarv/snakemake
MAINTAINER Oskar Vidarsson <oskar.vidarsson@uib.no>

# Install wget and bwa
RUN apt-get update && apt-get install -y \
bwa \
wget \
unzip \
samtools \
libgomp1 \
openjdk-8-jdk && \
rm -rf /var/lib/apt/lists/* && \
rm -rf /usr/share/locale/ /usr/share/man/ /root/.cache

# Download and set up gatk4
RUN mkdir /Jar && \
wget https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip -O /Jar/gatk4.zip && \
unzip -q /Jar/gatk4.zip -d /Jar/ && \
mv /Jar/gatk*/gatk* /Jar/ && \
rm -r /Jar/gatk*/ /Jar/gatk4.zip && \
cp /Jar/gatk /usr/local/bin/gatk

# export path for gatk jar
ENV GATK_LOCAL_JAR=/Jar/gatk-package-4.1.0.0-local.jar

RUN mkdir /tsd /net /work /projects /cluster /references
