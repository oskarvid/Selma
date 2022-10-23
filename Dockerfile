FROM samtools-gnuplot-bwa-snakemake-bcftools
MAINTAINER Oskar Vidarsson <oskarvidarsson@gmail.com>

# Download and set up gatk4
RUN mkdir /Jar && \
wget --no-check-certificate https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip -O /Jar/gatk4.zip && \
unzip -q /Jar/gatk4.zip -d /Jar/ && \
mv /Jar/gatk*/gatk* /Jar/ && \
rm -r /Jar/gatk*/ /Jar/gatk4.zip

ENV PATH="$PATH:/Jar/"

# export path for gatk jar
ENV GATK_LOCAL_JAR=/Jar/gatk-package-4.3.0.0-local.jar
