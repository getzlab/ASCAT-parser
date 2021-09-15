FROM ubuntu:16.04
RUN apt-get update && apt-get install -y apt-utils python3.5 python3-pip libssl-dev build-essential gcc make cmake libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev wget bzip2 && pip3 install --upgrade pip==20.3.4
RUN pip3 install pandas scipy numpy gcsfs seaborn pysam matplotlib
RUN apt-get install -y git tar && apt-get autoremove && apt-get autoclean && ln -s /usr/bin/python3 /usr/bin/python

# Copy the built tool and any supporting files into the image
COPY src /src

# Set the working directory
WORKDIR src
