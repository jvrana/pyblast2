FROM biopython/biopython:latest

RUN mkdir /home/code
ADD . /home/code
WORKDIR /home/code

RUN chmod 666 *

RUN pip3 install pip --upgrade
RUN pip install . --ignore-installed