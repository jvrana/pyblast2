FROM jvrana/blastdocker:latest


RUN mkdir /code
WORKDIR code
COPY . .
RUN pip install pip --upgrade
RUN pip install .
CMD ["python -m pytest tests"]
