FROM python:3.8

RUN pip install pycryptodome
COPY checker.py /root/
CMD ["/usr/local/bin/python3", "-u", "/root/checker.py"]
