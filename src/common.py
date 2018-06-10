import re


def clean_seqs(string):
    m = re.match(r'[aAtTcCgG]+', string)
    return m.group()
