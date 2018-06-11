import re


def clean_seqs(string):
    m = re.match(r'[aAtTcCgG]+', string)
    return m.group()


def replace_string(string, idx, replacement):
    return string[:idx] + replacement + string[idx + 1:]