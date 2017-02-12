# Author: Dongwook Kim / forestkeep21@naver.com

import click
import click_completion
import itertools
import os
import re
import subprocess

from Bio import SeqIO
from datetime import datetime
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError

click_completion.init()
COLLECTION_JOINED = 'joined'


class DuplicatedDBName(Exception):
    pass


class DBConnectionFailed(Exception):
    pass


def connect_to_mongodb():
    try:
        client = MongoClient('localhost', 27017)
        client.server_info()
    except ServerSelectionTimeoutError:
        raise DBConnectionFailed
    else:
        return client


def make_db_name(username):
    return datetime.now().strftime('%Y%m%d%H%M%S') + '-' + username


def init(db_name):
    client = connect_to_mongodb()

    if db_name in client.database_names():
            raise DuplicatedDBName

    db = client[db_name]
    db.create_collection('joined')
    db.get_collection('joined').ensure_index([('seq', 'text')])
    return db


def bulk_insert(db, data, chunk_size):
    while True:
        chunked_data = list(itertools.islice(data, chunk_size))
        if not chunked_data:
            break
        db[COLLECTION_JOINED].insert_many(chunked_data)
    return


def find_all_file_with_path(path, extension='.fastq'):
    file_map = {}
    for (path, _, files) in os.walk(path):
        for file in files:
            if not file.endswith(extension):
                continue

            file_map.setdefault(path, []).append(file)
    return file_map


def is_valid_seq_line(line):
    return True if re.match(r'^(A|T|C|G|a|t|c|g)+$', line) else False


def read_from_fastq_file(path):
    with open(path, 'rU') as handle:
        for line in SeqIO.parse(handle,'fastq'):
            if not is_valid_seq_line(str(line.seq)):
                continue

            yield {
                'id': str(line.id),
                'seq': str(line.seq)
            }


def get_latest_db_name():
    client = connect_to_mongodb()
    db_list = client.database_names()
    for sys_name in ('admin', 'local'):
        db_list.remove(sys_name)

    import pdb; pdb.set_trace()
    print(sorted(db_list))


@click.group()
def mongodb():
    pass


@mongodb.command()
@click.option('--chunk', default=100000, type=int, help='한번에 디비로 넣는 사이즈입니다. 컴퓨터 성능에 따라 조정하세요.')
def insert_joined_data(chunk):
    click.echo('디비명은 현재시간-사용자이름 으로 생성됩니다.')
    click.echo('예를 들어 사용자이름: dongwookkim / 현재시간: 2017-01-02 11:23:34 일 경우,')
    click.echo('20170102112334-dongwookkim 이라는 database 명으로 생성됩니다.')
    username = click.prompt('디비에서 사용할 사용자 이름', type=str)

    db_name = make_db_name(username)

    click.confirm(db_name + ' 으로 생성됩니다, 계속할까요?', abort=True)

    try:
        db = init(db_name)
    except DuplicatedDBName:
        click.echo(db_name + ' 는 이미 있는 이름입니다.')
    except DBConnectionFailed:
        click.echo('db 연결에 실패했습니다. mongoDB가 실행중인지 확인바랍니다.')
        return

    while True:
        path = click.prompt('대상 폴더', type=str)
        file_map = find_all_file_with_path(path)
        for path, files in file_map.items():
            click.echo('[{path}] fastq 파일이 {len} 개 있습니다.'.format(path=path, len=len(files)))
        if click.confirm('모두 맞나요?'):
            break

    files_with_path = []
    for path, files in file_map.items():
        for file in files:
            files_with_path.append('/'.join((path, file)))

    with click.progressbar(files_with_path) as files:
        for file in files:
            chunked_data = read_from_fastq_file(file)
            bulk_insert(db, chunked_data, chunk)


@mongodb.command()
@click.argument('db_name')
@click.argument('input', type=click.File('rb'))
def extract(db_name, input):
    """
    DB_NAME = 추출할 대상이 있는 디비 이름\n
    INPUT = 바코드 파일명
    """
    pass


if __name__ == '__main__':
    mongodb()
