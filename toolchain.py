# Author: Dongwook Kim / forestkeep21@naver.com

import asyncio
import click
import click_completion
import itertools
import os
import re
import motor.motor_asyncio

from Bio import SeqIO
from datetime import datetime
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError

click_completion.init()

ENV_MONGODB_KEY = 'CURRENT_MONGODB'

MAX_THREAD_SIZE = 20

loop = asyncio.get_event_loop()


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


def connect_to_mongodb_with_motor():
    try:
        client = motor.motor_asyncio.AsyncIOMotorClient('localhost', 27017)
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
    return db


def bulk_insert(db, data, chunk_size, created_at):
    async def async_bulk_insert(coll, data):
        await coll.insert_many(data)

    joined_collname = '{}-{}'.format(created_at.strftime('%Y%m%d%H%M%S'), 'joined')


    tasks = []
    while True:
        chunked_data = list(itertools.islice(data, chunk_size))
        if not chunked_data:
            break
        tasks.append(async_bulk_insert(db[joined_collname], chunked_data))
        if len(tasks) >= MAX_THREAD_SIZE:
            loop.run_until_complete(asyncio.wait(tasks))
            tasks = []
    return


def find_all_file_with_path(path, extensions=('.fastq', '.txt')):
    file_map = {}
    for (path, _, files) in os.walk(path):
        for file in files:
            if not any([file.endswith(extension) for extension in extensions]):
                continue

            file_map.setdefault(path, []).append(file)
    return file_map


def is_valid_seq_line(line):
    return True if re.match(r'^(A|T|C|G|a|t|c|g)+$', line) else False


def read_from_fastq_file(path):
    with open(path, 'rU') as handle:
        for line in SeqIO.parse(handle, 'fastq'):
            if not is_valid_seq_line(str(line.seq)):
                continue

            yield {
                'seq': str(line.seq)
            }


def read_from_txt_file(path):
    with open(path, 'r') as f:
        for line in f:
            if not is_valid_seq_line(str(line)):
                continue

            yield {
                'seq': str(line)
            }


def check_db_is_exist(client, dbname):
    return dbname in client.database_names()


@click.group()
def mongodb():
    pass


@mongodb.command()
@click.option('--chunk', default=10000, type=int, help='한번에 디비로 넣는 사이즈입니다. 컴퓨터 성능에 따라 조정하세요.')
def insert_joined_data(chunk):
    click.echo('data 를 mongodb에 집어넣는 작업을 시작합니다.')

    while True:
        db_name = click.prompt('디비명', type=str)
        click.confirm(db_name + ' 으로 생성됩니다, 계속할까요?', abort=True)

        try:
            db = init(db_name)
        except DuplicatedDBName:
            click.echo(db_name + ' 는 이미 있는 이름입니다.')
        except DBConnectionFailed:
            click.echo('db 연결에 실패했습니다. mongoDB가 실행중인지 확인바랍니다.')
            click.get_current_context().abort()
        else:
            break

    os.environ[ENV_MONGODB_KEY] = db_name

    while True:
        path = click.prompt('대상 폴더', type=str)
        file_map = find_all_file_with_path(path)
        for path, files in file_map.items():
            click.echo('[{path}] 파일이 {len} 개 있습니다.'.format(path=path, len=len(files)))
        if click.confirm('모두 맞나요?'):
            break

    files_with_path = []
    for path, files in file_map.items():
        for file in files:
            files_with_path.append('/'.join((path, file)))

    now = datetime.now()

    joined_collname = '{}-{}'.format(now.strftime('%Y%m%d%H%M%S'), 'joined')
    db[joined_collname].ensure_index([('seq', 'text')])

    db = connect_to_mongodb_with_motor()[db_name]
    with click.progressbar(files_with_path) as files:
        for file in files:
            if file.endswith('.fastq'):
                chunked_data = read_from_fastq_file(file)
            elif file.endswith('.txt'):
                chunked_data = read_from_txt_file(file)
            else:
                click.echo('잘못된 형식의 파일입니다.')
                click.get_current_context().abort()
                return

            bulk_insert(db, chunked_data, chunk, now)


def parse_barcode_file(file):
    barcode_maps = {}
    with open(file, 'r') as f:
        for line in f.readlines():
            elements = line.split(':')
            try:
                barcode_maps[elements[0].strip()] = elements[1].strip()
            except IndexError:
                click.echo('{} : {}'.format(file, line))
                click.echo('잘못된 형식의 바코드입니다.')
                click.get_current_context().abort()
    return barcode_maps


async def select_mongodb_by_barcode(source_coll, dest_coll, key, barcode):
    barcode = barcode.upper()
    # extracted_data = list(source_coll.find({'seq': {'$regex': barcode}}))
    data = []
    async for row in source_coll.find({'seq': {'$regex': barcode}}):
        data.append(row)
    dest_coll.insert_one({
        '_id': key,
        'barcode': barcode,
        'extracted_data': data
    })


def get_db_name_by_index(db_maps):
    while True:
        for idx, dbname in db_maps.items():
            click.echo('[{}] {}'.format(idx, dbname))
        click.echo('[] 안의 번호를 입력해주세요. [2] dwkim 을 선택한다면 2')
        db_idx = click.prompt('작업을 수행할 db 번호', type=int)
        try:
            dbname = db_maps[db_idx]
        except IndexError:
            click.echo('잘못입력했습니다. 종료하려면 Ctrl+D')
        else:
            return dbname


def get_coll_name_by_index(coll_maps):
    while True:
        for idx, collname in coll_maps.items():
            click.echo('[{}] {}'.format(idx, collname))
        click.echo('[] 안의 번호를 입력해주세요. [2] dwkim 을 선택한다면 2')
        coll_idx = click.prompt('작업을 수행할 collection 번호', type=int)
        try:
            collname = coll_maps[coll_idx]
        except IndexError:
            click.echo('잘못입력했습니다. 종료하려면 Ctrl+D')
        else:
            return collname


def get_doc_id_by_index(doc_list):
    while True:
        for idx, doc in enumerate(doc_list):
            click.echo('[{}] {}'.format(idx, doc))
        click.echo('[] 안의 번호를 입력해주세요. [2] dwkim 을 선택한다면 2')
        doc_idx = click.prompt('작업을 수행할 doc 번호', type=int)
        try:
            doc_id = doc_list[doc_idx]
        except IndexError:
            click.echo('잘못입력했습니다. 종료하려면 Ctrl+D')
        else:
            return doc_id


def get_barcode_file():
    while True:
        file = click.prompt('바코드 파일 위치', type=str)
        if os.path.exists(file):
            return file
        click.echo('해당 위치에 파일이 없습니다. 종료하려면 Ctrl+D')


def upsert_result(coll):
    for doc in coll.find({'_id': {'$ne': 'result_info'}}):
        coll.update({'_id': 'result_info'},
                    {'$set': {doc['_id']: len(doc['extracted_data'])}},
                    upsert=True)


@mongodb.command()
def extract():
    client = connect_to_mongodb()

    click.echo('mongodb에 있는 data를 파일에 있는 각 barcode 별로 나눕니다.')

    dbname = None
    if ENV_MONGODB_KEY in os.environ:
        dbname = os.environ[ENV_MONGODB_KEY]
        if not click.confirm('현재 작업하던 db인 {} 에 그대로 하실건가요?'.format(dbname)):
            dbname = None

    if not dbname:
        db_list = client.database_names()
        db_list.remove('admin')
        db_list.remove('local')
        db_maps = {idx + 1: dbname for idx, dbname in enumerate(db_list)}
        dbname = get_db_name_by_index(db_maps)

    revised_coll_list = [coll for coll in client[dbname].collection_names() if coll.endswith('joined')]
    coll_maps = {idx + 1: coll_name for idx, coll_name in enumerate(revised_coll_list)}
    collname = get_coll_name_by_index(coll_maps)
    file = get_barcode_file()

    barcode_maps = parse_barcode_file(file)
    now = datetime.now()
    client = connect_to_mongodb_with_motor()
    source_coll = client[dbname][collname]
    dest_coll = client[dbname]['{}-{}'.format(now.strftime('%Y%m%d%H%M%S'), 'extracted')]
    with click.progressbar(barcode_maps) as keys:
        tasks = []
        for key in keys:
            barcode = barcode_maps[key]
            tasks.append(select_mongodb_by_barcode(source_coll, dest_coll, key, barcode))
            if len(tasks) >= MAX_THREAD_SIZE:
                loop.run_until_complete(asyncio.wait(tasks))
                tasks = []

    client = connect_to_mongodb()
    dest_coll = client[dbname]['{}-{}'.format(now.strftime('%Y%m%d%H%M%S'), 'extracted')]
    upsert_result(dest_coll)


@mongodb.command()
def extract_to_file():
    client = connect_to_mongodb()

    click.echo('mongodb에 있는 extract된 데이터를 파일로 뽑습니다.')

    db_list = client.database_names()
    db_list.remove('admin')
    db_list.remove('local')
    db_maps = {idx + 1: dbname for idx, dbname in enumerate(db_list)}
    dbname = get_db_name_by_index(db_maps)

    revised_coll_list = [coll for coll in client[dbname].collection_names() if coll.endswith('extracted')]
    coll_maps = {idx + 1: coll_name for idx, coll_name in enumerate(revised_coll_list)}
    collname = get_coll_name_by_index(coll_maps)

    doc_list = [data['_id'] for data in client[dbname][collname].find() if data['_id'] != 'result_info']
    doc_id = get_doc_id_by_index(doc_list)

    with open('{}.txt'.format(doc_id), 'w') as f:
        data = client[dbname][collname].find_one({'_id': doc_id})['extracted_data']
        for row in data:
            f.write('{}\n'.format(row['seq']))

if __name__ == '__main__':
    mongodb()
    loop.close()
