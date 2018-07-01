import os
import re
from collections import namedtuple, defaultdict
from pathlib import Path

from Bio.pairwise2 import align
from Bio.SubsMat.MatrixInfo import blosum62


ref_datum = namedtuple('RefData', ['file', 'barcode', 'full', 'position'])
ref_data = []


def clean_seqs(string):
    m = re.match(r'[aAtTcCgG]+', string)
    return m.group()


def check_insertion(target_ref_seqs, except_position):
    index = 0
    for seq in target_ref_seqs.strip('-'):
        if except_position[0] <= index < except_position[1]:
            continue

        if seq != '-':
            index += 1
        else:
            return True
    return False


def check_deletion(sort_seqs, except_position, barcode):
    seqs = sort_seqs[sort_seqs.find(barcode):]
    index = 0
    for seq in seqs:
        if except_position[0] <= index < except_position[1]:
            continue

        if seq != '-':
            index += 1
        else:
            return True
    return False


def check_indel(sort_seqs, position, barcode):
    seqs = sort_seqs[sort_seqs.find(barcode):]
    index = 0
    for seq in seqs:
        if position[0] <= index <= position[1] and seq == '-':
            return True

        if seq != '-':
            index += 1
    return False


while True:
    ref_file = input('reference 파일 > ')
    if not ref_file.endswith('.txt'):
        ref_file += '.txt'
    if os.path.exists(ref_file):
        break
    print(ref_file + ' 이런 파일 없음')

while True:
    sort_base_folder = input('sort 폴더 (엔터 치면 reference 동일 위치) > ')
    if not sort_base_folder:
        sort_base_folder = os.path.dirname(ref_file)
        break

    sort_base_folder = sort_base_folder
    if os.path.exists(sort_base_folder):
        break
    print(sort_base_folder + ' 이런 폴더 없음')

while True:
    result_base_folder = input('result 폴더 > ')
    if os.path.exists(result_base_folder):
        break
    print(result_base_folder + ' 이런 폴더 없음')

while True:
    try:
        indel_start = int(input('InDel 포지션(뒤부터) > '))
        indel_len = int(input('InDel 몇개나? > '))
    except TypeError:
        print('오타 ㄴㄴ')
    else:
        break

with open(ref_file, 'r') as f:
    for line in f.readlines():
        sorting_file, barcode, full = line.split(':')
        sorting_file += '.txt'
        barcode, full = clean_seqs(barcode), clean_seqs(full)
        start = len(full) - indel_start
        end = start + indel_len
        ref_data.append(ref_datum(sorting_file, barcode, full, (start, end)))

counts = defaultdict(dict)
for datum in ref_data:
    indels = []
    normals = []
    removed = []
    total = 0

    try:
        with open(os.path.join(sort_base_folder, datum.file), 'r') as f:
            for line in f.readlines():
                seqs = clean_seqs(line)
                total += 1
                pref, psort, _, _, _ = align.localds(datum.full, seqs, blosum62, -10, -1)[0]
                # pref, psort, _, _, _ = align.localms(datum.full, seqs, 5, -4, -2, -0.5)[0]
                if (check_insertion(pref, datum.position)
                        or check_deletion(psort, datum.position, datum.barcode)):
                    removed.append(seqs)
                elif check_indel(psort, datum.position, datum.barcode):
                    indels.append(seqs)
                else:
                    start = seqs.find(datum.barcode)
                    end = start + (len(datum.full) - len(datum.barcode))
                    normals.append(seqs[start:end])

    except FileNotFoundError:
        print(datum.file + ' < 이런 이름의 파일 없음')
        continue
    else:
        count = defaultdict(dict)
        for seqs in normals:
            for idx, seq in enumerate(seqs):
                if seq in count[idx]:
                    count[idx][seq] += 1
                else:
                    count[idx][seq] = 1

        counts[datum.file]['count'] = count
        counts[datum.file]['ref'] = datum.full[datum.full.find(datum.barcode) + len(datum.barcode):]

    counts[datum.file]['indels'] = str(len(indels))
    counts[datum.file]['indel_seqs'] = indels
    counts[datum.file]['removed'] = str(len(removed))
    counts[datum.file]['normals'] = str(len(normals))
    counts[datum.file]['total'] = str(total)

result_folder = os.path.join(result_base_folder, 'result')
indel_folder = os.path.join(result_base_folder, 'result_indel')
for folder in [result_folder, indel_folder]:
    if not os.path.exists(folder):
        os.mkdir(folder)

for file, data in counts.items():
    with open(os.path.join(result_folder, file.split('.')[0] + '_count.txt'), 'w') as f:
        buffer = []
        buffer.append('-,' + ','.join(data['ref']))
        buffer.append('A,' + ','.join([str(counts.get('A', 0)) for idx, counts in data['count'].items()]))
        buffer.append('C,' + ','.join([str(counts.get('C', 0)) for idx, counts in data['count'].items()]))
        buffer.append('G,' + ','.join([str(counts.get('G', 0)) for idx, counts in data['count'].items()]))
        buffer.append('T,' + ','.join([str(counts.get('T', 0)) for idx, counts in data['count'].items()]))
        buffer = [line + '\n' for line in buffer]

        f.writelines(buffer)

    with open(os.path.join(result_folder, file.split('.')[0] + '_info.txt'), 'w') as f:
        f.writelines([
            'indel count: ' + data['indels'] + '\n',
            'remove count: ' + data['removed'] + '\n',
            'normal count: ' + data['normals'] + '\n',
            'total count: ' + data['total'] + '\n',
        ])

    with open(os.path.join(indel_folder, file.split('.')[0] + '.txt'), 'w') as f:
        data['indel_seqs'] = [datum + '\n' for datum in data['indel_seqs']]
        f.writelines(data['indel_seqs'])
