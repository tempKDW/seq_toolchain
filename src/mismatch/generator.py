from collections import namedtuple
from itertools import repeat, permutations

import os

from src.common import clean_seqs

SEQUENCE = {'A', 'T', 'C', 'G'}
SEQUENCE_STRING = ''.join(SEQUENCE)

BASE_FOLDER_NAME = ''
INPUT_FOLDER_NAME = 'data'
OUTPUT_FOLDER_NAME = 'results'

FIRST_SEQ = 4
TARGET_SEQ_LENGTH = 20
SECOND_SEQ = FIRST_SEQ + TARGET_SEQ_LENGTH
THIRD_SEQ = SECOND_SEQ + 3
FORTH_SEQ = THIRD_SEQ + 3

output_folder = os.path.join(os.getcwd(), BASE_FOLDER_NAME, OUTPUT_FOLDER_NAME)


def _replace_string(string, idx, replacement):
    return string[:idx] + replacement + string[idx + 1:]


def _select_two_index(length):
    results = []
    for i in range(0, length):
        for j in range(i + 1, length):
            results.append((i, j))

    return results


def _generate_sequence_combinations(exclude_first, exclude_second):
    first = SEQUENCE_STRING.replace(exclude_first, '')
    second = SEQUENCE_STRING.replace(exclude_second, '')
    results = set()
    for r, p in zip(repeat(first), permutations(second)):
        for pair in zip(r, p):
            results.add(pair)
    return results


def _put_file_extension(file):
    if not file.endswith('.txt'):
        return file + '.txt'


def generate_one_mismatch(seqs):
    results = []
    key = '1MM-{}position-{}to{}'
    for idx, seq in enumerate(seqs):
        for revised_seq in SEQUENCE:
            if revised_seq == seq:
                continue
            results.append((key.format(idx + 1, seq, revised_seq),
                            _replace_string(seqs, idx, revised_seq)))

    return results


def generate_two_mismatch(seqs):
    if len(seqs) < 2:
        raise Exception('too short sequences!')

    results = []
    key = '2MM-({},{})positions-({}to{},{}to{})'
    for index_pair in _select_two_index(len(seqs)):
        orig_seq_first = seqs[index_pair[0]]
        orig_seq_second = seqs[index_pair[1]]
        for candidate_pair in _generate_sequence_combinations(orig_seq_first, orig_seq_second):
            result = _replace_string(seqs, index_pair[0], candidate_pair[0])
            result = _replace_string(result, index_pair[1], candidate_pair[1])
            results.append((key.format(index_pair[0] + 1, index_pair[1] + 1,
                                       orig_seq_first, candidate_pair[0],
                                       orig_seq_second, candidate_pair[1]),
                            result))
    return results


def add_one_bulge(seqs):
    results = []
    key = 'DNA1bulge-{}position-{}'
    for idx in range(1, len(seqs)):
        for seq in SEQUENCE_STRING:
            results.append((key.format(idx, seq), seqs[:idx] + seq + seqs[idx:]))
    return results


def remove_one_bulge(seqs):
    results = []
    key = 'RNA1bulge-{}position-{}'
    for idx, orig_seq in enumerate(seqs):
        results.append((key.format(idx + 1, orig_seq), seqs[:idx] + seqs[idx + 1:]))

    return results


def get_input_file_and_wild(file_path):
    results = []
    input_data = namedtuple('Input', 'file, wild')
    with open(os.path.join(os.getcwd(), BASE_FOLDER_NAME, INPUT_FOLDER_NAME, file_path), 'r') as f:
        for line in f.readlines():
            try:
                file, wild_seq = line.split(':')
                wild_seq = clean_seqs(wild_seq)
                if not len(wild_seq) == 30:
                    raise BufferError(wild_seq)
            except ValueError:
                print('line : {} is wrong format! it should be FILE_NAME:WILD_SEQ.'.format(line))
                continue
            except BufferError as e:
                print(e.args[0] + ' should be 30 length.')
                continue
            else:
                results.append(input_data(file, wild_seq.upper()))
    return results


def get_partial_seqs(wild):
    return wild[:FIRST_SEQ], wild[FIRST_SEQ:SECOND_SEQ], wild[SECOND_SEQ:THIRD_SEQ], wild[THIRD_SEQ:FORTH_SEQ]


def check_duplicate(wild, output, third_seqs):
    candidates = []
    for seq in SEQUENCE_STRING:
        candidates.append(output + _replace_string(third_seqs, 0, seq))

    for candidate in candidates:
        if candidate in wild:
            return True
    return False


def handle():
    def bulk_write(data):
        for key, seqs in data:
            row_template = '{}:{}:{}'
            if check_duplicate(input_data.wild, seqs, third_seqs):
                row_template += ':wild'
            row_template += '\n'

            output_f.writelines(row_template.format(input_data.file, key, first_seqs + seqs + third_seqs + forth_seqs))

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for input_data in get_input_file_and_wild(input('input file > ')):
        with open(os.path.join(output_folder, _put_file_extension(input_data.file)), 'w') as output_f:
            # 4(first) + 20(19 or 21, second) + 3(third) + 4(forth)
            first_seqs, second_seqs, third_seqs, forth_seqs = get_partial_seqs(input_data.wild)

            bulk_write(generate_one_mismatch(second_seqs))
            bulk_write(generate_two_mismatch(second_seqs))
            bulk_write(add_one_bulge(second_seqs))
            bulk_write(remove_one_bulge(second_seqs))


if __name__ == '__main__':
    handle()
