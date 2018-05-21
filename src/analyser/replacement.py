from collections import namedtuple

import Levenshtein
from itertools import product

import os


def _replace_char(string, char, idx):
    return string[:idx] + char + string[idx + 1:]


def get_possible_seqs(wild_seqs, src_seq, dest_seq):
    results = []
    src_indexes = [idx for idx, ch in enumerate(wild_seqs) if ch == src_seq]
    for src_possibles in product([False, True], repeat=wild_seqs.count(src_seq)):
        possible_seq = wild_seqs
        for idx, possible in zip(src_indexes, src_possibles):
            if not possible:
                continue

            possible_seq = _replace_char(possible_seq, dest_seq, idx)

        if possible_seq != wild_seqs:
            results.append(possible_seq)

    return results


def is_insertion_or_deletion_in_seqs(wild_seqs, target_seqs):
    for data in Levenshtein.opcodes(wild_seqs, target_seqs):
        if data[0] in ('insert', 'delete') and any(data[1:]):
            return True
    return False


def get_reference_seq_file(file_path):
    def remove_linebreak(seqs):
        return seqs[:-1] if seqs.endswith('\n') else seqs

    results = []
    with open(file_path) as f:
        for line in f.readlines():
            try:
                file_name, barcode = line.split(':')
            except ValueError:
                pass
            else:
                reference_info = namedtuple('ReferenceInfo', 'file_name, barcode')
                results.append(reference_info(file_name, remove_linebreak(barcode)))

    return results


# def analyse(file_name, barcode, target_seqs_number, loc_target_seqs, src_seq, dest_seq):
#     with open(os.path.join(os.getcwd(), file_name + '.txt')) as f:
#         for line in f.readlines():
#             if not is_insertion_or_deletion_in_seqs()


def command():
    print('--- replacement 분석기 ---')

    target_seqs_number = input('비율 분석 시퀸스 갯수 : ')

    print('위치 ex. TTTNNNNNATCG barcode가 TTT고 조사 위치가 ATCG라면 5 입력')
    loc_target_seqs_after_barcode = input('비율 분석 시퀸스 위치(barcode 뒤로 몇칸) : ')

    print('변형 조사 (원본 -> 변이)')
    replacement_src_seq = input('변형 조사 원본 시퀸스 : ')
    replacement_dest_seq = input('변형 조사 변이 시퀸스 : ')

    file_name = input('reference seq file 위치 : ')

    for reference_infos in get_reference_seq_file(os.path.join(os.getcwd(), file_name)):
        print(reference_infos)


if __name__ == '__main__':
    command()
