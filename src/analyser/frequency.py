from Levenshtein import editops

from src.common import replace_string

PREFIX_AND_BARCODE_LEN = 19
INDEL_TARGET_LEN = 4


def _get_target_seqs(sorting_seqs, ref_seqs, barcode_with_prefix):
    target_from_sorting = sorting_seqs[sorting_seqs.find(barcode_with_prefix) + len(barcode_with_prefix):]
    target_from_ref = ref_seqs[ref_seqs.find(barcode_with_prefix) + len(barcode_with_prefix):]
    return target_from_sorting, target_from_ref


def is_indel_1(sorting_seqs, ref_seqs, barcode_with_prefix, position):
    target_from_sorting, target_from_ref = _get_target_seqs(sorting_seqs, ref_seqs, barcode_with_prefix)
    start = len(target_from_ref) - position
    end = start + INDEL_TARGET_LEN
    return len(editops(target_from_sorting[start:end],
                       target_from_ref[start:end])) >= 2


def is_indel_2(sorting_seqs, ref_seqs, barcode_with_prefix, position):
    sorting_seqs = put_margin_of_deletion(sorting_seqs, ref_seqs, barcode_with_prefix)
    target_from_sorting, target_from_ref = _get_target_seqs(sorting_seqs, ref_seqs, barcode_with_prefix)
    start = len(target_from_ref) - position
    end = start + INDEL_TARGET_LEN
    return target_from_sorting[start:end][1] == '-' or target_from_sorting[start:end][2] == '-'


def should_be_removed(sorting_seqs, ref_seqs, barcode_with_prefix):
    target_from_sorting, target_from_ref = _get_target_seqs(sorting_seqs, ref_seqs, barcode_with_prefix)

    for data in editops(target_from_ref, target_from_sorting):
        if data[0] in ('insert', 'delete'):
            return True
    return False


def put_margin_of_deletion(sorting_seqs, ref_seqs, barcode_with_prefix):
    target_from_sorting, target_from_ref = _get_target_seqs(sorting_seqs, ref_seqs, barcode_with_prefix)
    result_seqs = sorting_seqs.replace(target_from_sorting, '')

    for data in editops(target_from_ref, target_from_sorting):
        if not data[0] == 'delete':
            continue

        target_from_sorting = replace_string(target_from_sorting, data[1], '-')

    return result_seqs + target_from_sorting


if __name__ == '__main__':
    pass
