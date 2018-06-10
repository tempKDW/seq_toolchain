import unittest

from src.analyser.replacement import get_possible_seqs, is_insertion_or_deletion_in_seqs


class ReplacementTest(unittest.TestCase):
    def test_get_pairs_should_return_all_possible_seqs(self):
        wild_seqs = 'TACAG'
        src_seq = 'A'
        dest_seq = 'G'

        returned_reqs = get_possible_seqs(wild_seqs, src_seq, dest_seq)

        self.assertListEqual(sorted(returned_reqs), sorted([
            'TGCAG',
            'TACGG',
            'TGCGG',
        ]))

    def test_check_insertion_or_deletion_from_seqs(self):
        wild_seqs = 'TACAG'
        insertion_seqs = 'TACAAG'
        deletion_seqs = 'TACA'
        replacement_seqs = 'TAGAG'

        self.assertTrue(is_insertion_or_deletion_in_seqs(wild_seqs, insertion_seqs))
        self.assertTrue(is_insertion_or_deletion_in_seqs(wild_seqs, deletion_seqs))
        self.assertFalse(is_insertion_or_deletion_in_seqs(wild_seqs, replacement_seqs))
        self.assertFalse(is_insertion_or_deletion_in_seqs(wild_seqs, wild_seqs))


if __name__ == '__main__':
    unittest.main()
