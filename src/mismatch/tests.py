import unittest

from src.mismatch.generator import (generate_one_mismatch, generate_two_mismatch, _select_two_index,
                                    _generate_sequence_combinations, add_one_bulge, remove_one_bulge, get_partial_seqs,
                                    check_duplicate)


def _factorial(a):
    f = 1
    for i in range(1, a+1):
        f = f * i
    return f


def _combination(n, m):
    return _factorial(n) / (_factorial(m) * _factorial(n-m))


class TestCase(unittest.TestCase):
    def test_generate_one_mismatch_length_check(self):
        seqs = 'ATCG'
        self.assertEqual(len(set(generate_one_mismatch(seqs))), len(seqs) * 3)

    def test_generate_one_mismatch_contains(self):
        seqs = 'AA'

        results = generate_one_mismatch(seqs)

        self.assertCountEqual(results, [
            ('1MM-1position-AtoT', 'TA'),
            ('1MM-1position-AtoC', 'CA'),
            ('1MM-1position-AtoG', 'GA'),
            ('1MM-2position-AtoT', 'AT'),
            ('1MM-2position-AtoC', 'AC'),
            ('1MM-2position-AtoG', 'AG'),
        ])

    def test_select_two_index(self):
        self.assertCountEqual(_select_two_index(5),
                              [(0, 1), (0, 2), (0, 3), (0, 4),
                               (1, 2), (1, 3), (1, 4),
                               (2, 3), (2, 4),
                               (3, 4)
                               ])

    def test_generate_sequence_combinations(self):
        self.assertNotIn(('A', 'A'), _generate_sequence_combinations('A', 'C'))
        self.assertNotIn(('A', 'T'), _generate_sequence_combinations('A', 'C'))
        self.assertNotIn(('A', 'C'), _generate_sequence_combinations('A', 'C'))
        self.assertNotIn(('A', 'G'), _generate_sequence_combinations('A', 'C'))
        self.assertNotIn(('A', 'C'), _generate_sequence_combinations('A', 'C'))
        self.assertNotIn(('T', 'C'), _generate_sequence_combinations('A', 'C'))
        self.assertNotIn(('C', 'C'), _generate_sequence_combinations('A', 'C'))
        self.assertNotIn(('G', 'C'), _generate_sequence_combinations('A', 'C'))

    def test_generate_two_mismatch_length(self):
        seqs = 'ATCGATCGATCGATCGATCG'
        self.assertEqual(len(set(generate_two_mismatch(seqs))), _combination(len(seqs), 2) * 9)

    def test_generate_two_mismatch_contains(self):
        seqs = 'AA'
        expected_com = [('2MM-(1,2)positions-(AtoC,AtoC)', 'CC'),
                        ('2MM-(1,2)positions-(AtoC,AtoT)', 'CT'),
                        ('2MM-(1,2)positions-(AtoC,AtoG)', 'CG'),
                        ('2MM-(1,2)positions-(AtoG,AtoC)', 'GC'),
                        ('2MM-(1,2)positions-(AtoG,AtoT)', 'GT'),
                        ('2MM-(1,2)positions-(AtoG,AtoG)', 'GG'),
                        ('2MM-(1,2)positions-(AtoT,AtoC)', 'TC'),
                        ('2MM-(1,2)positions-(AtoT,AtoT)', 'TT'),
                        ('2MM-(1,2)positions-(AtoT,AtoG)', 'TG'),
                        ]

        self.assertCountEqual(generate_two_mismatch(seqs), expected_com)

    def test_add_one_bulge(self):
        self.assertCountEqual(add_one_bulge('AA'),
                              [('DNA1bulge-1position-A', 'AAA'),
                               ('DNA1bulge-1position-C', 'ACA'),
                               ('DNA1bulge-1position-G', 'AGA'),
                               ('DNA1bulge-1position-T', 'ATA'),
                               ])
        self.assertCountEqual(add_one_bulge('ATCG'),
                              [('DNA1bulge-1position-A', 'AATCG'),
                               ('DNA1bulge-1position-T', 'ATTCG'),
                               ('DNA1bulge-1position-C', 'ACTCG'),
                               ('DNA1bulge-1position-G', 'AGTCG'),
                               ('DNA1bulge-2position-A', 'ATACG'),
                               ('DNA1bulge-2position-T', 'ATTCG'),
                               ('DNA1bulge-2position-C', 'ATCCG'),
                               ('DNA1bulge-2position-G', 'ATGCG'),
                               ('DNA1bulge-3position-A', 'ATCAG'),
                               ('DNA1bulge-3position-T', 'ATCTG'),
                               ('DNA1bulge-3position-C', 'ATCCG'),
                               ('DNA1bulge-3position-G', 'ATCGG'),
                               ])

    def test_remove_one_bulge(self):
        self.assertCountEqual(remove_one_bulge('ATCG'),
                              [('RNA1bulge-1position-A', 'TCG'),
                               ('RNA1bulge-2position-T', 'ACG'),
                               ('RNA1bulge-3position-C', 'ATG'),
                               ('RNA1bulge-4position-G', 'ATC'),
                               ])

    def test_get_partial_seqs(self):
        seqs = 'AAAATTTTTTTTTTTTTTTTTTTTCCCGGG'
        self.assertCountEqual(get_partial_seqs(seqs),
                              ['AAAA',
                               'TTTTTTTTTTTTTTTTTTTT',
                               'CCC',
                               'GGG',
                               ])

    def test_check_duplicate(self):
        seqs = 'AAAATTTTTTTTTTTTTTTTTTTTCCCGGG'

        self.assertTrue(check_duplicate(seqs, 'ATTTTTTTTTTTTTTTTTTT', 'CCC'))
        self.assertFalse()