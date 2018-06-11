from unittest import TestCase

from src.analyser.frequency import is_indel_1, is_indel_2, put_margin_of_deletion


class FrequencyTestCase(TestCase):
    prefix = 'TTTG'
    barcode = 'GTGCACACACATATA'

    valid_ref_seqs = (prefix
                      + barcode
                      + 'ACTGGAACACAAAGCATAGACTGCGGGGCG')
    _pre_buffer_sorting = 'CTTGAAAAAGTGGCACCGAGTCGGTGCTTT'
    valid_indel1_sorting_seqs = (_pre_buffer_sorting
                                 + prefix
                                 + barcode
                                 + 'ACTGGAACACAAAGCATAGCGGGGCGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA')
    valid_indel2_sorting_seqs = (_pre_buffer_sorting
                                 + prefix
                                 + barcode
                                 + 'ACTGGAACACAAAGCATAGCGGGGCG')
    invalid_indel1_sorting_seqs = (_pre_buffer_sorting
                                   + prefix
                                   + barcode
                                   + 'ACTGGAACACAAAGCATAGCGGGGCGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA')

    def test_is_indel_1(self):
        self.assertTrue(is_indel_1(self.valid_indel1_sorting_seqs, self.valid_ref_seqs, self.prefix + self.barcode, 11))

    def test_is_indel_2(self):
        self.assertTrue(is_indel_2(self.valid_indel2_sorting_seqs, self.valid_ref_seqs, self.prefix + self.barcode, 11))
