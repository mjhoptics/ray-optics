import unittest
from rayoptics.codev import reader


class TestReadSeqBuffer(unittest.TestCase):
    def test_continuation_proc(self):
        inpt = ['ab&', 'cd&', 'ef', 'gh']
        outpt = reader.read_seq_buffer(inpt)
        self.assertEqual(outpt, [['abcdef'], ['gh']])

    def test_eol_comment_precedence(self):
        inpt = ['ab; cd!ef; gh']
        outpt = reader.read_seq_buffer(inpt)
        self.assertEqual(outpt, [['ab'], ['cd']])

    def test_end_of_cmd(self):
        inpt = ['ab; cd; ef', 'gh;ij& ', 'kl; mn']
        outpt = reader.read_seq_buffer(inpt)
        self.assertEqual(outpt, [['ab'], ['cd'], ['ef'], ['gh'], ['ijkl'],
                                 ['mn']])

    def test_spaces_as_delimiters(self):
        inpt = ['s &', '.07 10 ', 's & ', ' 0 .001']
        outpt = reader.read_seq_buffer(inpt)
        self.assertEqual(outpt, [['s', '.07', '10'], ['s', '0', '.001']])

    def test_strip_comment_lines(self):
        inpt = ['tit "this is a title"', ' ! this is a comment line', 'dim m',
                'so 0 1e11 ! infinite object']
        outpt = reader.read_seq_buffer(inpt)
        self.assertEqual(outpt, [['tit', 'this is a title'], ['dim', 'm'],
                                 ['so', '0', '1e11']])

    def test_empty_buffer(self):
        inpt = []
        outpt = reader.read_seq_buffer(inpt)
        self.assertEqual(outpt, [])

    def test_empty_cmd_strings(self):
        inpt = ['', '! this comment will be stripped out', '', '    ! so will this one']
        outpt = reader.read_seq_buffer(inpt)
        self.assertEqual(outpt, [])


class TestNextLine(unittest.TestCase):
    def test_next_line_with_cont(self):
        inpt = ['ab&', 'cd&', 'ef']
        it = iter(inpt)
        outpt = reader.next_line(it)
        self.assertEqual(outpt, 'abcdef')

    def test_next_line(self):
        inpt = ['ab&', 'cd', 'ef']
        it = iter(inpt)
        outpt = reader.next_line(it)
        self.assertEqual(outpt, 'abcd')
        outpt = reader.next_line(it)
        self.assertEqual(outpt, 'ef')
        with self.assertRaises(StopIteration):
            outpt = reader.next_line(it)


class TestStripComments(unittest.TestCase):
    def test_strip_comment_line(self):
        inpt = '!abcd'
        outpt = reader.strip_comments(inpt)
        self.assertEqual(outpt, '')
        inpt = ' !abcd'
        outpt = reader.strip_comments(inpt)
        self.assertEqual(outpt, ' ')

    def test_strip_eol_comment(self):
        inpt = 'ab!cd'
        outpt = reader.strip_comments(inpt)
        self.assertEqual(outpt, 'ab')

    def test_strip_comment_no_comment(self):
        inpt = 'abcd'
        outpt = reader.strip_comments(inpt)
        self.assertEqual(outpt, 'abcd')


if __name__ == '__main__':
    unittest.main(verbosity=2)
