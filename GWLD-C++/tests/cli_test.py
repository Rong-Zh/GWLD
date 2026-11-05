import math
import gzip
import pathlib
import subprocess
import sys
import tempfile
import unittest

EXE, EXAMPLES = map(pathlib.Path, sys.argv[1:3])
OPENMP = len(sys.argv) < 4 or sys.argv[3] == 'ON'
CONVERTER = pathlib.Path(sys.argv[4]) if len(sys.argv) > 4 else None
sys.argv = sys.argv[:1]

class CLI(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.temp.name)
        self.vcf = self.root / 'input.vcf'
        self.vcf.write_text('##fileformat=VCFv4.2\n'
            '##contig=<ID=1,length=6000000000>\n##contig=<ID=2,length=1000000>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n'
            '1\t2000\ta\tA\tG\t.\t.\t.\tDP:GT\t8:0/0\t9:0/1\t9:1|1\t4:./.\n'
            '2\t100\tb\tC\tT\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t0/0\n'
            '1\t1000\tc\tG\tT\t.\t.\t.\tGT\t0/0\t1/0\t1/1\t./.\n'
            '1\t5000000000\td\tA\tC\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/0\n')
    def tearDown(self):
        self.temp.cleanup()
    def run_cli(self, *args, success=True):
        result = subprocess.run([str(EXE), *map(str, args)], capture_output=True, text=True, cwd=self.root)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result
    def test_help_version_errors(self):
        self.assertIn('2.0.0', self.run_cli('--version').stdout)
        self.run_cli('--help')
        for args in [('--unknown',), ('--vcf', self.vcf, '--threads', '0'),
                     ('--vcf', self.vcf, '--method', 'wrong'), ('--vcf', self.vcf, '--max-dist', '-1'),
                     ('--vcf', self.vcf, '--bfile', 'other'), ('--vcf', 'missing.vcf')]:
            self.run_cli(*args, success=False)
    def test_heatmap_and_streaming(self):
        for method in ['MI', 'RMI', 'r2', 'r^2', 'D', "D'"]:
            args = ('--vcf', self.vcf, '--method', method, '--max-dist', '1', '--decay', '--circos')
            self.run_cli(*args, '--threads', '1', '-o', 'one')
            self.assertFalse((self.root / 'one.heatmap').exists())
            self.run_cli(*args, '--threads', '2' if OPENMP else '1', '--heatmap', '-o', 'two')
            for suffix in ['decay', 'circos']:
                self.assertEqual((self.root / f'one.{suffix}').read_text(), (self.root / f'two.{suffix}').read_text())
            row = (self.root / 'one.decay').read_text().strip().split('\t')
            self.assertEqual(row[:3], ['a', 'c', '1000'])
            matrix = [[float(v) for v in line.split()] for line in (self.root/'two.heatmap').read_text().splitlines()]
            self.assertEqual(len(matrix), 4)
            self.assertTrue(math.isnan(matrix[0][0]))
            self.assertAlmostEqual(matrix[0][2], float(row[3]))
    def test_formats_and_roundtrip(self):
        for flag, name in [('--bfile', 'example'), ('--file', 'example-2')]:
            self.run_cli(flag, EXAMPLES/name, '--2vcf', '--allele-freq', '-o', 'convert')
            target = self.root/'convert.vcf'
            self.assertTrue(target.read_text().startswith('##fileformat=VCFv4.2'))
            self.run_cli('--vcf', target, '--method', 'MI', '-o', 'reread')
            self.assertTrue((self.root/'reread.heatmap').exists())
        self.run_cli('--vcf', EXAMPLES/'example-3.vcf', '--allele-freq')
    def test_compressed_and_bcf(self):
        # BCF's coordinate representation is narrower than text VCF's.
        self.vcf.write_text(self.vcf.read_text().replace('5000000000', '500000'))
        self.run_cli('--vcf', self.vcf, '--code012', '-o', 'plain')
        expected = (self.root/'plain.vcf').read_text()
        compressed = self.root/'ordinary.vcf.gz'
        compressed.write_bytes(gzip.compress(self.vcf.read_bytes()))
        self.run_cli('--vcf', compressed, '--code012', '-o', 'gzip')
        self.assertEqual((self.root/'gzip.vcf').read_text(), expected)
        if CONVERTER:
            for suffix, mode in [('bgzf.vcf.gz', 'wz'), ('binary.bcf', 'wb')]:
                path = self.root/suffix
                subprocess.run([str(CONVERTER), str(self.vcf), str(path), mode], check=True)
                self.run_cli('--vcf', path, '--code012', '--threads', '2' if OPENMP else '1', '-o', 'converted')
                self.assertEqual((self.root/'converted.vcf').read_text(), expected)
    def test_genotype_edge_cases(self):
        header = self.vcf.read_text().split('#CHROM')[0] + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n'
        self.vcf.write_text(header + '1\t1\ta\tA\tC\t.\t.\t.\tDP:GT\t2:1/0\t3:0|.\t.\t2:1|1')
        self.run_cli('--vcf', self.vcf, '--code012')
        self.assertEqual((self.root/'gwld.vcf').read_text().splitlines()[-1].split('\t')[9:], ['1', '-1', '-1', '2'])
        self.run_cli('--vcf', self.vcf, '--allele-freq', '-o', 'frequency')
        self.assertEqual(list(map(float,(self.root/'frequency.freq').read_text().split())), [0.25,0.75,0.5])
        for genotypes in ['0\t0/1\t1/1\t./.', '0/1/1\t0/1\t1/1\t./.',
                          '0/2\t0/1\t1/1\t./.', '0/0\t0/1\t1/1',
                          '0/0\t0/1\t1/1\t./.\t0/0']:
            with self.subTest(genotypes=genotypes):
                self.vcf.write_text(header + '1\t1\ta\tA\tC\t.\t.\t.\tGT\t' + genotypes + '\n')
                self.run_cli('--vcf', self.vcf, '--code012', success=False)
    def test_chromosomes_and_empty(self):
        self.run_cli('--vcf', self.vcf, '--by-chrom', '-o', 'chrom')
        self.assertEqual(len((self.root/'chrom.chr1.heatmap').read_text().splitlines()), 3)
        self.assertEqual(len((self.root/'chrom.chr2.heatmap').read_text().splitlines()), 1)
        self.run_cli('--vcf', self.vcf, '--decay', '--max-dist', '0', '-o', 'empty')
        self.assertEqual((self.root/'empty.decay').read_text(), '')
    def test_malformed_inputs(self):
        self.vcf.write_text('##fileformat=VCFv4.2\n')
        self.run_cli('--vcf', self.vcf, success=False)
        self.vcf.write_text('#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT s\n1 1 a A C,G . . . GT 0/2\n')
        self.run_cli('--vcf', self.vcf, success=False)
        (self.root/'bad.fam').write_text('f s 0 0 1 -9\n')
        (self.root/'bad.bim').write_text('1 a 0 1 A G\n')
        for data in [b'bad', bytes([0x6c, 0x1b, 0x00, 0]), bytes([0x6c, 0x1b, 0x01])]:
            (self.root/'bad.bed').write_bytes(data)
            self.run_cli('--bfile', self.root/'bad', success=False)

if __name__ == '__main__':
    unittest.main()
