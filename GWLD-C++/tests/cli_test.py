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
        self.assertIn('convert', self.run_cli('--help').stdout)
        for command in ['calculate', 'convert', 'summarize']:
            self.assertIn('--output-prefix', self.run_cli(command, '--help').stdout)
            self.assertIn('--output-prefix', self.run_cli('help', command).stdout)
        for args in [('unknown',), ('calculate', '--unknown'), ('calculate', '--vcf', self.vcf, '--threads', '0'),
                     ('calculate', '--vcf', self.vcf, '--method', 'wrong'),
                     ('calculate', '--vcf', self.vcf, '--decay', '--max-distance-kb', '-1'),
                     ('calculate', '--vcf', self.vcf, '--bed-prefix', 'other'), ('calculate', '--vcf', 'missing.vcf'),
                     ('calculate', '--vcf', self.vcf, '--max-distance-kb', '1'),
                     ('calculate', '--vcf', self.vcf, '--circos-threshold', '0.2'),
                     ('convert', '--vcf', self.vcf, '--format', 'bad'),
                     ('convert', '--vcf', self.vcf, '--method', 'mi'),
                     ('summarize', '--vcf', self.vcf, '--decay'),
                     ('calculate', '--vcf', self.vcf, '--output-prefix', '')]:
            self.run_cli(*args, success=False)
    def test_removed_options(self):
        for command in ['ld', 'freq', 'compute', 'frequency']:
            self.run_cli(command, '--vcf', self.vcf, success=False)
        self.run_cli('--vcf', self.vcf, success=False)
        for option, values in [('--bfile', ['example']), ('--file', ['example-2']),
                               ('--2vcf', []), ('--code012', []), ('--allele-freq', []),
                               ('--by-chrom', []), ('--max-dist', ['1']),
                               ('--threshold', ['0.1']), ('--out-prefix', ['old'])]:
            self.run_cli('calculate', '--vcf', self.vcf, option, *values, success=False)
    def test_short_options(self):
        self.run_cli('calculate', '-v', self.vcf, '-m', 'dprime', '-t', '1', '-o', 'short')
        self.run_cli('calculate', '--vcf', self.vcf, '--method', 'dprime', '--threads', '1', '--output-prefix', 'long')
        self.assertEqual((self.root/'short.heatmap').read_text(), (self.root/'long.heatmap').read_text())
        self.run_cli('convert', '-b', EXAMPLES/'example', '-f', 'dosage', '-o', 'short-bed')
        self.run_cli('convert', '-p', EXAMPLES/'example-2', '-f', 'vcf', '-o', 'short-ped')
        self.assertTrue((self.root/'short-bed.dosage.tsv').exists())
        self.assertTrue((self.root/'short-ped.vcf').exists())
    def test_heatmap_and_streaming(self):
        for method in ['mi', 'rmi', 'r2', 'd', 'dprime']:
            args = ('calculate', '--vcf', self.vcf, '--method', method, '--max-distance-kb', '1', '--decay', '--circos')
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
        for flag, name in [('--bed-prefix', 'example'), ('--ped-prefix', 'example-2')]:
            self.run_cli('convert', flag, EXAMPLES/name, '--format', 'vcf', '-o', 'convert')
            target = self.root/'convert.vcf'
            self.assertTrue(target.read_text().startswith('##fileformat=VCFv4.2'))
            self.run_cli('convert', flag, EXAMPLES/name, '--format', 'dosage', '-o', 'original')
            self.run_cli('convert', '--vcf', target, '--format', 'dosage', '-o', 'roundtrip')
            self.assertEqual((self.root/'original.dosage.tsv').read_text(), (self.root/'roundtrip.dosage.tsv').read_text())
            self.run_cli('calculate', '--vcf', target, '--method', 'mi', '-o', 'reread')
            self.assertTrue((self.root/'reread.heatmap').exists())
            self.run_cli('summarize', flag, EXAMPLES/name)
            self.assertTrue((self.root/'gwld.freq.tsv').exists())
        self.run_cli('summarize', '--vcf', EXAMPLES/'example-3.vcf')
        self.assertFalse((self.root/'gwld.heatmap').exists())
    def test_compressed_and_bcf(self):
        # BCF's coordinate representation is narrower than text VCF's.
        self.vcf.write_text(self.vcf.read_text().replace('5000000000', '500000'))
        self.run_cli('convert', '--vcf', self.vcf, '--format', 'dosage', '-o', 'plain')
        expected = (self.root/'plain.dosage.tsv').read_text()
        self.assertFalse((self.root/'plain.vcf').exists())
        self.assertFalse((self.root/'plain.heatmap').exists())
        compressed = self.root/'ordinary.vcf.gz'
        compressed.write_bytes(gzip.compress(self.vcf.read_bytes()))
        self.run_cli('convert', '--vcf', compressed, '--format', 'dosage', '-o', 'gzip')
        self.assertEqual((self.root/'gzip.dosage.tsv').read_text(), expected)
        if CONVERTER:
            for suffix, mode in [('bgzf.vcf.gz', 'wz'), ('binary.bcf', 'wb')]:
                path = self.root/suffix
                subprocess.run([str(CONVERTER), str(self.vcf), str(path), mode], check=True)
                self.run_cli('convert', '--vcf', path, '--format', 'dosage', '--threads', '2' if OPENMP else '1', '-o', 'converted')
                self.assertEqual((self.root/'converted.dosage.tsv').read_text(), expected)
    def test_genotype_edge_cases(self):
        header = self.vcf.read_text().split('#CHROM')[0] + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n'
        self.vcf.write_text(header + '1\t1\ta\tA\tC\t.\t.\t.\tDP:GT\t2:1/0\t3:0|.\t.\t2:1|1')
        self.run_cli('convert', '--vcf', self.vcf, '--format', 'dosage')
        dosage = (self.root/'gwld.dosage.tsv').read_text().splitlines()
        self.assertEqual(dosage[0].split('\t'), ['CHROM', 'POS', 'ID', 'REF', 'ALT', 's1', 's2', 's3', 's4'])
        self.assertEqual(dosage[-1].split('\t')[5:], ['1', 'NA', 'NA', '2'])
        self.run_cli('summarize', '--vcf', self.vcf, '-o', 'frequency')
        freq = (self.root/'frequency.freq.tsv').read_text().splitlines()
        self.assertEqual(freq[0].split('\t'), ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'REF_FREQ', 'ALT_FREQ', 'MISSING_RATE'])
        self.assertEqual(freq[1].split('\t')[:5], ['1', '1', 'a', 'A', 'C'])
        self.assertEqual(list(map(float, freq[1].split('\t')[5:])), [0.25,0.75,0.5])
        for genotypes in ['0\t0/1\t1/1\t./.', '0/1/1\t0/1\t1/1\t./.',
                          '0/2\t0/1\t1/1\t./.', '0/0\t0/1\t1/1',
                          '0/0\t0/1\t1/1\t./.\t0/0']:
            with self.subTest(genotypes=genotypes):
                self.vcf.write_text(header + '1\t1\ta\tA\tC\t.\t.\t.\tGT\t' + genotypes + '\n')
                self.run_cli('convert', '--vcf', self.vcf, '--format', 'dosage', success=False)
    def test_chromosomes_and_empty(self):
        self.run_cli('calculate', '--vcf', self.vcf, '--by-chromosome', '-o', 'chrom')
        self.assertEqual(len((self.root/'chrom.chr1.heatmap').read_text().splitlines()), 3)
        self.assertEqual(len((self.root/'chrom.chr2.heatmap').read_text().splitlines()), 1)
        self.run_cli('calculate', '--vcf', self.vcf, '--decay', '--max-distance-kb', '0', '-o', 'empty')
        self.assertEqual((self.root/'empty.decay').read_text(), '')
    def test_malformed_inputs(self):
        self.vcf.write_text('##fileformat=VCFv4.2\n')
        self.run_cli('calculate', '--vcf', self.vcf, success=False)
        self.vcf.write_text('#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT s\n1 1 a A C,G . . . GT 0/2\n')
        self.run_cli('calculate', '--vcf', self.vcf, success=False)
        (self.root/'bad.fam').write_text('f s 0 0 1 -9\n')
        (self.root/'bad.bim').write_text('1 a 0 1 A G\n')
        for data in [b'bad', bytes([0x6c, 0x1b, 0x00, 0]), bytes([0x6c, 0x1b, 0x01])]:
            (self.root/'bad.bed').write_bytes(data)
            self.run_cli('calculate', '--bed-prefix', self.root/'bad', success=False)

    def test_sample_identity_roundtrip(self):
        prefix = self.root/'families'
        samples = ['f1 sample 0 0 1 -9', 'f2 sample 0 0 1 -9', 'f:3 sample% 0 0 1 -9']
        prefix.with_suffix('.fam').write_text('\n'.join(samples)+'\n')
        prefix.with_suffix('.bim').write_text('1 a 0 10 G A\n')
        prefix.with_suffix('.bed').write_bytes(bytes([0x6c, 0x1b, 0x01, 35]))
        prefix.with_suffix('.map').write_text('1 a 0 10\n')
        prefix.with_suffix('.ped').write_text('\n'.join(s+' '+g for s,g in zip(samples,['A A','G G','A G']))+'\n')
        for flag in ['--bed-prefix','--ped-prefix']:
            self.run_cli('convert', flag, prefix, '--format', 'vcf', '-o', 'identity')
            header = (self.root/'identity.vcf').read_text().splitlines()[2]
            self.assertEqual(header.split('\t')[9:], ['f1:sample','f2:sample','f%3A3:sample%25'])
            self.run_cli('convert', flag, prefix, '--format', 'dosage', '-o', 'before')
            self.run_cli('convert', '--vcf', self.root/'identity.vcf', '--format', 'dosage', '-o', 'after')
            self.assertEqual((self.root/'before.dosage.tsv').read_text(), (self.root/'after.dosage.tsv').read_text())
        prefix.with_suffix('.fam').write_text(samples[0]+'\n'+samples[0]+'\n'+samples[2]+'\n')
        prefix.with_suffix('.ped').write_text(samples[0]+' A A\n'+samples[0]+' G G\n'+samples[2]+' A G\n')
        for flag in ['--bed-prefix','--ped-prefix']:
            self.run_cli('convert', flag, prefix, '--format', 'vcf', success=False)

    def test_multimodal_likelihood(self):
        table = [[0,2,0],[2,262,5],[2,5,2]]
        pairs = [(i,j) for i in range(3) for j in range(3) for _ in range(table[i][j])]
        header = ('##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n'
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
        lines = [header+'\t'.join(f's{i}' for i in range(len(pairs)))]
        genotypes = ['0/0','0/1','1/1']
        for locus in range(2):
            lines.append(f'1\t{locus+1}\tv{locus}\tA\tG\t.\t.\t.\tGT\t'+'\t'.join(genotypes[p[locus]] for p in pairs))
        self.vcf.write_text('\n'.join(lines)+'\n')
        self.run_cli('calculate','--vcf',self.vcf,'--method','d')
        value = float((self.root/'gwld.heatmap').read_text().splitlines()[0].split()[1])
        self.assertAlmostEqual(value,-0.23361366206657,places=8)

if __name__ == '__main__':
    unittest.main()
