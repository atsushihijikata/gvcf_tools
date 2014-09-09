"""
gvcf_tools.py

@version: 0.1

@description: Converting gVCF format into several formats.
 Currently two formats are supported.
 1. multi-FASTA format with variations in the gVCF file.
 2. conventional VCF format.

@author:  Atsushi Hijikata 
          atsushi.hijikata@gmail.com
@created: July 21, 2014

Usage:
    gvcf = gVCF(gvcf_file)
    gvcf.parse()
    gvcf.write_fasta()

Prerequisite:
    1. Python 2.6 or later,
    2. genome data in FASTA format (please download from public databases,
       e.g. NCBI, EBI or UCSC).

"""
import os, sys

VERSION='0.1'
# change as your environment
GENOME='/path/to/genome_fasta/'
PREFIX='chr' # ex. 'hs_ref_GRCh37.p13_'
SUFFIX='.fa' # ex. '.fa'

# parameters for analysis
MIN_DP   = 10 # lower limit for base call
MIN_QUAL = 50 # lower limit for variant call

def find_max(value_list):
    i = 0
    max_v = 0
    max_i = 0
    for i in range(len(value_list)):
        v = int(value_list[i])
        if v > max_v:
            max_i = i
            max_v = v
    return max_i

def ambiguous_code(ref, alt):
    # use if a location contains heterozygous alleles 
    amb_code = {'K': ('G','T'), # Keto
                'M': ('A','C'), # aMino
                'R': ('A','G'), # puRine
                'Y': ('C','T'), # pYrimidine
                'S': ('C','G'), # Strong
                'W': ('A','T')  # Weak
                }

    for c, t in amb_code.items():
        c1, c2 = t
        if ref == c1 and alt == c2 or \
           ref == c2 and alt == c1:
               return c

def get_genome_sequence(chrom=None):
    genome_fasta = GENOME + PREFIX + chrom + SUFFIX
    file = open(genome_fasta,'r')
    gseq = ''
    for line in file:
        if line.startswith('>'): continue
        line = line.rstrip()
        gseq += line
    file.close()
    return gseq

def get_genome_subsequence(gseq=None, start=None, end=None):
    """
    The letters of nucleotides in the reference are shown in upper case.
    """
    subseq = gseq[start-1:end]
    return subseq.upper()

class VCF:
    def __init__(self, data):
        # data is expected each variant line in VCF
        cols = data.split()
        self.chrom    = cols[0]
        self.start    = int(cols[1])
        self.end      = int(cols[1])
        self.ref      = cols[3]
        self.alt      = cols[4]
        self.qual     = cols[5]
        self.filt     = cols[6]
        self.info     = cols[7]
        self.frmt     = cols[8]
        self.vals     = cols[9]
        self.is_var   = False # Variant flag
        self.has_two  = False # Two or more alternative alleles
        self.vtype    = None  # Variation type
        self.AD       = '-'   # Depth per allele
        self.GT       = '-'   # Genotype
        self.DP       = 0     # sequence depth
        self.MQ       = 0     # mapping quality 
        self.GQX      = 0     # genotype quality
        self.alts     = []    # cases for two or more variants in a site
        self.is_valid = True
        if self.info.startswith('END='):
            vs = self.info.split(';')
            en = vs[0].replace('END=','')
            self.end = int(en)

        #if self.alt != '.':
        self.alt = self.alt.replace(',<NON_REF>','')
        if self.alt != '<NON_REF>': # line for variant(s)
            self.is_var = True

        try:
            head = self.frmt.split(':')
            vals = self.vals.split(':')
            for i in range(len(head)):
                setattr(self, head[i], vals[i])
        except:
            pass

    def assign(self):
        if not self.is_var:
            return

        if self.alt.find(',') != -1: # two or more non-ref alleles
            self.alts = self.alt.split(',')
            self.has_two = True
        else: # a single non-ref allele
            self.alts.append(self.alt)

        if len(self.alts) > 2: # unrealistic case, three or more variants exist
            sys.stderr.write("Warning: %s has three or more variant.\n" % data)
            self.is_valid = False

        if self.has_two:
            if self.GT == '1/1': # seems to be homozygous while two non-ref
                n = find_max(self.AD.split(','))
                self.alt = self.alts[n-1] # 0 is for reference allele
                self.has_two = False
   
        # filter by quality
        if float(self.qual) < MIN_QUAL:
            self.filt = 'LowQual'
        else:
            self.filt = 'PASS'

        # assign variant type (SNV/MNV/DIV)
        if self.has_two:
            na = self.alts[0]
            nb = self.alts[1]
        else:
            na = self.ref
            nb = self.alt

        if len(na) == len(nb):
            if len(na) == 1:
                self.vtype = 'SNV'
            else:
                self.vtype = 'MNV'
        else:
            self.vtype = 'DIV'

class gVCF:
    def __init__(self, gvcf_file):
        self.gvcf_file = gvcf_file
        self.vcf_list = []
        self.seqregions = []
        self.homo_variants = 0
        self.hetero_variants = 0

        dummy = self.gvcf_file.split('/')[-1]
        self.title = dummy.split('.')[0]

    def parse(self):
        """
        parse the gVCF file and convert each line into the VCF class object.
        """
        file = open(self.gvcf_file,'r')
        n = 0
        sys.stderr.write("Parsing gVCF file: %s..." % self.gvcf_file)
        for line in file:
            #if line.startswith('chr'):
            if not line.startswith('#'):
                self.vcf_list.append(VCF(line.rstrip()))
                n += 1
        sys.stderr.write("\nDone.\n%d lines were read.\n" % n)
        self.seqregion()
        file.close()

    def seqregion(self):
        """
        compile the sequenced regions separated by line 
        """
        flag = 0
        for sr in self.vcf_list:
            sr.assign()
            if flag == 0:
                st = sr.start
                en = sr.end
                ch = sr.chrom
                flag = 1
                continue
            elif sr.start - en == 1 and sr.chrom == ch:
                en = sr.end

            elif sr.start - en > 1 or sr.chrom != ch:
                self.seqregions.append((ch,st,en))
                ch = sr.chrom
                st = sr.start
                en = sr.end
        self.seqregions.append((ch,st,en))

    def add_variation(self, ch, st, en, sq):
        vsq = ''
        n_v = 0
        n_n = 0
        for v in self.vcf_list:
            if ch != v.chrom or st > v.start or en < v.start:
                continue
            if v.is_var:
                sys.stderr.write("%s %s %s DP:%s, qual:%s, %s " % \
                        (v.chrom, v.start, v.end, v.DP, v.qual, v.vtype))
                if float(v.qual) < MIN_QUAL:
                    vsq += 'N' * len(v.ref)
                    sys.stderr.write('LOW QUALITY, REJECTED.\n')
                    n_n += 1
                    continue
                else:
                    sys.stderr.write("PASS\n")
                n_v += 1

                if v.vtype == 'SNV' or v.vtype == 'MNV':
                    if v.GT == '1/1':   # homozygous
                        vsq += v.alt
                        self.homo_variants += 1
                    elif v.GT == '0/1': # heterozygous
                        ambcode = ''
                        if v.has_two:
                            for ai in range(len(v.alts[0])):
                                ambcode += ambiguous_code(v.alts[0][ai],v.alts[1][ai])
                        else:
                            for ai in range(len(v.alts[0])):
                                ambcode += ambiguous_code(v.ref[ai],v.alt[ai])
                        vsq += ambcode
                        self.hetero_variants += 1

                else: # DIV
                    if v.GT == '1/1': # homozygous
                        vsq += v.alt
                        self.homo_variants += 1
                    else: # heterozygous 
                        vsq += v.alt.lower()
                        self.hetero_variants += 1
            else:
                sys.stderr.write("%s %s %s DP:%s\n" % \
                        (v.chrom, v.start, v.end, v.DP))
                if int(v.DP) >= MIN_DP: 
                    vsq += sq[v.start-st:v.end-st+1]
                else:
                    vsq += 'N' * len(sq[v.start-st:v.end-st+1])
        return vsq, n_v, n_n

    def to_fasta(self,reference=False):
        """
        convert the sequence region into fasta format sequence
        with variation
        chromosome data is required for this conversion.
        """
        oner = 70 # residues per line
        fasta = ''

        sys.stderr.write("%d regions are converting to FASTA.\n" % len(self.seqregions))
        gseq = None
        pchr = None
        detected = rejected = 0
        for r in self.seqregions:
            chrom, start, end = r
            if gseq == None or chrom != pchr:
                gseq = get_genome_sequence(chrom=chrom)
            pchr = chrom
            sequence = get_genome_subsequence(gseq=gseq, start=start, end=end)
            if not reference:
                sequence, n_v, n_n = self.add_variation(chrom, start, end, sequence)

            # output in fasta
            fasta += '>%s_%s_%d_%d [%d variants, %d rejected]\n' % (self.title, chrom,\
                    start, end, n_v, n_n)
            detected += n_v
            rejected += n_n
            for i in range(len(sequence)/oner+1):
                if len(sequence) <= (i+1)*oner:
                    fasta += sequence[i*oner:] + '\n'
                else:
                    fasta += sequence[i*oner:(i+1)*oner] + '\n'

        sys.stderr.write("%d variants were detected.\n" % detected)
        sys.stderr.write("%d variants were rejected.\n" % rejected)

        return fasta

    def write_fasta(self):
        out_file = self.title + ".fa"
        fasta = self.to_fasta()
        out = open(out_file,'w')
        out.write(fasta)
        out.close()

    def to_vcf(self):
        vcf_out = ''
        for v in self.vcf_list:
            if v.is_var:
                vcf_out += '%s\t%d\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                        (v.chrom,v.start,v.ref,v.alt,v.qual,v.filt,v.info,v.frmt,v.vals)
        return vcf_out

    def write_vcf(self):
        out_file = self.title + '.vcf'
        vcf_out = ''
        for v in self.vcf_list:
            if v.is_var:
                vcf_out += '%s\t%d\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                        (v.chrom,v.start,v.ref,v.alt,v.qual,v.filt,v.info,v.frmt,v.vals)
        out = open(out_file,'w')
        out.write(vcf_out)
        out.close()


if __name__ == '__main__':

    usage = "usage: %s gVCF" % sys.argv[0]

    """ for Python 2.7.x or later
    import argparse
    p = argparse.ArgumentParser(description="gVCF tools")
    p.add_argument('gvcf_file')
    p.add_argument('-g','--genome',default='./',help="Genome Fasta Directory")
    p.add_argument('-p','--prefix',default='chr',help="Prefix of chromosome fasta files")
    p.add_argument('-s','--suffix',default='.fa',help="Suffix of chromosome fasta files")
    p.add_argument('--min_depth',default=10,
            help="The minimum depth to accept base calls in a region.")
    p.add_argument('--min_qual', default=50,
            help="The minimum quality value to accept a variant call.")
    p.add_argument('--version',action='version',version='%(prog)s '+VERSION)
    args = p.parse_args()
    GENOME = args.genome
    PREFIX = args.prefix
    SUFFIX = args.suffix
    MIN_DP   = args.min_depth
    MIN_QUAL = args.min_qual

    gvcf = gVCF(args.gvcf_file)
    """

    from optparse import OptionParser
    p = OptionParser(usage=usage, version='%prog ' + VERSION)
    p.add_option('-g','--genome',dest='GENOME',default=GENOME,type='string',
                 
                 help="Genome fasta directory (default: %s)" % GENOME)
    p.add_option('-p','--prefix',dest='PREFIX',default=PREFIX,type='string',
            help="Prefix of chromosome fasta (defalut: %s)" % PREFIX)
    p.add_option('-s','--suffix',dest='SUFFIX',default=SUFFIX,type='string',
            help="Suffix of chromosome fasta (default: %s)" % SUFFIX)
    p.add_option('--min_depth',dest='MIN_DP',default=MIN_DP,type='int',
            help="The minimum depth to accept base calls in a region. (default: %d)" % MIN_DP)
    p.add_option('--min_qual',dest='MIN_QUAL',default=MIN_QUAL,type='int',
            help="The minimum quality value to accept a variant call. (default: %d)" % MIN_QUAL)

    options, args = p.parse_args()

    if len(args) != 1:
        print usage
        sys.exit()

    GENOME   = options.GENOME
    PREFIX   = options.PREFIX
    SUFFIX   = options.SUFFIX
    MIN_DP   = options.MIN_DP
    MIN_QUAL = options.MIN_QUAL

    # read and parse gVCF file
    gvcf = gVCF(args[0])
    gvcf.parse()

    # convert to fasta
    #print gvcf.to_fasta()
    gvcf.write_fasta()

    # convert to conventional VCF
    #print gvcf.to_vcf()

    # write VCF to a file
    gvcf.write_vcf()
