import argparse
import time

parser = argparse.ArgumentParser(description='This script filter out rows from the syncronized file with more than bi-allelic counts.')
parser.add_argument('--input', dest='input', type=str, required=True, help='Input file containing allele counts in a synchronized file.')

'''
input format:

3L      8724429 G       63:0:0:39:0:0   279:0:0:179:0:0 147:0:0:85:0:0  84:0:0:107:0:0  ...
3L      8724539 A       21:0:137:0:0:0  59:0:337:0:0:0  38:0:193:1:0:0  65:0:158:0:0:0  ...
3L      8736299 A       129:0:0:124:0:0 225:0:0:221:0:0 94:0:0:110:0:0  98:1:0:0:0:0     ...

Note: The 2nd and 3rd row will be filtered out.

'''

parser.add_argument('--output', dest='output', type=str, required=True, help='Output file containing the filtered data in the same format.')
args = parser.parse_args()

infile=open(args.input,'r')
outfile=open(args.output,'w')

start = time.time()

for l in infile:	
    a=l.split()
    lists=[]
    for item in a[3:]:
        lists.append(sum([map(int,x) for x in [item.split(":")]],[]))
    if(sum([not v==0 for v in map(sum, zip(*lists))]) < 3):
        print >> outfile, "\t".join(a)

print 'Time elapsed:', time.time()-start

infile.close()
outfile.close()

