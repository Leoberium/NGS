import sys


for line in sys.stdin:
    chr, left, right, strand = line.rstrip().split()
    print(f'{chr}\t{int(left)+1}\t{int(right)-1}\t{strand}')
    
