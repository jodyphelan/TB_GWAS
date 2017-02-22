import sys
from collections import defaultdict
from tqdm import tqdm

l_snps = defaultdict(set)
for l in open("lineage_snps.txt"):
	arr = l.rstrip().split()
	l_snps[arr[0]].add(arr[1])

with open(sys.argv[2],"w") as o:
	for l in tqdm(open(sys.argv[1])):
		chrom,pos = l.rstrip().split()[:2]
		if pos not in l_snps[chrom]:
			o.write(l)
