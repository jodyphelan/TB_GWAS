import sys
from collections import defaultdict


def matrix_filter(pos_dict,matrix_file):
	import subprocess
	tabix = "/home/jody/bin/tabix"
	bed_file = "temp.bed"
	with open(bed_file,"w") as out:
		for chrom in pos_dict:
			for pos in sorted([int(x) for x in pos_dict[chrom]]):
				out.write("%s\t%s\t%s\n" % (chrom,pos,pos))
	cmd = "%s %s -h -T %s" % (tabix,matrix_file,bed_file)
	result = []
	call = subprocess.Popen(cmd,shell=True,stdout = subprocess.PIPE)
	for l in call.stdout:
		result.append(l)
	return result


def load_barcode(filename):
	temp = {}
	for l in open(filename):
		arr = l.rstrip().split()
		ref,alt = arr[3].split("/")
		temp[arr[1]] = {"lin":arr[0],"ref":ref,"alt":alt}
	return temp


barcode = load_barcode("/home/jody/refgenome/barcode.txt")
positions = {"Chromosome":[]}
for x in barcode:
	positions["Chromosome"].append(x)
results = defaultdict(list)
for l in matrix_filter(positions,sys.argv[1]):
	arr = l.rstrip().split()
	if l[0]=="#":
		header = arr
		continue
	chrom,pos = arr[:2]
	for i in xrange(5,len(arr)):
		if arr[i]==barcode[pos]["alt"]:
			results[header[i]].append(barcode[pos]["lin"])
for s in results:
	temp =  sorted(results[s])
	main = []
	for x in temp:
		if x in ["lineage1","lineage2","lineage3","lineage4","lineage5","lineage6","lineage7","lineageBOV","lineageBOV_AFRI"]:
			main.append(x)
	sub = sorted(list(set(temp)-set(main)))
	if len(main)==0:
		if len(sub)>0:
			main.append(sub[0].split(".")[0])
		else:
			main.append("-")
	if len(sub)==0:
		sub.append("-")
	print "%s\t%s\t%s" % (s,"%".join(main),"%".join(sub))
	
