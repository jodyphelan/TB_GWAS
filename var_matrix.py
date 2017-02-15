import sys
from collections import defaultdict
from tqdm import tqdm

class varMat:
	header = []
	samples = []
	matrix = defaultdict(dict)
	positions = []
	bedlines = []
	numlines = 0
	ann = defaultdict(lambda : defaultdict(dict))
	ref = defaultdict(dict)
	def __init__(self,filename):
		print "Loading matrix"
		with open(filename) as f:
			self.header = f.readline().rstrip().split()
			self.samples = self.header[5:]
			self.numlines+=1
			for l in tqdm(f):
				self.numlines+=1
				line = l.rstrip()
				arr = line.split()
				self.positions.append((arr[0],arr[1]))
				self.matrix[arr[0]][arr[1]] = line
				self.ref[arr[0]][arr[1]] = arr[2]
				if arr[3]!=".":
					for x in arr[3].split(";"):
						for a,b in x.split("="):
							self.ann[arr[0]][arr[1]][a] = b

	def write_bed_lines(self,outfile):
		print "Writing BED file"
		with open(outfile,"w") as o:
			for chrom,pos in tqdm(self.positions):
				o.write("%s\t%s\t%s\n" % (chrom,pos,pos))
	
	def binarise(self,outfile):
	        print "Converting to binary"
	        with open(outfile,"w") as o:
			o.write("\t".join(self.header)+"\n")
			for chrom,pos in tqdm(self.positions):
                                arr = self.matrix[chrom][pos].split()
                                for i in range(5,len(arr)):
                                        if arr[i]=="NA":
                                                continue
                                        elif arr[i]!=arr[2]:
                                                arr[i] = "1"
                                        else:
                                                arr[i] = "0"
                                o.write("\t".join(arr)+"\n")

	def load_variant_ann(self,infile,meta_name,col_num=3):
		col_num -= 1
		print "Loading variant annotation"
		for l in tqdm(open(infile)):
			arr = l.rstrip().split()
			chrom,pos = arr[:2]
			self.ann[arr[0]][arr[1]][meta_name] = arr[col_num]
			temp_arr = self.matrix[arr[0]][arr[1]].split()
			temp_ann_arr = []
			for x in sorted(self.ann[chrom][pos].keys()):
				temp_ann_arr.append("%s=%s" % (x,self.ann[chrom][pos][x]))
			ann_str = ";".join(temp_ann_arr)
			temp_arr[3] = ann_str
			self.matrix[chrom][pos] = "\t".join(temp_arr)

	def write_to_file(self,outfile):
		print "Writing matrix to file"
		with open(outfile,"w") as o:
			o.write("\t".join(self.header)+"\n")
			for chrom,pos in tqdm(self.positions):
				o.write("%s\n" % (self.matrix[chrom][pos]))

	def locus_sum(self,outfile,meta_name,fmt="matrix"):
		print "Aggregating mutations over loci"
		locus_sum = defaultdict(lambda :defaultdict(int))
		locus_pos = defaultdict(tuple)
		for chrom,pos in tqdm(self.positions):
			arr = self.matrix[chrom][pos].rstrip().split()
			locus_pos[self.ann[chrom][pos][meta_name]] = (arr[0],arr[1])
			for i in range(5,len(arr)):
				if arr[i]=="NA":
					continue
				elif arr[i]!=self.ref[chrom][pos]:
					locus_sum[self.ann[chrom][pos][meta_name]][self.header[i]]+=1
				else:
					locus_sum[self.ann[chrom][pos][meta_name]][self.header[i]]==0
		with open(outfile,"w") as o:	
			if fmt=="matrix":
				o.write("\t".join(self.header)+"\n")
				for locus in sorted(locus_pos.keys()):
					chrom,pos = locus_pos[locus]
					ref = "."
					info = "locus=%s" % locus
					type = "locus_sum"
					calls = "\t".join([str(locus_sum[locus][x]) for x in self.samples])
					o.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom,pos,ref,info,type,calls))
			elif fmt=="bimbam":
				for locus in sorted(locus_pos.keys()):
					temp_arr = [locus,"NA","NA"]
					for s in self.samples:
						temp_arr.append(str(locus_sum[locus][s]))
					o.write(",".join(temp_arr)+"\n")
	
	

	def write_bimbam(self,outfile):
		print "Writing bimbam"
		genofile = outfile 
		with open(genofile,"w") as o:
			for chrom,pos in tqdm(self.positions):
				temp_arr = []
				arr = self.matrix[chrom][pos].split()
				temp_arr.append("%s-%s" % (chrom,pos))
				temp_arr.append("NA")
				temp_arr.append("NA")
				for i in range(5,len(arr)):
					if arr[i]=="NA":
                                               	temp_arr.append("0") 
                                        elif arr[i]!=arr[2]:
                                                temp_arr.append("1")
                                        else:
                                                temp_arr.append("0")
				o.write(",".join(temp_arr)+"\n")


