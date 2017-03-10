import sys
import var_matrix
import subprocess
from tqdm import tqdm
import	random 
from collections import defaultdict
import argparse

def mat2genes(matObj,outfile):

	tempfile = "temp."+str(random.randint(0,1000))+".bed"
	matObj.write_bed_lines(tempfile)
	cmd = subprocess.Popen("tabix /opt/storage/mtuberculosis/tb_scripts/MTB-h37rv_asm19595v2-eg18.tab.ann.gz -T %s" % tempfile, shell=True,stdout=subprocess.PIPE)
	print "Annotating positions"
	with open(outfile,"w") as o:
		for l in tqdm(cmd.stdout):
			arr = l.rstrip().split()
			o.write("\t".join([arr[0],arr[1],arr[15],arr[16]])+"\n")
	subprocess.call("rm %s" % tempfile,shell=True)

def mat2vartype(matObj,outfile):
	tempfile = "temp."+str(random.randint(0,1000))+".bed"
	matObj.write_bed_lines(tempfile)
	cmd = subprocess.Popen("tabix /opt/storage/mtuberculosis/tb_scripts/MTB-h37rv_asm19595v2-eg18.tab.ann.gz -T %s" % tempfile, shell=True,stdout=subprocess.PIPE)
	print "Loading annotation"
	ann_dict = defaultdict(dict)
	for l in tqdm(cmd.stdout):
		arr = l.rstrip().split()
		nuc_obj = {}
		nuc_obj[arr[2]] = {"codon":arr[6],"aa":arr[11]}
		nuc_obj[arr[3]] = {"codon":arr[8],"aa":arr[12]}
		nuc_obj[arr[4]] = {"codon":arr[9],"aa":arr[13]}
		nuc_obj[arr[5]] = {"codon":arr[10],"aa":arr[14]}
		ann_dict[arr[0]][arr[1]] = {"ref_nt":arr[2],"ref_codon":arr[6],"ref_aa":arr[11],"chr":arr[0],"pos":arr[1],"ann":nuc_obj,"rv":arr[15],"gene":arr[16],"gene_syn":arr[17],"ncr":arr[18],"start":arr[19],"end":arr[20],"strand":arr[21],"drug":arr[22],"ppe":arr[23],"codon_num":arr[24],"gene_nt":arr[25],"operon":arr[26]}
	with open(outfile,"w") as o:
		print "Assigning variant types"
		for chrom,pos in tqdm(matObj.positions):
			arr = matObj.matrix[chrom][pos].split()	
			calls = set(arr[5:])
			if "N" in calls:
				calls.remove("N")
			if "NA" in calls:
				calls.remove("NA")
			if "-" in calls:
				calls.remove("-")
			if ann_dict[chrom][pos]["ref_nt"] in calls:
				calls.remove(ann_dict[chrom][pos]["ref_nt"])
			metaset = set()
			for c in calls:
				if "+" in c or "-" in c:
					if ann_dict[chrom][pos]["ncr"]=="CDS":
						metaset.add((c,"CDS_NS"))
					else:
						metaset.add((c,"inter"))
				elif ann_dict[chrom][pos]["ann"][c]["aa"] != ann_dict[chrom][pos]["ref_aa"]:
					metaset.add((c,"CDS_NS"))
				elif ann_dict[chrom][pos]["ann"][c]["aa"] == ann_dict[chrom][pos]["ref_aa"]:
					if ann_dict[chrom][pos]["ncr"]=="CDS":
						metaset.add((c,"CDS_S"))
					else:
						metaset.add((c,"inter"))
			o.write("%s\t%s\t%s\n" % (chrom,pos,",".join([x+":"+y for x,y in metaset])))
	subprocess.call("rm %s" % tempfile,shell=True)

	
def loadTSV(infile):
        tsv_dict = defaultdict(dict)
        with open(infile) as f:
                header = f.readline().split()
                for l in f:
                        arr = l.rstrip().split()
                        for i in range(1,len(arr)):
                                tsv_dict[arr[0]][header[i]] = arr[i]
        return tsv_dict,header[1:]

def write_pheno(samples,metafile,meta_name,outfile):
	tsv_dict,colnames = loadTSV(metafile)
	with open(outfile,"w") as o:
		for s in samples:
			o.write("%s\n" % tsv_dict[s][meta_name])
			
def write_samples(matObj,outfile):
	with open(outfile,"w") as o:
		for s in matObj.samples:
			o.write("%s\n" % s)

def load_samples(samplefile):
	return [x.rstrip() for x in open(samplefile).readlines()]

def parse_variant_assoc(prefix,filename):
	temp_ann = defaultdict(dict)
	for l in open(prefix+".genes.txt"):
		arr = l.rstrip().split()
		temp_ann[arr[0]][arr[1]] = (arr[2],arr[3])
	outlines = []
	for l in open(filename):
		arr = l.rstrip().split()
		if arr[0]=="chr":
			arr.append("locus_tag")
			arr.append("gene")
				
		else:
			temp_arr = arr[1].split("-")
			arr.append(temp_ann[temp_arr[0]][temp_arr[1]][0])
			arr.append(temp_ann[temp_arr[0]][temp_arr[1]][1])
		outlines.append(arr)
	with open(filename,"w") as o:
		for arr in outlines:
			o.write("%s\n" % ("\t".join(arr)))

	
	
####################### Main functions ############################

def main_preprocess(args):

	infile = args.matrix
	prefix = args.prefix

	binary_file = prefix+".mat.bin"
	genes_file = prefix+".genes.txt"
	vartype_file = prefix+".vartype.txt"
	variant_geno = prefix+".variants.geno"
	locus_sum_geno = prefix+".locus_sum.geno"
	sample_file = prefix+".samples"

	mat_obj =  var_matrix.varMat(infile)
	write_samples(mat_obj,sample_file)
	mat_obj.binarise(binary_file)
	mat2genes(mat_obj,genes_file)
	mat2vartype(mat_obj,vartype_file)
	mat_obj.load_variant_ann(genes_file,"locus_tag")
	mat_obj.load_variant_ann(genes_file,"gene_id",4)
	mat_obj.load_variant_ann(vartype_file,"var_type")
	mat_obj.locus_sum(locus_sum_geno,"locus_tag","bimbam")
	mat_obj.write_bimbam(variant_geno)


def main_assoc(args):

	prefix = args.prefix
	outprefix = args.outprefix
	metafile = args.metafile
	selected_meta = args.metacolumn

	pheno_file = outprefix+".pheno"
	sample_file = prefix+".samples"
	variant_geno = prefix+".variants.geno"
	locus_sum_geno = prefix+".locus_sum.geno"
	
	samples = load_samples(sample_file)
	write_pheno(samples,metafile,selected_meta,pheno_file)
#	subprocess.call("gemma -gk 1 -g %s -p %s -o %s" % (variant_geno,pheno_file,prefix),shell=True)
#	subprocess.call("gemma -g %s -p %s -k output/%s.cXX.txt -lmm 1 -o %s" % (variant_geno,pheno_file,prefix,outprefix+".variants"),shell=True)
#	subprocess.call("gemma -g %s -p %s -k output/%s.cXX.txt -lmm 1 -o %s" % (locus_sum_geno,pheno_file,prefix,outprefix+".locus_sum"),shell=True)
	parse_variant_assoc(prefix,"output/%s.variants.assoc.txt" % outprefix)

####################### Argument Parser ###########################

parser = argparse.ArgumentParser(description='Perform GWAS')
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('preprocess', help='Generate files for association', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('matrix', help='Initial matrix')
parser_sub.add_argument('prefix', help='Prefix for output files')
parser_sub.set_defaults(func=main_preprocess)

parser_sub = subparsers.add_parser('assoc', help='Generate files for association', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('prefix', help='Prefix for input files')
parser_sub.add_argument('outprefix', help='Prefix for output files')
parser_sub.add_argument('metafile', help='Metadata TSV file with first column as the sample name and the rest phenotypes')
parser_sub.add_argument('metacolumn', help='Column to use in association')
parser_sub.set_defaults(func=main_assoc)


args = parser.parse_args()
args.func(args)

