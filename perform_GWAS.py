import sys
import var_matrix
import subprocess
from tqdm import tqdm
import	random 
from collections import defaultdict

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
			

if len(sys.argv)!=5:
	print "perform_GWAS.py <matrix> <metatfile> <meta_column> <prefix>"
	quit()

infile = sys.argv[1]
metafile = sys.argv[2]
selected_meta = sys.argv[3]
prefix = sys.argv[4]

binary_file = prefix+".mat.bin"
genes_file = prefix+".genes.txt"
variant_geno = prefix+".variants.geno"
locus_sum_geno = prefix+".locus_sum.geno"
pheno_file = prefix+".pheno"


mat_obj =  var_matrix.varMat(infile)
mat_obj.binarise(binary_file)
mat2genes(mat_obj,genes_file)
mat_obj.load_variant_ann(genes_file,"locus_tag")
mat_obj.load_variant_ann(genes_file,"gene_id",4)
mat_obj.locus_sum(locus_sum_geno,"locus_tag","bimbam")
mat_obj.write_bimbam(variant_geno)
write_pheno(mat_obj.samples,metafile,selected_meta,pheno_file)
subprocess.call("gemma -gk 1 -g %s -p %s -o %s" % (variant_geno,pheno_file,prefix),shell=True)
subprocess.call("gemma -g %s -p %s -k output/%s.cXX.txt -lmm 1 -o %s" % (variant_geno,pheno_file,prefix,prefix+".variants"),shell=True)
subprocess.call("gemma -g %s -p %s -k output/%s.cXX.txt -lmm 1 -o %s" % (locus_sum_geno,pheno_file,prefix,prefix+".locus_sum"),shell=True)
