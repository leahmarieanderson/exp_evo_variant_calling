"""
annotation.py

adapted from Chris Large, udated by Leah Anderson and Nicole Wang Sept 2023

takes input of VCF formatted file with mutation locations (from compare.py) and annotates SNPs

python3 annotation.py -f <VCF file containing mutations> -s <ORF sequences in fasta format> -n <non-coding GFF file for annotations> -g <genome fasta>

"""

#Chris edits were to make it accept common VCF4.2 headers, and accept variants with N (from Caiti's edits)

def main(vcf, orfs, noncoding_file, genome_file):

	"""GET DATA TOGETHER"""

	#gather data about genes:
	#orf id, start, stop (both 5'-3'), exons, introns, protein length, chr
	genes = {}
	for record in SeqIO.parse( open(orfs, 'r'), 'fasta' ):
		start = 0
		stop = 0
		exons = list(map(lambda x: x.split('-'), record.description.split(', ')[1].split(' ')[-1].split(',')))
		ch = chromosome_conversion(record.description.split('Chr ')[1].split(' ')[0])

		#reverse start and stop position of genes on Crick strand
		introns = []
		if record.id.split('-')[0][-1] == 'C':
			exons = list(exons)[::-1]
			start = int(exons[0][0])
			stop = int(exons[0][1])
			#find introns in genes that have them
			if len(exons) > 1:
				for e in range(len(exons)-1):
					introns.append([int(exons[e][1])+1, int(exons[e+1][0])-1])
		#repeat for Watston strand
		else:
			start = int(exons[0][0])
			stop = int(exons[0][1])
			if len(exons) > 1:
				for e in range(len(exons)-1):
					introns.append([int(exons[e][1])+1, int(exons[e+1][0])-1])
		
		genes[record.id] = [record.id] + [start, stop] + [exons] + [introns] + [len(record.seq)/3.0] + [ch]
		
	print ("done reading file")	
	#populate list of chromosomes in genome. genome fasta should have an identifier format of ">chrI" ">chrII" ">chrIII" etc
	genome = {}
	for record in SeqIO.parse(open(genome_file, 'r'), 'fasta'):
		genome[chromosome_conversion(record.description[:])] = str(record.seq)
	
	#populate second dictionary of non-coding annotations
	#noncoding[ID] = [ ID, chrom, regiontype, start, stop ]
	noncoding = {}
	key_int = 0
	for line in open(noncoding_file, 'r'):
		if line.startswith('#') != True:
			l = line.strip().split('\t')
			key_int = key_int + 1
			noncoding[key_int] = [ l[8].split(';')[0].split('=')[1], chromosome_conversion(l[0][3:]), l[2], int(l[3]), int(l[4]) ]
		
	"""ANNOTATE SNPS"""

	#open output file
	f_out = open(vcf[0:-21]+"annotated_vcf.txt", 'w')

	#start reading vcf file
	#chr, start, stop, ref, obs
	Chrom = 0
	Pos = 1
	Ref = 3
	Alt = 4
	for line in open(vcf, 'r'):
		#if it's a header, then print out more header to the next line
		if 'POS' in line:
			f_out.write(line.strip() + '\t' + '\t'.join(['ANNOTATION', 'REGION', 'PROTEIN']) + '\n')
			l = line.strip().split('\t')
			continue
		elif '##' in line:
			f_out.write(line.strip()+'\n')
			continue
		else:
			l = line.strip().split('\t')
			if len(l[Alt]) > 1:
				base = 'X'
				print ('indel found!')
				print (l[Alt])
			elif len(l[Ref]) >1:
				base = 'X'
				print ('indel found!')
				print (l[Alt])
			elif len(l[Alt]) == 1:
				base = l[Alt]
				print ('SNP found')
				print (l[3])
		print (str(l[Chrom]) + " " + str(l[Pos]))
		snp_pos = int(l[Pos])
		
# 		#make copy of chromosome with mutation
		mut_chr = list(copy(genome[chromosome_conversion(l[Chrom])]))
		mut_chr[snp_pos-1] = base
				
		annotation = False
		
		#loop through genes, trying to find one containing snp
		for g in genes:
			#if gene on correct chromosome
			if genes[g][6] == chromosome_conversion(l[Chrom]):
				#CRICK STRAND
				if genes[g][0].split('-')[0][-1] == 'C':
					#found gene containing snp
					if snp_pos <= genes[g][1] and snp_pos >= genes[g][2]:
						#check if snp is in an intron:
						for intron in genes[g][4]:
							if snp_pos <= intron[0] and snp_pos >= intron[1]:
								#found an intronic snp - check if it's near a splice site
								if snp_pos in range(intron[0]-2, intron[0]+1) or snp_pos in range(intron[1],intron[1]+3):
									# NICOLE'S CHANGES: fixed below syntax
									print('\t'.join(l + ['splice-site', genes[g][0], 'NA']), file=f_out)
								#otherwise it's just intronic (boring!)
								else:
									# NICOLE'S CHANGES: fixed below syntax
									print('\t'.join(l + ['intron', genes[g][0], 'NA']), file=f_out)
 						#if the snp is within gene start-stop, but not in an intron, then it's in an exon
 						#remove intronic sequences from gene
						mut_gene = []
						wt_gene = []
						
						for exon in genes[g][3]:
							#check first if the snp is near a splice site
							if snp_pos in range(int(exon[0])-2,int(exon[0])+1) or snp_pos in range(int(exon[1]),int(exon[1])+3):
								#if it's not in the start/stop regions (these aren't splice-sites)
								if snp_pos not in range(genes[g][1]-2,genes[g][1]+1) and snp_pos not in range(genes[g][2],genes[g][2]+3):
									# NICOLE'S CHANGES: fixed below syntax
									print('\t'.join(l + ['splice-site', genes[g][0], 'NA']), file=f_out)
								
							mut_gene += list(complement(mut_chr[ int(exon[1])-1:int(exon[0]) ])[::-1])
							wt_gene += list(complement(genome[genes[g][6]][ int(exon[1])-1:int(exon[0]) ])[::-1])
													
						#loop through codons, find mismatch
						for codon in range(0, len(mut_gene), 3):
							if mut_gene[codon:codon+3] != wt_gene[codon:codon+3]:
								mut_aa = lookup_codon(''.join(mut_gene[codon:codon+3]))
								wt_aa = lookup_codon(''.join(wt_gene[codon:codon+3]))
								if mut_aa != wt_aa:
									# NICOLE'S CHANGES: changed below line to match python 3 syntax
									print('\t'.join(l + ['coding-nonsynonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa]), file=f_out)
								else:
									# NICOLE'S CHANGES: fixed below syntax
									print('\t'.join(l + ['coding-synonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa]), file=f_out)
					
						annotation = True
						print ('Crick')
					#if snp_pos isn't in a gene, check if it's upstream of a gene
					elif snp_pos in range(genes[g][1]+1,genes[g][1]+201):
						# NICOLE'S CHANGES: changed below line to match python 3 syntax
						print('\t'.join(l + ["5'-upstream", genes[g][0], 'NA']), file=f_out)
						annotation = True	
						print ('crick 5\'')
				#WATSON STRAND
				else:
					#found gene containing snp
					if snp_pos >= genes[g][1] and snp_pos <= genes[g][2]:

						#check if snp is in an intron:
						for intron in genes[g][4]:
							if snp_pos >= intron[0] and snp_pos <= intron[1]:
								#found an intronic snp - check if it's splice site
								if snp_pos in range(intron[0],intron[0]+3) or snp_pos in range(intron[1]-2,intron[1]+1):
									# NICOLE'S COMMENTS: fixed below syntax
									print('\t'.join(l + ['splice-site', genes[g][0], 'NA']), file=f_out)
								#otherwise, it's just intronic
								else:
									# NICOLE'S COMMENTS: fixed below syntax
									print('\t'.join(l + ['intron', genes[g][0], 'NA']), file=f_out)
								annotation = True	
								print ('watson')
						#if the snp is within the gene start-stop, but not in an intron, then it's in an exon
						#remove intronic sequences
						mut_gene = []
						wt_gene = []
						
						for exon in genes[g][3]:
							#first, check if in a splice-site
							if snp_pos in range(int(exon[0]),int(exon[0])+3) or snp_pos in range(int(exon[1])-2,int(exon[1])+1):
								#if it's not the start/stop
								if snp_pos not in range(genes[g][1],genes[g][1]+3) and snp_pos not in range(genes[g][2]-2,genes[g][2]+1):
									# NICOLE'S COMMENTS: fixed below syntax
									print('\t'.join(l + ['splice-site', genes[g][0], 'NA']), file=f_out)
								
							mut_gene += list(mut_chr[ int(exon[0])-1:int(exon[1]) ])
							wt_gene += list(genome[genes[g][6]][int(exon[0])-1:int(exon[1])])

						#loop through codons, find mismatch
						for codon in range(0, len(mut_gene), 3):
							if mut_gene[codon:codon+3] != wt_gene[codon:codon+3]:
								mut_aa = lookup_codon(''.join(mut_gene[codon:codon+3]))
								wt_aa = lookup_codon(''.join(wt_gene[codon:codon+3]))
								#check if synonymous, non-synonymous
								if mut_aa != wt_aa:
									# NICOLE'S CHANGES: below line is the updated python 3 syntax
									print('\t'.join(l + ['coding-nonsynonymous', genes[g][0], wt_aa + str(codon /3+1) + mut_aa]), file=f_out)
								else:
									# NICOLE'S CHANGES: below line is the updated python 3 syntax
									print('\t'.join(l + ['coding-synonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa]), file=f_out)

								annotation = True
								print ('watson')
						
						annotation = True
									
					#if snp isn't in a gene, check if it's in the upstream region of a gene
					elif snp_pos in range(genes[g][1]-200,genes[g][1]):
						# NICOLE'S COMMENTS: fixed below syntax
						print('\t'.join(l + ["5'-upstream", genes[g][0], 'NA']), file=f_out)
						annotation = True						
						print ('5')
		#after checking all genes for genic or non-coding snps, check non-coding elements
		if annotation == False:
			print ('watson and crick')	
			#if snp isn't upstream of a gene, check if it's in a non-coding region of the genome
			for n in noncoding:
				#if correct chromosome
				if noncoding[n][1] == chromosome_conversion(l[0]):
					if snp_pos >= noncoding[n][3] and snp_pos <= noncoding[n][4]:
						# NICOLE'S COMMENTS: fixed below syntax
						print('\t'.join(l + [noncoding[n][2], noncoding[n][0], 'NA']), file=f_out)
						annotation = True		
						print ('non-coding')
			#if annotation is still false after check non-coding elements, then the snp is just intergenic
			if annotation == False:
				# NICOLE'S COMMENTS: fixed below syntax
				print('\t'.join(l + ['intergenic', 'NA', 'NA']), file=f_out)
									
	"""OTHER USEFUL FUNCTIONS"""

def chromosome_conversion(chrom_number):
	chrom_conv = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 
				'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16, 'Mito':17, 'mitochondrion':17, 'M':17, '0M':17}
	if chrom_number.startswith('chr'):
		chrom_number = chrom_number[3:]
	
	try:
		if int(chrom_number) in chrom_conv.values():
			return int(chrom_number)
	except ValueError:			
		return chrom_conv[chrom_number]

def complement(base):
	compbase = []
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'X':'X', 'N':'N' }
	for i in range(len(base)):
		compbase.append(comp[base[i].upper()])
	return ''.join(compbase)

def lookup_codon(codon):
	lookup = { 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C', 'tnt': 'missing',
			'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C', 'tnc': 'missing',
			'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*', 'tna': 'missing',
			'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W', 'tng': 'missing',
			'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R', 'cnt': 'missing',
			'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R', 'cnc': 'missing',
			'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R', 'cna': 'missing',
			'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R', 'cng': 'missing',
			'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S', 'ant': 'missing',
			'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S', 'anc': 'missing',
			'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R', 'ana': 'missing',
			'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R', 'ang': 'missing',
			'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G', 'gnt': 'missing',
			'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G', 'gnc': 'missing',
			'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G', 'gna': 'missing',
			'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G' , 'gng': 'missing',
			'tnn': 'missing', 'cnn': 'missing', 'ann': 'missing', 'gnn': 'missing', 'nnn': 'missing',
			'nnt': 'missing', 'nnc': 'missing', 'nna': 'missing', 'nng': 'missing',
			'ttn': 'missing', 'tcn': 'S', 'tan': 'missing', 'tgn': 'missing',
			'ctn': 'L', 'ccn': 'P', 'can': 'missing', 'cgn': 'R',
			'atn': 'missing', 'acn': 'T', 'aan': 'missing', 'agn': 'missing',
			'gtn': 'V', 'gcn': 'A', 'gan': 'missing', 'ggn': 'G',
			'aax': 'indel', 'agx': 'indel', 'acx': 'indel', 'anx': 'indel', 'axn': 'indel',
			'atx': 'indel', 'axa': 'indel', 'axg': 'indel', 'nxa': 'indel', 'xna': 'indel',
			'axc': 'indel', 'axt': 'indel', 'axx': 'indel', 'cnx': 'indel', 'cxn': 'indel',
			'gax': 'indel', 'ggx': 'indel', 'gcx': 'indel', 'nxc': 'indel', 'xnc': 'indel',
			'gtx': 'indel', 'gxa': 'indel', 'gxg': 'indel', 'gnx': 'indel', 'gxn': 'indel',
			'gxc': 'indel', 'gxt': 'indel', 'gxx': 'indel', 'nxg': 'indel', 'xng': 'indel',
			'cax': 'indel', 'cgx': 'indel', 'ccx': 'indel', 'tnx': 'indel', 'txn': 'indel',
			'ctx': 'indel', 'cxa': 'indel', 'cxg': 'indel', 'nxt': 'indel', 'xnt': 'indel',
			'cxc': 'indel', 'cxt': 'indel', 'cxx': 'indel', 'xan': 'indel', 'nax': 'indel',
			'tax': 'indel', 'tgx': 'indel', 'tcx': 'indel', 'xcn': 'indel', 'ncx': 'indel',
			'ttx': 'indel', 'txa': 'indel', 'txg': 'indel', 'xgn': 'indel', 'ngx': 'indel',
			'txc': 'indel', 'txt': 'indel', 'txx': 'indel', 'xtn': 'indel', 'ntx': 'indel',
			'xaa': 'indel', 'xag': 'indel', 'xac': 'indel',
			'xat': 'indel', 'xax': 'indel', 'xga': 'indel',
			'xgg': 'indel', 'xgc': 'indel', 'xgt': 'indel', 'xnn': 'indel',
			'xgx': 'indel', 'xca': 'indel', 'xcg': 'indel', 'xxn': 'indel',
			'xcc': 'indel', 'xct': 'indel', 'xcx': 'indel', 'xnx': 'indel',
			'xta': 'indel', 'xtg': 'indel', 'xtc': 'indel', 'nxx': 'indel',
			'xtt': 'indel', 'xtx': 'indel', 'xxa': 'indel', 'nnx': 'indel',
			'xxg': 'indel', 'xxc': 'indel', 'xxt': 'indel', 'xxx': 'indel'}
	return lookup[codon.lower()]

# translate DNA -> amino acid
def translate_sequence(seq):
	translated_seq = ''
	i = 0
	while i <= len(seq)-3:
		translated_seq += lookup_codon(seq[i:i+3])
		i += 3
	return translated_seq

if __name__ == "__main__":
	from optparse import OptionParser
	import sys
	import Bio
	from Bio import SeqIO
	from copy import copy
	from random import randint
	
	parser = OptionParser()
	parser.add_option('-f', '--input', action = 'store', type = 'string', dest = 'inputfile', help = 'file with mutations')
	parser.add_option('-s', '--sequences', action = 'store', type = 'string', dest = 'sequences', help = 'fasta file of coding sequences')
	parser.add_option('-n', '--non-coding', action = 'store', type = 'string', dest = 'noncoding', help = 'gff file containing non-coding regions')
	parser.add_option('-g', '--genome', action = 'store', type = 'string', dest = 'genome', help = 'fasta file containing genome sequence')
	(option, args) = parser.parse_args()

	main(option.inputfile, option.sequences, option.noncoding, option.genome)