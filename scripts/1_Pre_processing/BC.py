"""
Sasha Levy, public domain. 2016
Code for parsing fastq files and counting barcode frequencicies 
Bartender must be installed (https://github.com/LaoZZZZZ/bartender-1.1)
"""
__author__ = "Sasha Levy"
__date__ = "Jan 2016"
__version__ = "1.0.0"


#Formats of input variables
directory = "/Users/sasha/Sequencing/" #directory where files are read and saved
f_gzipped_fastqfile = "ZLPCA06_NoIndex_L007_R1_all.fastq.gz" #The forward reads, gzip file
r_gzipped_fastqfile = "ZLPCA06_NoIndex_L007_R3_all.fastq.gz" #The reverse reads, gzip file
q = "fastq" #the type of fastq file coming off of the sequencer, what format are the quality scores?
f_seqtag_length = 8 #the length of the sequencing tag on the first read (UMI1)
r_seqtag_length = 8 #the length of the sequencing tag on the second read (UMI2)
f_multitag_length = 6 #the length of the multiplexing tag on the first read
r_multitag_length = 6 #the length of the multiplexing tag on the second read
f_lintag_length = 38 #the length of the lineage tag on the first read (first barcode)
r_lintag_length = 38 #the length of the lineage tag on the second read (second barcode)
f_spacer_length = 43 #distance to first barcode in forward read, not including the multitag and the seqtag
r_spacer_length = 29 #distance second barcode in reverse read, not including the multitag and the seqtag
min_qs = 30 #the minimum avareage quality score for both lineage tags
lintag_grep_filter1 ='\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*' #first barcode
lintag_grep_filter2 ='\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*' #second barcode
clip_ends = 1 #logical of whether or not to clip the front and back ends off of lintag1 and lintag2
lintag1_front_clipper = '(.ACC|T.CC|TA.C|TAC.)' #only report lintag1 after this sequence
lintag2_front_clipper = '(.ACC|T.CC|TA.C|TAC.)' #only report lintag2 after this sequence
lintag1_rear_clipper = '(.ATA|A.TA|AA.A|AAT.)' #only report lintag1 before this sequence, this must be the COMPLIMENT of the true sequence
lintag2_rear_clipper = '(.ATA|A.TA|AA.A|AAT.)' #only report lintag2 before this sequence, this must be the COMPLIMENT of the true sequence
multitags = ["TAGCTTGCGTAC", "CGATGTGAGACG"] #concatenated multiplexing tags from the first and second reads that uniquely identify a sample, currently must have 2 or more multitags
write_multitags = False #True will write multitags for all reads, otherwise only multitags for unmatched reads will be written
merged_output = True #True will merge double barcodes into one barcode and treat these as single barcodes
cores = 8 #the number of processors
frequency_cutoff = 1 #clusters less than this number of reads under this cutoff will not be reported by bartender


#must add this module to the python path if not already added
import sys
sys.path.append("/Users/sasha/Dropbox/code/BC")
import BC

#The main function
def count(directory, f_gzipped_fastqfile, r_gzipped_fastqfile,
	q = "fastq",
	f_seqtag_length = 8,
	r_seqtag_length = 8,
	f_multitag_length = 6,
	r_multitag_length = 6,
	f_lintag_length = 38,
	r_lintag_length = 38, 
	f_spacer_length = 43, #distance to first barcode in forward read (ignoring the length the thie multitag and the seqtag)
	r_spacer_length = 29, #distance second barcode in reverse read (ignoring the length the thie multitag and the seqtag)
	min_qs = 30, #the minimum avareage quality score for both lineage tags
	lintag_grep_filter1 ='\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*', #first barcode
	lintag_grep_filter2 ='\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*', #second barcode
	clip_ends = 1, #logical of whether or not to clip the front and back ends off of lintag1 and lintag2
	lintag1_front_clipper = '(.ACC|T.CC|TA.C|TAC.)', #only report lintag1 after this sequence
	lintag2_front_clipper = '(.ACC|T.CC|TA.C|TAC.)', #only report lintag2 after this sequence
	lintag1_rear_clipper = '(.TAA|A.AA|AT.A|ATA.)', #only report lintag1 before this sequence, this must be the COMPLIMENT of the true sequence
	lintag2_rear_clipper = '(.TAA|A.AA|AT.A|ATA.)', #only report lintag2 before this sequence, this must be the COMPLIMENT of the true sequence
	multitags = '', #concatenated multiplexing tags from the first and second reads that uniquely identify a sample, currently must have 2 or more multitags
	write_multitags = False, #True will write multitags for all reads, otherwise only multitags for unmatched reads will be written
	merged_output = True, #True will merge double barcodes into one barcode and treat these as single barcodes
	cores = 8, #the number of processors
	frequency_cutoff = 1): #clusters less than this number of reads under this cutoff will not be reported by bartendder
    
        import os
        from itertools import izip
        import BC
        if len(multitags) < 2:
            parse_fastq(directory, f_gzipped_fastqfile, r_gzipped_fastqfile,
            q = q, 	f_seqtag_length = f_seqtag_length, r_seqtag_length = r_seqtag_length,
            f_multitag_length = f_multitag_length, r_multitag_length = r_multitag_length,
            f_lintag_length = f_lintag_length, r_lintag_length = r_lintag_length,
            f_spacer_length = f_spacer_length, r_spacer_length = r_spacer_length,
            min_qs = min_qs, lintag_grep_filter1 = lintag_grep_filter1,
            lintag_grep_filter2 = lintag_grep_filter2, 
            clip_ends = clip_ends, lintag1_front_clipper = lintag1_front_clipper,
            lintag2_front_clipper = lintag2_front_clipper,
            lintag1_rear_clipper = lintag1_rear_clipper,
            lintag2_rear_clipper = lintag2_rear_clipper)
            #merge barcodes for bartender input
            l1 = directory + '_lintag1.txt'
            l2 = directory  + '_lintag2.txt'
            umi = directory + '_seqtag.txt'
            i = ''
            if merged_output:
                os.system("paste -d \'\\0,\' " + l1 + " " + l2 + " " + umi + " > " + directory + i + "_merged.txt")
            else:
                os.system("perl -lne \'print $_, \",\", $.\' " + l1 + " > " + directory + i + "_BC1.txt")
                os.system("perl -lne \'print $_, \",\", $.\' " + l2 + " > " + directory + i + "_BC2.txt")
            #run bartenter
            if merged_output:
                os.system("bartender_single_com -f " + directory + i + "_merged.txt" +  " -o " + directory + i + "_merged -c " + str(frequency_cutoff) + " -t " + str(cores) + " > " + directory + i + "_bartentder__merged_stdout.txt")
            else:
                os.system("bartender_single_com -f " + directory + i + "_BC1.txt" +  " -o " + directory + i + "_BC1 -c " + str(frequency_cutoff) + " -t " + str(cores) + " > " + directory + i + "_bartentder_BC1_stdout.txt")
                os.system("bartender_single_com -f " + directory + i + "_BC2.txt" +  " -o " + directory + i + "_BC2 -c " + str(frequency_cutoff) + " -t " + str(cores) + " > " + directory + i + "_bartentder_BC2_stdout.txt")
                        
            #For split runs, merge DBCs, and dedupe
            import csv
            if not merged_output:
                with open(directory +  '_BC1_barcode.csv', 'rb') as f:
                    reader = csv.reader(f)
                    b1 = list(reader)
                with open(directory +  '_BC2_barcode.csv', 'rb') as f:
                    reader = csv.reader(f)
                    b2 = list(reader)
                with open(directory + i +  '_BC1_cluster.csv', 'rb') as f:
                    reader = csv.reader(f)
                    b3= list(reader)
                with open(directory + i +  '_BC2_cluster.csv', 'rb') as f:
                    reader = csv.reader(f)
                    b4= list(reader)
        
                #make dictionaries
                d1 = dict(zip([a for a,b,c in b1], [c for a,b,c in b1]))
                d2 = dict(zip([a for a,b,c in b2], [c for a,b,c in b2]))
                d3 = dict(zip([a for a,b,c,d in b3], [b for a,b,c,d in b3]))
                d4 = dict(zip([a for a,b,c,d in b4], [b for a,b,c,d in b4]))
                #count
                bc1_file = open(directory + i + "_lintag1.txt", "r")
                bc2_file = open(directory + i + "_lintag2.txt", "r")
                umi_file = open(directory + i + "_seqtag.txt", "r")
        
                #merge barcodes and count duplicates
                d = dict(); e = dict(); #d is the count of each DBC, e is a list of UMIs for each DBC
                for bc1, bc2, umi in izip(bc1_file, bc2_file, umi_file):
                    bc1 = bc1.strip()
                    bc2 = bc2.strip()
                    umi = umi.strip()
                    if bc1 in d1 and bc2 in d2:
                        y = d3[d1[bc1]] + '_' + d4[d2[bc2]]
                        d[y] = d.get(y,0)+1
                        j = e.get(y)
                        if j == None:
                            e[y] = [umi]
                        else:
                            j.append(umi)
                            e[y] = j
                
                #write to file
                m = list([['BC1cluster_BC2cluster', 'all_reads', 'deduped_reads']])
                for k, v in d.iteritems():
                    temp = [k,str(v), str(len(set(e[k])))]
                    m.append(temp)
        
                with open(directory +  '_DBCcounts.csv', 'w') as file:
                    file.writelines(','.join(a) + '\n' for a in m)
        else:
            parse_fastq_by_multitag(directory, f_gzipped_fastqfile, r_gzipped_fastqfile,
            q = q, 	f_seqtag_length = f_seqtag_length, r_seqtag_length = r_seqtag_length,
            f_multitag_length = f_multitag_length, r_multitag_length = r_multitag_length,
            f_lintag_length = f_lintag_length, r_lintag_length = r_lintag_length,
            f_spacer_length = f_spacer_length, r_spacer_length = r_spacer_length,
            min_qs = min_qs, lintag_grep_filter1 = lintag_grep_filter1,
            lintag_grep_filter2 = lintag_grep_filter2, 
            clip_ends = clip_ends, lintag1_front_clipper = lintag1_front_clipper,
            lintag2_front_clipper = lintag2_front_clipper,
            lintag1_rear_clipper = lintag1_rear_clipper,
            lintag2_rear_clipper = lintag2_rear_clipper, multitags = multitags, 
            write_multitags = write_multitags)
            #merge barcodes for bartender input
            for i in multitags:
                l1 = directory + i + '_lintag1.txt'
                l2 = directory + i + '_lintag2.txt'
                umi = directory + i + '_seqtag.txt'
                if merged_output:
                    os.system("paste -d \'\\0,\' " + l1 + " " + l2 + " " + umi + " > " + directory + i + "_merged.txt")
                else:
                    #os.system("paste -d \',\' " + l1 + " " + umi + " > " + directory + i + "_BC1.txt")
                    #os.system("paste -d \',\' " + l2 + " " + umi + " > " + directory + i + "_BC2.txt")
                    os.system("perl -lne \'print $_, \",\", $.\' " + l1 + " > " + directory + i + "_BC1.txt")
                    os.system("perl -lne \'print $_, \",\", $.\' " + l2 + " > " + directory + i + "_BC2.txt")
            #run bartenter
            for i in multitags:
                if merged_output:
                    os.system("bartender_single_com -f " + directory + i + "_merged.txt" +  " -o " + directory + i + "_merged -c " + str(frequency_cutoff) + " -t " + str(cores) + " > " + directory + i + "_bartentder__merged_stdout.txt")
                else:
                    os.system("bartender_single_com -f " + directory + i + "_BC1.txt" +  " -o " + directory + i + "_BC1 -c " + str(frequency_cutoff) + " -t " + str(cores) + " > " + directory + i + "_bartentder_BC1_stdout.txt")
                    os.system("bartender_single_com -f " + directory + i + "_BC2.txt" +  " -o " + directory + i + "_BC2 -c " + str(frequency_cutoff) + " -t " + str(cores) + " > " + directory + i + "_bartentder_BC2_stdout.txt")
                    
            #For split runs, merge DBCs, and dedupe
            import csv
            if not merged_output:
                for i in multitags:
                    with open(directory + i +  '_BC1_barcode.csv', 'rb') as f:
                        reader = csv.reader(f)
                        b1 = list(reader)
                    with open(directory + i +  '_BC2_barcode.csv', 'rb') as f:
                        reader = csv.reader(f)
                        b2 = list(reader)
                    with open(directory + i +  '_BC1_cluster.csv', 'rb') as f:
                        reader = csv.reader(f)
                        b3= list(reader)
                    with open(directory + i +  '_BC2_cluster.csv', 'rb') as f:
                        reader = csv.reader(f)
                        b4= list(reader)
                        
                    #make dictionaries
                    d1 = dict(zip([a for a,b,c in b1], [c for a,b,c in b1]))
                    d2 = dict(zip([a for a,b,c in b2], [c for a,b,c in b2]))
                    d3 = dict(zip([a for a,b,c,d in b3], [b for a,b,c,d in b3]))
                    d4 = dict(zip([a for a,b,c,d in b4], [b for a,b,c,d in b4]))
            
                    #count
                    bc1_file = open(directory + i + "_lintag1.txt", "r")
                    bc2_file = open(directory + i + "_lintag2.txt", "r")
                    umi_file = open(directory + i + "_seqtag.txt", "r")
            
                    #merge barcodes and count duplicates
                    d = dict(); e = dict(); #d is the count of each DBC, e is a list of UMIs for each DBC
                    for bc1, bc2, umi in izip(bc1_file, bc2_file, umi_file):
                        bc1 = bc1.strip()
                        bc2 = bc2.strip()
                        umi = umi.strip()
                        if bc1 in d1 and bc2 in d2:
                            y = d3[d1[bc1]] + '_' + d4[d2[bc2]]
                            d[y] = d.get(y,0)+1
                            j = e.get(y)
                            if j == None:
                                e[y] = [umi]
                            else:
                                j.append(umi)
                                e[y] = j
                    #write to file
                    m = list([['BC1cluster_BC2cluster', 'all_reads', 'deduped_reads']])
                    for k, v in d.iteritems():
                        temp = [k,str(v), str(len(set(e[k])))]
                        m.append(temp)
            
                    with open(directory + i +  '_DBCcounts.csv', 'w') as file:
                        file.writelines(','.join(a) + '\n' for a in m)

#Small support functions = gzip_mapcount, grep, best_match, mismatches
def gzip_mapcount(filename):
	"""counts the number of lines in a gz file"""
	import mmap
	import gzip
	f = gzip.open(filename, "r+")
	lines = 0
	buf_size = 1024 * 1024
	read_f = f.read # loop optimization
	buf = read_f(buf_size)
	while buf:
	   lines += buf.count('\n')
	   buf = read_f(buf_size)

	return lines

def grep(tag, regex):
	"""returns 1 if grep match or 0 if no match"""
	import re
	r = re.compile(regex) 
	s = r.match(tag)
	if s == None:
		return 0
	else:
		return 1

def mismatches(SEQ1, SEQ2, MAX = float("inf"), IGNORE_N = 0 ):
	"""Returns the number of mismatches between two strings.
	MAX sets the number of max number of mismatches that are reported.
	Lowering MAX increases performance.
	IGNORE_N = 1 will ignore mismatches with N."""
	mismatches = 0
	if SEQ1 != SEQ2: #first check for exact match
		if IGNORE_N == 0:
			for i in range(len(SEQ1)):
				if SEQ1[i] != SEQ2[i]:
					mismatches += 1
				if mismatches >= MAX:
					break
			return mismatches
		else:
			for i in range(len(SEQ1)):
				if SEQ1[i] != 'N' and SEQ2[i] != 'N':
					if SEQ1[i] != SEQ2[i]:
						mismatches += 1
				if mismatches >= MAX:
					break
			return mismatches
	else:
		return mismatches

def best_match(SEQ1, LIST, MAX = float("inf"), IGNORE_N = 0, PRINT = 0 ):
	"""finds the best match for a sequence in a list of sequences.
	MAX sets the number of max number of mismatches before it moves on.
	Lowering MAX increases performance.
	IGNORE_N = 1 will ignore mismatches with N."""
	x = []
	xcount = []
	y = MAX
	no_exact_match = 1
	#first search for exact matach
	for i in range(len(LIST)):
		if SEQ1 ==  LIST[i]:
			no_exact_match = 0
			return i
			break 
	if no_exact_match:
		for i in range(len(LIST)):
		   z = BC.mismatches(SEQ1, LIST[i], y, IGNORE_N)
		   if z < y:
			  y = z
			  x.append(i)
			  xcount.append(z)
		   if z == 0:
			  break
		if len(x) > 0:
		   comp = "==" + str(min(xcount))
		   best =  [a for a,b in enumerate(xcount) if eval(str(b) + comp)]
		   if PRINT == 1:
			  print SEQ1
			  print LIST[x[0]]
			  print x[best[0]], xcount[best[0]]
		   return x[best[0]]
		else:
		   return -1

#Big support functions
def parse_fastq(directory, f_gzipped_fastqfile, r_gzipped_fastqfile,
		q = "fastq", #the type of fastq file coming off of the sequencer
		f_seqtag_length = 8, #the length of the sequencing tag on the first read (UMI1)
		r_seqtag_length = 8, #the length of the sequencing tag on the second read (UMI2)
		f_multitag_length = 6, #the length of the multiplexing tag on the first read
		r_multitag_length = 6, #the length of the multiplexing tag on the second read
		f_lintag_length = 38, #the length of the lineage tag on the first read (first barcode)
		r_lintag_length = 38,  #the length of the lineage tag on the second read (second barcode)
		f_spacer_length = 43, #distance to first barcode in forward read (ignoring the length the thie multitag and the seqtag)
		r_spacer_length = 29, #distance second barcode in reverse read (ignoring the length the thie multitag and the seqtag)
		min_qs = 30, #the minimum avareage quality score for both lineage tags
		lintag_grep_filter1 ='\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*', #first barcode
		lintag_grep_filter2 ='\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*', #second barcode
		clip_ends = 1, #logical of whether or not to clip the front and back ends off of lintag1 and lintag2
		lintag1_front_clipper = '(.ACC|T.CC|TA.C|TAC.)', #only report lintag1 after this sequence
		lintag2_front_clipper = '(.ACC|T.CC|TA.C|TAC.)', #only report lintag2 after this sequence
		lintag1_rear_clipper = '(.ATA|A.TA|AA.A|AAT.)', #only report lintag1 before this sequence, this must be the COMPLIMENT of the true sequence
		lintag2_rear_clipper = '(.ATA|A.TA|AA.A|AAT.)'):#only report lintag1 before this sequence, this must be the COMPLIMENT of the true sequence
		
	"""
	Parses a F and R gzipped FastQ files and saves the UMIs, multiplexing tags, and barcodes
	Removes reeads where the mean quality score for each lineage tag is not greater than min_qs
	Removes reeads where both lineage tags do not match the regular expression 
	"""

	from Bio import SeqIO
	import os
	import gzip
	import numpy
	import BC
	import re
	from itertools import izip
	os.chdir(directory)
	print("Loading " + f_gzipped_fastqfile + " and " + r_gzipped_fastqfile + " and parsing")
	print( "Saving the combined forward and reverse sequencing tags as seqtag.txt")
	print( "Saving the combined forward and reverse multiplexing tags  as multitag.txt")
	print( "Saving the first lineage tag as lintag1.txt")
	print( "Saving the first lineage tag as lintag2.txt")
	
	
	#assign boundries
	f_boundries = (0, f_seqtag_length , f_multitag_length + f_seqtag_length,
			f_multitag_length + f_seqtag_length + f_spacer_length,
			f_multitag_length + f_seqtag_length + f_spacer_length + f_lintag_length)
	r_boundries = (0, r_seqtag_length , r_multitag_length + r_seqtag_length,
			r_multitag_length + r_seqtag_length + r_spacer_length,
			r_multitag_length + r_seqtag_length + r_spacer_length + r_lintag_length)
	
	
	#open files for writing
	seqtag = open(directory + '_seqtag.txt', 'w')
	multitag = open(directory + '_multitag.txt', 'w')
	lintag1 = open(directory + '_lintag1.txt', 'w')
	lintag2 = open(directory + '_lintag2.txt', 'w')
	
	
	#open files for reading by SeqIO
	f_file = SeqIO.parse(gzip.open(directory + f_gzipped_fastqfile, "rU"), q)
	r_file = SeqIO.parse(gzip.open(directory + r_gzipped_fastqfile, "rU"), q)
	
	#eliminate low quality reads and reads that don't pass a quality filter, optionally clip off ends of lintags
	quality_reads = 0
	total_reads = 0
	for f, r in izip(f_file, r_file):
		fq = f.letter_annotations["phred_quality"]
		rq = r.letter_annotations["phred_quality"]
		total_reads = total_reads + 1
		if numpy.mean(fq[f_boundries[3]:f_boundries[4]]) > min_qs and numpy.mean(rq[r_boundries[3]:r_boundries[4]]) > min_qs: #checks that the quality scores of forward and reverse lintags are OK
			fr = str(f.seq)
			rr = str(r.seq)
			if BC.grep(fr[f_boundries[3]:f_boundries[4]], lintag_grep_filter1) and BC.grep(rr[r_boundries[3]:r_boundries[4]], lintag_grep_filter2):
				#checks the both lineage tags meet the regular expression filter
				quality_reads = quality_reads + 1
				seqtag.write(fr[f_boundries[0]:f_boundries[1]] + rr[r_boundries[0]:r_boundries[1]] +'\n')
				multitag.write(fr[f_boundries[1]:f_boundries[2]] + rr[r_boundries[1]:r_boundries[2]] +'\n')
				ftag = fr[f_boundries[3]:f_boundries[4]]
				rtag = rr[r_boundries[3]:r_boundries[4]] 
				if (clip_ends):
					fstart = re.search(lintag1_front_clipper, ftag).span()[1]
					fend = re.search(lintag1_rear_clipper, ftag[::-1]).span()[1]*-1
					if fend == 0: fend = len(ftag)
					ftag = ftag[fstart:fend]
					rstart = re.search(lintag2_front_clipper, rtag).span()[1]
					rend = re.search(lintag2_rear_clipper, rtag[::-1]).span()[1]*-1
					if rend == 0: rend = len(rtag)
					rtag = rtag[rstart:rend]
				lintag1.write(ftag + '\n')
				lintag2.write(rtag + '\n')
			
	print ( str(quality_reads) + " out of " + str(total_reads) +" reads  passed grep and quality filters")
	seqtag.close()
	multitag.close()
	lintag1.close()
	lintag2.close()
	f_file.close()
	r_file.close()

def parse_fastq_by_multitag(directory, f_gzipped_fastqfile, r_gzipped_fastqfile,
		q = "fastq",
		f_seqtag_length = 8,
		r_seqtag_length = 8,
		f_multitag_length = 6,
		r_multitag_length = 6,
		f_lintag_length = 38,
		r_lintag_length = 38, 
		f_spacer_length = 43, #distance to first barcode in forward read (ignoring the length the thie multitag and the seqtag)
		r_spacer_length = 29, #distance second barcode in reverse read (ignoring the length the thie multitag and the seqtag)
		min_qs = 30, #the minimum avareage quality score for both lineage tags
		lintag_grep_filter1 ='\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*', #first barcode
		lintag_grep_filter2 ='\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*', #second barcode
		clip_ends = 1, #logical of whether or not to clip the front and back ends off of lintag1 and lintag2
		lintag1_front_clipper = '(.ACC|T.CC|TA.C|TAC.)', #only report lintag1 after this sequence
		lintag2_front_clipper = '(.ACC|T.CC|TA.C|TAC.)', #only report lintag2 after this sequence
		lintag1_rear_clipper = '(.TAA|A.AA|AT.A|ATA.)', #only report lintag1 before this sequence, this must be the COMPLIMENT of the true sequence
		lintag2_rear_clipper = '(.TAA|A.AA|AT.A|ATA.)', #only report lintag2 before this sequence, this must be the COMPLIMENT of the true sequence
		multitags = ["TAGCTTGCGTAC", "CGATGTGAGACG"], #concatenated multiplexing tags from the first and second reads that uniquely identify a sample, currently must have 2 or more multitags
		write_multitags = False): #write multitags to file
		
	"""
	Parses a F and R gzipped FastQ files and saves the UMIs, multiplexing tags, and barcodes
	Removes reeads where the mean quality score for each lineage tag is not greater than min_qs
	Removes reeads where both lineage tags do not match the regular expression 
	"""

	from Bio import SeqIO
	import os
	import gzip
	import numpy
	import BC
	import re
	from itertools import izip
	os.chdir(directory)
	print("Loading " + f_gzipped_fastqfile + " and " + r_gzipped_fastqfile + " and parsing")
	print( "Saving the combined forward and reverse sequencing tags as seqtag.txt")
	print( "Saving the combined forward and reverse multiplexing tags  as multitag.txt")
	print( "Saving the first lineage tag as lintag1.txt")
	print( "Saving the first lineage tag as lintag2.txt")
	
	#assign boundries
	f_boundries = (0, f_seqtag_length , f_multitag_length + f_seqtag_length,
			f_multitag_length + f_seqtag_length + f_spacer_length,
			f_multitag_length + f_seqtag_length + f_spacer_length + f_lintag_length)
	r_boundries = (0, r_seqtag_length , r_multitag_length + r_seqtag_length,
			r_multitag_length + r_seqtag_length + r_spacer_length,
			r_multitag_length + r_seqtag_length + r_spacer_length + r_lintag_length)
	
	
	#open files for writing
	#reads that sort to a multiplexing tag
	for i in multitags:
		vars()[i+'_seqtag'] = open(directory + i + '_seqtag.txt', 'w')
		vars()[i+'_lintag1'] = open(directory + i + '_lintag1.txt', 'w')
		vars()[i+'_lintag2'] = open(directory + i + '_lintag2.txt', 'w')
		if write_multitags: vars()[i+'_multitag'] = open(directory + i + '_multitag.txt', 'w')
	
	#reads that do not sort to a multiplexing tag
	unmatched_seqtag = open(directory + 'unmatched_seqtag.txt', 'w')
	unmatched_lintag1 = open(directory + 'unmatched_lintag1.txt', 'w')
	unmatched_lintag2 = open(directory + 'unmatched_lintag2.txt', 'w')
	unmatched_multitag = open(directory + 'unmatched_multitag.txt', 'w')
	
	#reads that do not grep to a barcode
	#grepfail_seqtag = open(directory + 'grepfail_seqtag.txt', 'w')
	#grepfail_lintag1 = open(directory + 'grepfail_lintag1.txt', 'w')
	#grepfail_lintag2 = open(directory + 'grepfail_lintag2.txt', 'w')
	#grepfail_multitag = open(directory + 'grepfail_multitag.txt', 'w')
	#grepfail_all_1= open(directory + 'grepfail_all_1.txt', 'w')
	#grepfail_all_2= open(directory + 'grepfail_all_2.txt', 'w') 
	
	#reads that do not pass quality filter, only output the multitag
	#qualityfail_multitag= open(directory + 'qualityfail_multitag.txt', 'w')
	
	#open files for reading by SeqIO
	f_file = SeqIO.parse(gzip.open(directory + f_gzipped_fastqfile, "rU"), q)
	r_file = SeqIO.parse(gzip.open(directory + r_gzipped_fastqfile, "rU"), q)
	
	#eliminate low quality reads and reads that don't pass a quality filter, optionally clip off ends of lintags
	# sort by multiplexing tags
	quality_reads = 0
	total_reads = 0
	grep_matching_quality_reads = 0
	passing_reads_that_dont_match_a_multitag = 0
	grep_failures = 0
	
	for f, r in izip(f_file, r_file):
		fq = f.letter_annotations["phred_quality"]
		rq = r.letter_annotations["phred_quality"]
		total_reads = total_reads + 1
		if numpy.mean(fq[f_boundries[3]:f_boundries[4]]) > min_qs and numpy.mean(rq[r_boundries[3]:r_boundries[4]]) > min_qs:
			#checks that the quality scores of forward and reverse lintags are OK
			#print "quality ok"
			fr = str(f.seq)
			#print fr
			rr = str(r.seq)
			#print rr
			quality_reads = quality_reads + 1 #these are reads where both lintags pass the quality filters
			if BC.grep(fr[f_boundries[3]:f_boundries[4]], lintag_grep_filter1) and BC.grep(rr[r_boundries[3]:r_boundries[4]], lintag_grep_filter2):
				#checks the both lineage tags meet the regular expression filter
				#print "grep ok"
				grep_matching_quality_reads = grep_matching_quality_reads + 1 #these are reads where both lintags pass the quality and grep filters
				#next, find the closest matching multitag
				m = fr[f_boundries[1]:f_boundries[2]] + rr[r_boundries[1]:r_boundries[2]] #the concatintated multiplexing tag
				#print m
				j = BC.best_match(m, multitags, MAX = (f_multitag_length + r_multitag_length + 1)/3) #best matched multiplexing tag
				#print j
				if j > -1:
					tm = BC.mismatches(m, multitags[j]) #distance to this tag
				else:
					tm = 1000
				if tm < (f_multitag_length + r_multitag_length + 1)/4: #A multitag match has been found
				#if tm < 1: #A multitag match has been found
					ftag = fr[f_boundries[3]:f_boundries[4]]
					rtag = rr[r_boundries[3]:r_boundries[4]] 
					if(clip_ends):
						fstart = re.search(lintag1_front_clipper, ftag).span()[1]
						fend = re.search(lintag1_rear_clipper, ftag[::-1]).span()[1]*-1
						if fend == 0: fend = len(ftag)
						ftag = ftag[fstart:fend]
						rstart = re.search(lintag2_front_clipper, rtag).span()[1]
						rend = re.search(lintag2_rear_clipper, rtag[::-1]).span()[1]*-1
						if rend == 0: rend = len(rtag)
						rtag = rtag[rstart:rend]
					vars()[multitags[j]+'_lintag1'].write(ftag + '\n')
					vars()[multitags[j]+'_lintag2'].write(rtag + '\n')
					vars()[multitags[j]+'_seqtag'].write(fr[f_boundries[0]:f_boundries[1]] + rr[r_boundries[0]:r_boundries[1]] + '\n')
					if write_multitags: vars()[multitags[j]+'_multitag'].write(fr[f_boundries[1]:f_boundries[2]] + rr[r_boundries[1]:r_boundries[2]] + '\n')
					#if (len(fr[f_boundries[1]:f_boundries[2]] + rr[r_boundries[1]:r_boundries[2]]) < 12
					#or len(fr[f_boundries[0]:f_boundries[1]] + rr[r_boundries[0]:r_boundries[1]]) < 16
					#or len(ftag) < 20
					#or len(rtag) < 20): 
					#    print rea
					#    print "match to " + multitags[j]
					#    print "multitag = " + fr[f_boundries[1]:f_boundries[2]] + rr[r_boundries[1]:r_boundries[2]] + " " + str(len(fr[f_boundries[1]:f_boundries[2]] + rr[r_boundries[1]:r_boundries[2]]))
					#    print "seqtag = " + fr[f_boundries[0]:f_boundries[1]] + rr[r_boundries[0]:r_boundries[1]]  + " " + str(len(fr[f_boundries[0]:f_boundries[1]] + rr[r_boundries[0]:r_boundries[1]]))
					#    print "lintag1 = " + ftag + " " + str(len(ftag))
					#    print "lintag2 = " + rtag + " " + str(len(rtag))
					#    break
				else:
					passing_reads_that_dont_match_a_multitag = passing_reads_that_dont_match_a_multitag + 1
					ftag = fr[f_boundries[3]:f_boundries[4]]
					rtag = rr[r_boundries[3]:r_boundries[4]] 
					if(clip_ends):
						fstart = re.search(lintag1_front_clipper, ftag).span()[1]
						fend = re.search(lintag1_rear_clipper, ftag[::-1]).span()[1]*-1
						if fend == 0: fend = len(ftag)
						ftag = ftag[fstart:fend]
						rstart = re.search(lintag2_front_clipper, rtag).span()[1]
						rend = re.search(lintag2_rear_clipper, rtag[::-1]).span()[1]*-1
						if rend == 0: rend = len(rtag)
						rtag = rtag[rstart:rend]
					unmatched_lintag1.write(ftag + '\n')
					unmatched_lintag2.write(rtag + '\n')
					unmatched_seqtag.write(fr[f_boundries[0]:f_boundries[1]] + rr[r_boundries[0]:r_boundries[1]] +'\n')
					unmatched_multitag.write(fr[f_boundries[1]:f_boundries[2]] + rr[r_boundries[1]:r_boundries[2]] +'\n')
						
			else:
				grep_failures = grep_failures + 1
			
				#ftag = fr[f_boundries[3]:f_boundries[4]]
				#rtag = rr[r_boundries[3]:r_boundries[4]] 
				#grepfail_lintag1.write(ftag + '\n')
				#grepfail_lintag2.write(rtag + '\n')
				#grepfail_all_1.write(fr + '\n')
				#grepfail_all_2.write(rr + '\n')
				#grepfail_seqtag.write(fr[f_boundries[0]:f_boundries[1]] + rr[r_boundries[0]:r_boundries[1]] +'\n')
				#grepfail_multitag.write(fr[f_boundries[1]:f_boundries[2]] + rr[r_boundries[1]:r_boundries[2]] +'\n')
		#else:
			#fr = str(f.seq)
			#print fr
			#rr = str(r.seq)
			#qualityfail_multitag.write(fr[f_boundries[1]:f_boundries[2]] + rr[r_boundries[1]:r_boundries[2]] +'\n')		
	print ( str(quality_reads) + " out of " + str(total_reads) +" reads passed quality filters")
	print ( str(grep_failures) + " out of " + str(total_reads) +" did not match the barcode pattern")
	print ( str(grep_matching_quality_reads) + " out of " + str(total_reads) +" reads passed grep and quality filters")
	print ( str(passing_reads_that_dont_match_a_multitag) + " out of " + str(total_reads) +" reads passed grep and quality filters but did not match a multitag")
	for i in multitags:
		vars()[str(i)+'_seqtag'].close() 
		vars()[str(i)+'_lintag1'].close() 
		vars()[str(i)+'_lintag2'].close()
		if write_multitags: vars()[str(i)+'_multitag'].close()
	
	unmatched_seqtag.close()
	unmatched_lintag1.close()
	unmatched_lintag2.close()
	unmatched_multitag.close()
	#grepfail_seqtag.close()
	#grepfail_lintag1.close()
	#grepfail_lintag2.close()
	#grepfail_all_1.close()
	#grepfail_all_2.close()
	#grepfail_multitag.close()
	#qualityfail_multitag.close()
	f_file.close()
	r_file.close()
	
