#!/usr/bin/python

'''
Yue July Xing
06/27/2018

'''

import random
import os
import subprocess
import math
import sys
import time

def log_print(message):
	print '[SimulateCNVs ' + time.asctime( time.localtime(time.time()) ) + '] ' + message
	sys.stdout.flush()

# Function for read in fasta file
def read_fasta(fname):
	with open(fname, "r") as fh:
		name = None
		seqs = {}
		seqs2 = {}
		for line in fh.readlines():
			if line[0]=='>':
				line = line.rstrip()
				name = line[1:]
				name = name.split()[0]
				name = name.split("_")[0]
				seqs[name] = []
			else:
				seqs[name].append(line.rstrip())
	chrs = []
	for key in seqs:
		chrs.append(key)
		seqs2[key] = ''.join(seqs[key])
	fh.close()
	return(seqs2, chrs)

def read_cnv(cnvname, chrs):
	# this file must be pre-sorted!
	cnv_list_st = {}
	cnv_list_ed = {}
	cn = {}
	for ch in chrs:
		cnv_list_st[ch] = []
		cnv_list_ed[ch] = []
		cn[ch] = []
	with open(cnvname, "r") as cnh:
		next(cnh)
		for line in cnh:
			line = line.rstrip()
			line = line.split()
			if any(line[0] in s for s in chrs):
				cnv_list_st[line[0]].append(int(line[1])-1)
				cnv_list_ed[line[0]].append(int(line[2])-1)
				cn[line[0]].append(int(line[4]))
	cnh.close()
	return(cnv_list_st, cnv_list_ed, cn)

def read_cnv_len(lname, rcl_chrs):
	listl = {}
	for ch in rcl_chrs:
		listl[ch] = [] 
	with open(lname, "r") as ln:
		for line in ln:
			line = line.rstrip()
			line = line.split()
			for i in range(int(line[2])):
				listl[line[0]].append(line[1])
	ln.close()
	return(listl)

def make_num_cnv_list(num_cnv, tol_cnv, cnv_len_file, chrs, seqs):
	num_cnv_list = {}
	cnv_listl = None
	m_tol = 0
	if num_cnv:
		m_tol = num_cnv*len(chrs)  
		for ch in chrs:
			num_cnv_list[ch] = num_cnv
	elif tol_cnv:
		tolbp = 0
		for ch in chrs:
			tolbp += len(seqs[ch])		
		for ch in chrs:
			num_cnv_list[ch] = int(math.ceil(len(seqs[ch]) / tolbp * tol_cnv))
			if num_cnv_list[ch] == 0:
				num_cnv_list[ch] += 1
			m_tol += num_cnv_list[ch]
	elif cnv_len_file:
		log_print('Reading CNV length file...')
		if not os.path.exists(cnv_len_file):
			log_print('Error: The CNV length file does not exist!')
			exit(1)
		cnv_listl = read_cnv_len(cnv_len_file, chrs)
		for ch in chrs:
			num_cnv_list[ch] = len(cnv_listl[ch])
			m_tol += num_cnv_list[ch]
	return num_cnv_list, m_tol, cnv_listl

def intersect(a, b):
	return list(set(a) & set(b))

def find_missing(m_opt, m_chrs, m_seqs):
	ms = {}
	for ch in m_chrs:
		ms[ch] = []
		if m_opt:
			for pos, char in enumerate(m_seqs[ch]):
				if char == 'N':
					ms[ch].append(pos)
	return ms

def find_gauss(g_lg=None, g_cnv_min_len=None, g_cnv_max_len=None):
	s = -10
	while s<-4 or s>4:
		s = random.gauss(0,1)
	if g_lg:
		fg = int(0 + (s+4)/8 * (g_lg-1))
	elif g_cnv_min_len and g_cnv_max_len:
		fg = int((s+4)/8 * (g_cnv_max_len-g_cnv_min_len) + g_cnv_min_len)
	return fg

def assign_copy_numbers(chrs, tl, p_ins, min_cn, max_cn, cnv_list_st):
	copy_num = []
	cn = {}
	for ch in chrs:
		cn[ch] = []
	num_ins = int(tl * p_ins)
	num_del = tl - num_ins
	for i in range(num_del):
		copy_num.append(0)
	for i in range(num_ins):
		copy_num.append(random.randrange(min_cn,max_cn+1))
	random.shuffle(copy_num)
	j = 0
	for ch in chrs:
		for i in range(len(cnv_list_st[ch])):
			cn[ch].append(copy_num[j])
			j += 1
	return cn

def write_cnv(chrs, cnv_list_file, cnv_list_st, cnv_list_ed, cn):
	with open(cnv_list_file, 'w') as f:
		line = 'chr\tstart\tend\tlength\tcopy_number\n'
		f.write(line)
		for ch in chrs:
			for i in range(len(cnv_list_st[ch])):
				start = cnv_list_st[ch][i] + 1
				end = cnv_list_ed[ch][i] + 1
				length_cnv = end - start + 1
				line = ch + '\t' + str(start) + '\t' + str(end) + '\t' + \
				str(length_cnv) + '\t' + str(cn[ch][i]) + '\n'
				f.write(line)
	f.close()
	return

def write_cnv_genome(cnv_genome_file, chrs, seqs):
	n = 50
	with open(cnv_genome_file, 'w') as f:
		for ch in chrs:
			header = ">" + ch
			f.write(header + "\n")
			for i in range(0, len(seqs[ch]), n):
				line = seqs[ch][i:i+n]
				f.write(line + "\n")
	f.close()

def call_art(genome, cover, output, read_length, frag_size, stdev, paired, ql, qu):
	#subprocess.call(['samtools','faidx',genome], stderr=None)
	dirn = os.getcwd()
	os.chdir(os.path.dirname(os.path.realpath(__file__)))
	subprocess.call(['chmod', 'u=rwx', 'art_illumina'], stderr=None)
	subprocess.call(['chmod', 'u=rwx', 'call_art.sh'], stderr=None)
	if paired:
		log_print("Paired-end sequencing.")
		pp = "p"
	else:
		log_print("Single-end sequencing.")
		pp = "s"
	os.system('./call_art.sh ' + genome + ' ' + str(read_length) + ' ' + str(cover) + \
		' ' + str(frag_size) + ' ' + str(stdev) + ' ' + output + ' ' + pp + ' ' + str(ql) + ' ' + str(qu))
	os.chdir(dirn)

def make_bam(path_to_picard, path_to_GATK, output_name, out_dir, tmp_dir, paired):
	dirn = os.getcwd()
	os.chdir(os.path.dirname(os.path.realpath(__file__)))
	subprocess.call(['chmod', 'u=rwx', 'Create_bam.sh'])
	if paired:
		pp = "p"
	else:
		pp = "s"
	subprocess.call(['./Create_bam.sh', path_to_picard, path_to_GATK, output_name, out_dir, tmp_dir, pp])
	os.chdir(dirn)
	#os.system('./Create_bam.sh ' + path_to_picard + ' ' + path_to_GATK + ' ' \
	#	+ output_name + ' ' + out_dir + ' ' + tmp_dir)