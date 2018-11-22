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

from Common import *

def g_assign_cnv_pos(chrs, num_cnv_list, cnv_min_len, cnv_max_len, \
	seqs, method_s, method_l, cnv_listl, ran_m, flank):
	c_len = None
	cnv_list_st = {}
	cnv_list_ed = {}
	tol_cnv = 0
	for ch in chrs:
		cnv_list_st[ch] = []
		cnv_list_ed[ch] = []
	for ch in chrs:
		if num_cnv_list[ch] == 0:
			mes = "Chromosome " + str(ch) + " has 0 CNVs."
			log_print(mes)
			continue
		iter_n = num_cnv_list[ch] * 100
		count = 0
		j = 0
		lg = len(seqs[ch])
		while j < iter_n:
			if method_s == 'random':
				cnv_st = random.randint(0,lg-1)
			elif method_s == 'uniform':
				cnv_st = int(random.uniform(0,lg-1))
			else:
				cnv_st = find_gauss(lg)

			if method_l == 'random':
				cnv_ed = cnv_st + random.randint(cnv_min_len,cnv_max_len) - 1
			elif method_l == 'uniform':
				cnv_ed = cnv_st + int(random.uniform(cnv_min_len,cnv_max_len)) - 1
			elif method_l == 'gauss':
				cnv_ed = cnv_st + find_gauss(None,cnv_min_len,cnv_max_len) - 1
			else:
				c_len = random.choice(cnv_listl[ch])
				cnv_ed = cnv_st + int(c_len) - 1
			
			if cnv_ed > (lg-1):
				j += 1
				continue
			cnv_list = range(cnv_st, cnv_ed+1)
			if intersect(cnv_list, ran_m[ch]):
				j += 1
				continue
			flag = 0
			for i in range(len(cnv_list_st[ch])):
				s2 = intersect(cnv_list,range(cnv_list_st[ch][i]-flank, \
					cnv_list_ed[ch][i]+flank+1))
				if s2:
					flag = 1
			if flag == 1:
				j += 1
				continue
			cnv_list_st[ch].append(cnv_st)
			cnv_list_ed[ch].append(cnv_ed) 
			count = count + 1
			mes = "Chromosome " + str(ch) + ": CNV " + str(count)
			log_print(mes)
			if c_len:
				cnv_listl[ch].remove(c_len)
				#c_len = None
			if count == num_cnv_list[ch]:
				break
		if j == iter_n:
			mes = "Chromosome " + str(ch) + " is too small or there are too many CNVs to be added into it."
			log_print(mes)
			mes = "There will be fewer CNVs on chromosome " + str(ch) + ": " + str(count) + " instead of " + str(num_cnv_list[ch]) + "."
			log_print(mes)
		cnv_list_st[ch].sort()
		cnv_list_ed[ch].sort()
		mes = "Generated " + str(count) + " CNV(s) on chromosome " + str(ch) + "."
		log_print(mes)
		tol_cnv = tol_cnv + count		
	mes = "Total CNV(s) generated: " + str(tol_cnv)
	log_print(mes)
	return(cnv_list_st, cnv_list_ed)

# Function to generate rearranged genome
def g_gen_rearranged_genome(chrs, n_cnv_list_st, n_cnv_list_ed, cn, n_seqs):
	cnv_list_st = dict(n_cnv_list_st)
	cnv_list_ed = dict(n_cnv_list_ed)
	seqs = dict(n_seqs)
	for ch in chrs:
		for i in range(len(cn[ch])):
			cnv_st = cnv_list_st[ch][i]
			cnv_ed = cnv_list_ed[ch][i]
			length = cnv_ed - cnv_st + 1
			pro_cnv_list_st = cnv_list_st[ch][(i+1):]
			pro_cnv_list_ed = cnv_list_ed[ch][(i+1):]
			if cn[ch][i] == 0:
				seqs[ch] = seqs[ch][:cnv_st] + seqs[ch][(cnv_ed+1):]
				for k in range(len(pro_cnv_list_st)):
					pro_cnv_list_st[k] -= length
					pro_cnv_list_ed[k] -= length 
			elif cn[ch][i] > 0:
				seqs[ch] = seqs[ch][:cnv_st] + seqs[ch][cnv_st:(cnv_ed+1)]*cn[ch][i] \
				+ seqs[ch][(cnv_ed+1):]     
				for k in range(len(pro_cnv_list_st)):
					pro_cnv_list_st[k] += length * (cn[ch][i] - 1)
					pro_cnv_list_ed[k] += length * (cn[ch][i] - 1) 
			cnv_list_st[ch] = cnv_list_st[ch][:(i+1)] + pro_cnv_list_st
			cnv_list_ed[ch] = cnv_list_ed[ch][:(i+1)] + pro_cnv_list_ed
	return seqs


# Simulation
def simulate_WGS(sim_params, gin_seqs, gin_chrs, sim_control):
	in_ran_m = None
	in_seqs = dict(gin_seqs)
	in_chrs = list(gin_chrs)
	ori_seqs = dict(gin_seqs)
	in_genome_file = sim_params['genome_file']
	in_cnvname = sim_params['g_cnv_list']
	in_num_cnv = sim_params['g_cnv_chr']
	in_tol_cnv = sim_params['g_cnv_tol']
	in_cnv_min_len = sim_params['cnv_min_len'] 
	in_cnv_max_len = sim_params['cnv_max_len'] 
	in_p_ins = sim_params['p_ins'] 
	in_min_cn = sim_params['min_cn']
	in_max_cn = sim_params['max_cn']
	in_cnv_list_file = os.path.join(sim_params['out_dir'], sim_params['rearranged_out']+".cnv.bed")
	rearranged_out_name = os.path.join(sim_params['out_dir'], sim_params['rearranged_out'])
	control_out_name = os.path.join(sim_params['out_dir'], 'control')
	out_cnv_genome_file = rearranged_out_name + '.fa'
	in_cover = sim_params['coverage']
	in_read_length = sim_params['read_length']
	in_frag_size = sim_params['frag_size']
	in_stdev = sim_params['stdev']
	in_paired_end = sim_params['paired_end']
	in_ql = sim_params['ql']
	in_qu = sim_params['qu']
	#in_sim_control = sim_params['sim_control']
	in_sim_control = sim_control
	in_sim_short_reads = sim_params['sim_short_reads']
	in_sim_bam = sim_params['sim_bam']
	in_method_s = sim_params['method_s']
	in_method_l = sim_params['method_l']
	in_cnv_len_file = sim_params['g_cnv_len_file']
	in_flank = sim_params['flank']
	opt = sim_params['opt']

	'''
	log_print('Reading genome file...')
	in_seqs, in_chrs = g_read_fasta(in_genome_file)
	'''
	
	# Generate CNVs that are randomly distributed in the genome
		# If #CNVs given for the whole genome, #CNVs on each chromosome is 
		# proportional to the length of the chromosome 
		# If #CNVs given for each chromosome, #CNVs on each chromosome is 
		# the same
	
	# Generate CNVs overlapping with exons and not overlapping with exons
	# for CNVs overlapping with target regions
	if in_cnvname:
		log_print('Reading CNVs from provided file...')
		if not os.path.exists(in_cnvname):
			log_print('Error: The provided CNV list does not exist!')
			exit(1)
		else:
			in_cnv_list_st, in_cnv_list_ed, in_cn = read_cnv(in_cnvname, in_chrs)			
	else:
		log_print('Exclude missing sequences in the genome...')
		in_ran_m = find_missing(opt, in_chrs, in_seqs)
		log_print('Generating CNVs...')
		in_num_cnv_list, tol, in_cnv_listl = make_num_cnv_list(in_num_cnv, \
			in_tol_cnv, in_cnv_len_file, in_chrs, in_seqs)
		in_cnv_list_st, in_cnv_list_ed = g_assign_cnv_pos(in_chrs, in_num_cnv_list, \
			in_cnv_min_len, in_cnv_max_len, in_seqs, in_method_s, in_method_l, \
			in_cnv_listl, in_ran_m, in_flank)
		in_cn = assign_copy_numbers(in_chrs, tol, in_p_ins, in_min_cn, in_max_cn, \
			in_cnv_list_st)

	# write CNV lists into files
	if not in_cnvname:
		log_print('Writing CNVs to file...')
		write_cnv(in_chrs, in_cnv_list_file, in_cnv_list_st, in_cnv_list_ed, in_cn)
	
	# Generate rearranged genome
	log_print('Generating rearranged genome...')
	seqs_new = g_gen_rearranged_genome(in_chrs, in_cnv_list_st, in_cnv_list_ed, \
		in_cn, in_seqs)


	# Write rearranged genome and targets for Wessim
	log_print('Writing rearranged genome to file...')
	write_cnv_genome(out_cnv_genome_file, in_chrs, seqs_new)    
	
	if in_sim_control:
		subprocess.call(['cp', in_genome_file, os.path.join(sim_params['out_dir'],'control.fa')])
		#shutil.copy2(in_genome_file, os.path.join(sim_params['out_dir'],'control.fa'))
		#write_cnv_genome(os.path.join(sim_params['out_dir'],'control.fa'), in_chrs, ori_seqs)

	
	# Simulation with ART_illumina
	if in_sim_short_reads:
		log_print('Simulation with ART_illumina for rearranged genome...')
		call_art(out_cnv_genome_file, in_cover, rearranged_out_name, \
			in_read_length, in_frag_size, in_stdev, in_paired_end, in_ql, in_qu)
		if in_sim_control:
			log_print('Simulation with ART_illumina for control genome...')
			call_art(in_genome_file, in_cover, control_out_name, \
				in_read_length, in_frag_size, in_stdev, in_paired_end, in_ql, in_qu)

	# Simulate bam files
	if in_sim_bam:
		ref_genome = os.path.join(sim_params['out_dir'], 'control.fa')
		if not os.path.exists(ref_genome):
			subprocess.call(['cp', in_genome_file, ref_genome])
			#shutil.copy2(in_genome_file, ref_genome)
			#write_cnv_genome(ref_genome, in_chrs, ori_seqs)

		if os.path.exists(os.path.join(sim_params['out_dir'], 'control.dict')) and \
		os.path.exists(os.path.join(sim_params['out_dir'], 'control.fa.sa')) and \
		os.path.exists(os.path.join(sim_params['out_dir'], 'control.fa.fai')):
			log_print('Index files already exist. Skip creating index files...')
		else:
			log_print('Create index files for the reference...')
			if not os.path.exists(os.path.join(sim_params['out_dir'], 'control.dict')):
				picard_path = sim_params['path_to_picard'] + "/picard.jar"
				ref_genome_dict = os.path.join(sim_params['out_dir'], 'control.dict')
				subprocess.call(['java', '-jar', picard_path, 'CreateSequenceDictionary', \
					'REFERENCE=' + ref_genome, 'OUTPUT=' + ref_genome_dict])
			if not os.path.exists(os.path.join(sim_params['out_dir'], 'control.fa.sa')):
				subprocess.call(['bwa', 'index', ref_genome])
			if not os.path.exists(os.path.join(sim_params['out_dir'], 'control.fa.fai')):
				subprocess.call(['samtools','faidx',ref_genome])
			if not os.path.exists(os.path.join(sim_params['out_dir'], 'control.dict')):
				log_print('Error: Fail to create index (.dict) for the control.')
				exit(1)
			if not os.path.exists(os.path.join(sim_params['out_dir'], 'control.fa.fai')):
				log_print('Error: Fail to create index (.fai) for the control.')
				exit(1)
			if not os.path.exists(os.path.join(sim_params['out_dir'], 'control.fa.sa')):
				log_print('Error: Fail to create bwa indexes for the control.')
				exit(1)

		log_print('Simulating bam file for rearranged genome...')
		make_bam(sim_params['path_to_picard'], sim_params['path_to_GATK'], sim_params['rearranged_out'], \
			sim_params['out_dir'], sim_params['tmp_dir'], in_paired_end)
		if in_sim_control:
			log_print('Simulating bam file for control genome...')
			make_bam(sim_params['path_to_picard'], sim_params['path_to_GATK'], 'control', \
				sim_params['out_dir'], sim_params['tmp_dir'], in_paired_end)