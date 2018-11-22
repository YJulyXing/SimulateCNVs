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

# Function for read in target region file
def read_target(tname, chrs):
	st = {}
	ed = {}
	for ch in chrs:
		st[ch] = []
		ed[ch] = []
	with open(tname, "r") as th:
		for line in th:
			line = line.rstrip()
			line = line.split()
			if any(line[0] in s for s in chrs):
				st[line[0]].append(int(line[1])-1)
				ed[line[0]].append(int(line[2])-1)
	# Sort by chrs
	for ch in chrs:
		st[ch].sort()
		ed[ch].sort()
	th.close()
	return(st, ed)

def assign_cnv_pos(chrs, st, ed, num_cnv_list, cnv_min_len, cnv_max_len, \
	overlap_bp, seqs, method_s, method_l, cnv_listl, ran_m, flank):
	c_len = None
	cnv_list_st = {}
	cnv_list_ed = {}
	tol_cnv = 0
	for ch in chrs:
		cnv_list_st[ch] = []
		cnv_list_ed[ch] = []
	for ch in chrs:
		if len(st[ch]) == 0:
			mes =  "Chromosome " + str(ch) + " has no target regions. No CNVs will be generated on it."
			log_print(mes)
			continue
		if num_cnv_list[ch] == 0:
			mes = "Chromosome " + str(ch) + " has 0 CNVs."
			log_print(mes)
			continue
		iter_n = num_cnv_list[ch] * 100
		count = 0
		j = 0
		lg = len(seqs[ch])

		target = []
		for i in range(len(st[ch])):
			target += range(st[ch][i],ed[ch][i]+1)
		
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

			if len(intersect(target, cnv_list)) < overlap_bp:
				j += 1
				continue
			
			flag2 = 0
			for i in range(len(cnv_list_st[ch])):
				s2 = intersect(cnv_list,range(cnv_list_st[ch][i]-flank, \
					cnv_list_ed[ch][i]+flank+1))
				if s2:
					flag2 = 1
			if flag2 == 1:
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
			mes = "Chromosome " + str(ch) + " is too small or there are too many CNVs to be generated or -ol is too large."
			log_print(mes)
			mes = "There will be fewer CNVs on chromosome " + str(ch) + ": " + str(count) + " instead of " + str(num_cnv_list[ch]) + "."
			log_print(mes)
		cnv_list_st[ch].sort()
		cnv_list_ed[ch].sort()
		mes = "Generated " + str(count) + " CNV(s) on chromosome " + str(ch) + "."
		log_print(mes)
		tol_cnv = tol_cnv + count		
	mes = "Total CNV(s) overlapping with exons generated: " + str(tol_cnv)
	log_print(mes)
	return(cnv_list_st, cnv_list_ed)

def assign_out_cnv_pos(chrs, st, ed, num_cnv_list, cnv_min_len, cnv_max_len, \
	seqs, cnv_ex_list_st, cnv_ex_list_ed, method_s, method_l, cnv_listl, ran_m, flank):
	cnv_list_st = {}
	cnv_list_ed = {}
	c_len = None
	tol_cnv = 0
	for ch in chrs:
		cnv_list_st[ch] = []
		cnv_list_ed[ch] = []
	for ch in chrs:
		if len(st[ch]) == 0:
			mes =  "Chromosome " + str(ch) + " has no target regions. No CNVs will be generated on it."
			log_print(mes)
			continue
		if num_cnv_list[ch] == 0:
			mes = "Chromosome " + str(ch) + " has 0 CNVs."
			log_print(mes)
			continue
		iter_n = num_cnv_list[ch] * 100
		count = 0
		j = 0
		lg = len(seqs[ch])

		ran_cnv = []
		for t in range(len(cnv_ex_list_st[ch])):
			ran_cnv += range(cnv_ex_list_st[ch][t]-flank, cnv_ex_list_ed[ch][t]+flank+1)

		target = []
		for i in range(len(st[ch])):
			target += range(st[ch][i],ed[ch][i]+1)

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
			if intersect(cnv_list, ran_cnv):
				j += 1
				continue
			if intersect(cnv_list, target):
				j += 1
				continue

			flag2 = 0
			for i in range(len(cnv_list_st[ch])):
				s2 = intersect(cnv_list, range(cnv_list_st[ch][i]-flank, \
					cnv_list_ed[ch][i]+flank+1))
				if s2:
					flag2 = 1
			if flag2 == 1:
				j += 1
				continue
			flag3 = 0
			for i in range(len(cnv_ex_list_st[ch])):
				s3 = intersect(cnv_list, range(cnv_ex_list_st[ch][i]-flank, \
					cnv_ex_list_ed[ch][i]+flank+1))
				if s3:
					flag3 = 1
			if flag3 == 1:
				j += 1
				continue
			cnv_list_st[ch].append(cnv_st)
			cnv_list_ed[ch].append(cnv_ed)
			count = count + 1
			mes = "Chromosome " + str(ch) + ": CNV " + str(count)
			log_print(mes)
			if c_len:
				cnv_listl[ch].remove(c_len)
			if count == num_cnv_list[ch]:
				break
		if j == iter_n:
			mes = "Chromosome " + str(ch) + " is too small or there are too many CNVs to be generated on it."
			log_print(mes)
			mes = "There will be fewer CNVs on chromosome " + str(ch) + ":" + str(count) + " instead of " + str(num_cnv_list[ch]) + "."
			log_print(mes)
		cnv_list_st[ch].sort()
		cnv_list_ed[ch].sort()
		mes = "Generated " + str(count) + " CNV(s) on chromosome " + str(ch) + "."
		log_print(mes)
		tol_cnv = tol_cnv + count
	mes = "Total CNV(s) outside of exons generated: " + str(tol_cnv)
	log_print(mes)
	return(cnv_list_st, cnv_list_ed)

# Function to generate rearranged genome
def gen_rearranged_genome(chrs, n_cnv_list_st, n_cnv_list_ed, cn, n_st, n_ed, n_seqs):
	cnv_list_st = dict(n_cnv_list_st)
	cnv_list_ed = dict(n_cnv_list_ed)
	st = dict(n_st) 
	ed = dict(n_ed) 
	seqs = dict(n_seqs)
	for ch in chrs:
		for i in range(len(cn[ch])):
			cnv_st = cnv_list_st[ch][i]
			cnv_ed = cnv_list_ed[ch][i]
			length = cnv_ed - cnv_st + 1
			#pre_cnv_list_st = cnv_list_st[ch][:i]
			#pre_cnv_list_ed = cnv_list_ed[ch][:i]
			pro_cnv_list_st = cnv_list_st[ch][(i+1):]
			pro_cnv_list_ed = cnv_list_ed[ch][(i+1):]
			pre_cnv_st = []
			pre_cnv_ed = []
			in_cnv_st = []
			in_cnv_ed = []
			pro_cnv_st = []
			pro_cnv_ed = []
			for j in range(len((st[ch]))):
				if st[ch][j] < cnv_st and ed[ch][j] >= cnv_st and ed[ch][j] <= cnv_ed:
					pre_cnv_st.append(st[ch][j])
					pre_cnv_ed.append(cnv_st-1)
					in_cnv_st.append(cnv_st)
					in_cnv_ed.append(ed[ch][j])
				elif st[ch][j] < cnv_st and ed[ch][j] > cnv_ed:
					pre_cnv_st.append(st[ch][j])
					pre_cnv_ed.append(cnv_st-1)
					in_cnv_st.append(cnv_st)
					in_cnv_ed.append(cnv_ed)
					pro_cnv_st.append(cnv_ed+1)
					pro_cnv_ed.append(ed[ch][j])
				elif st[ch][j] >= cnv_st and ed[ch][j] <= cnv_ed:
					in_cnv_st.append(st[ch][j])
					in_cnv_ed.append(ed[ch][j])
				elif st[ch][j] <= cnv_ed and ed[ch][j] > cnv_ed:
					in_cnv_st.append(st[ch][j])
					in_cnv_ed.append(cnv_ed)
					pro_cnv_st.append(cnv_ed+1)
					pro_cnv_ed.append(ed[ch][j])
				elif ed[ch][j] < cnv_st:
					pre_cnv_st.append(st[ch][j])
					pre_cnv_ed.append(ed[ch][j])
				elif st[ch][j] > cnv_ed:
					pro_cnv_st.append(st[ch][j])
					pro_cnv_ed.append(ed[ch][j])
			if cn[ch][i] == 0:
				in_cnv_st = []
				in_cnv_ed = []
				for k in range(len(pro_cnv_st)):
					pro_cnv_st[k] -= length 
					pro_cnv_ed[k] -= length
				seqs[ch] = seqs[ch][:cnv_st] + seqs[ch][(cnv_ed+1):]
				for k in range(len(pro_cnv_list_st)):
					pro_cnv_list_st[k] -= length
					pro_cnv_list_ed[k] -= length 
			elif cn[ch][i] > 0:
				in_cnv_st_new = []
				in_cnv_ed_new = []
				for k in range(cn[ch][i]):
					for s in in_cnv_st:
						in_cnv_st_new.append(length * k + s)
					for s in in_cnv_ed:
						in_cnv_ed_new.append(length * k + s)
				in_cnv_st = in_cnv_st_new
				in_cnv_ed = in_cnv_ed_new
				for k in range(len(pro_cnv_st)):
					pro_cnv_st[k] += length * (cn[ch][i] - 1)
					pro_cnv_ed[k] += length * (cn[ch][i] - 1)
				seqs[ch] = seqs[ch][:cnv_st] + seqs[ch][cnv_st:(cnv_ed+1)]*cn[ch][i] \
				+ seqs[ch][(cnv_ed+1):]           
				for k in range(len(pro_cnv_list_st)):
					pro_cnv_list_st[k] += length * (cn[ch][i] - 1)
					pro_cnv_list_ed[k] += length * (cn[ch][i] - 1) 
			if (len(pre_cnv_ed) != 0) and (len(in_cnv_st) != 0):
				if (pre_cnv_ed[-1] == in_cnv_st[0] - 1):
					del pre_cnv_ed[-1]
					del in_cnv_st[0]
			if (len(in_cnv_ed) != 0) and (len(pro_cnv_st) != 0):
				if (in_cnv_ed[-1] == pro_cnv_st[0] -1):
					del in_cnv_ed[-1]
					del pro_cnv_st[0]
			st[ch] = pre_cnv_st + in_cnv_st + pro_cnv_st
			ed[ch] = pre_cnv_ed + in_cnv_ed + pro_cnv_ed
			cnv_list_st[ch] = cnv_list_st[ch][:(i+1)] + pro_cnv_list_st
			cnv_list_ed[ch] = cnv_list_ed[ch][:(i+1)] + pro_cnv_list_ed
	return st, ed, seqs

def write_targets(targets_file, chrs, w_st, w_ed, inter):
	st = dict(w_st)
	ed = dict(w_ed)
	st_w = {}
	ed_w = {}
	if inter:
		for ch in chrs:
			rag = range(len(st[ch]))
			if rag:
				del rag[-1]
				st_w[ch] = [st[ch][0]]
				ed_w[ch] = []
				for t in rag:
					if (st[ch][t+1] - ed[ch][t] > inter):
						st_w[ch].append(st[ch][t+1])
						ed_w[ch].append(ed[ch][t])
				ed_w[ch].append(ed[ch][-1])
			else:
				st_w[ch] = []
				ed_w[ch] = []
	else:
		st_w = dict(st)
		ed_w = dict(ed)
	with open(targets_file, 'w') as f:
		for ch in chrs:
			for i in range(len(st[ch])):
				start = st[ch][i] + 1
				end = ed[ch][i] + 1
				line = ch + '\t' + str(start) + '\t' + str(end) + '\n'
				f.write(line)
	f.close()
	return(st_w, ed_w)

def write_exon_fatsta(exon_fasta_file, seqs, chrs, st, ed, fl):
	n = 50
	with open(exon_fasta_file, 'w') as f:
		for ch in chrs:
			ln = len(seqs[ch])
			for i in range(len(st[ch])):
				if st[ch][i]-fl >= 0:
					n_st = st[ch][i]-fl
					start = st[ch][i] - fl + 1
				else:
					n_st = 0
					start = 1			
				if ed[ch][i]+fl+1 <= ln:
					seq_i = seqs[ch][n_st:(ed[ch][i]+fl+1)]
					end = ed[ch][i] + fl + 1
				else:
					seq_i = seqs[ch][n_st:]
					end = ln
				header_i = ">" + ch + '_' + str(start) + '_' + str(end) + '\n'
				f.write(header_i)
				for t in range(0, len(seq_i), n):
					line = seq_i[t:t+n]
					f.write(line + "\n")
	f.close()


# Simulation
def simulate_WES(sim_params, ein_seqs, ein_chrs, ein_st, ein_ed, sim_control, eflag):
	in_out_cn = None
	in_ran_m = None
	in_seqs = dict(ein_seqs)
	in_chrs = list(ein_chrs)
	in_st = dict(ein_st)
	in_ed = dict(ein_ed)
	ori_st = dict(ein_st)
	ori_ed = dict(ein_ed)
	ori_seqs = dict(ein_seqs)
	in_genome_file = sim_params['genome_file']
	in_targets_file = sim_params['target_region_file']
	in_cnvname = sim_params['e_cnv_list']
	in_num_cnv = sim_params['e_cnv_chr']
	in_tol_cnv = sim_params['e_cnv_tol']
	in_out_cnvname = sim_params['o_cnv_list']
	in_num_cnv_out = sim_params['o_cnv_chr']
	in_tol_cnv_out = sim_params['o_cnv_tol']
	in_cnv_min_len = sim_params['cnv_min_len'] 
	in_cnv_max_len = sim_params['cnv_max_len'] 
	in_overlap_bp = sim_params['overlap_bp'] 
	in_p_ins = sim_params['p_ins'] 
	in_min_cn = sim_params['min_cn']
	in_max_cn = sim_params['max_cn']
	in_cnv_list_file = os.path.join(sim_params['out_dir'], sim_params['rearranged_out']+".cnv.overlap_exon.bed")
	in_cnv_out_list_file = os.path.join(sim_params['out_dir'], sim_params['rearranged_out']+".cnv.out_of_exon.bed")
	out_cnv_targets_file = os.path.join(sim_params['out_dir'], sim_params['rearranged_out']+".target_regions_for_gen_short_reads.bed")
	out_exon_fasta_file = os.path.join(sim_params['tmp_dir'], sim_params['rearranged_out']+".target_region_fasta.fa")
	out_control_fasta_file = os.path.join(sim_params['tmp_dir'], "control_target_region_fasta.fa")
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
	in_cnv_len_file = sim_params['e_cnv_len_file']
	in_cnv_len_file_out = sim_params['o_cnv_len_file']
	in_flank = sim_params['flank']
	opt = sim_params['opt']
	in_fl = sim_params['fl']
	in_inter = sim_params['inter']

	'''
	log_print('Reading genome file...')
	in_seqs, in_chrs = read_fasta(in_genome_file)
	log_print('Reading target region file...')
	in_st, in_ed = read_target(in_targets_file, in_chrs)
	'''

	# Generate CNVs that are randomly distributed in the genome
		# If #CNVs given for the whole genome, #CNVs on each chromosome is 
		# proportional to the length of the chromosome 
		# If #CNVs given for each chromosome, #CNVs on each chromosome is 
		# the same
	
	# Generate CNVs overlapping with exons and not overlapping with exons
	# for CNVs overlapping with target regions
	if in_cnvname:
		log_print('Reading CNVs overlapping with exons from provided file...')
		if not os.path.exists(in_cnvname):
			log_print('Error: The provided CNV list does not exist!')
			exit(1)
		else:
			in_cnv_list_st, in_cnv_list_ed, in_cn = read_cnv(in_cnvname, in_chrs)			
	else:
		if opt:
			log_print('Exclude missing sequences in the genome...')
		in_ran_m = find_missing(opt, in_chrs, in_seqs)
		log_print('Generating CNVs overlapping with exons...')
		in_num_cnv_list, tol, in_cnv_listl = make_num_cnv_list(in_num_cnv, \
			in_tol_cnv, in_cnv_len_file, in_chrs, in_seqs)
		in_cnv_list_st, in_cnv_list_ed = assign_cnv_pos(in_chrs, in_st, in_ed, in_num_cnv_list, \
			in_cnv_min_len, in_cnv_max_len, in_overlap_bp, in_seqs, in_method_s, in_method_l, \
			in_cnv_listl, in_ran_m, in_flank)
		in_cn = assign_copy_numbers(in_chrs, tol, in_p_ins, in_min_cn, in_max_cn, \
			in_cnv_list_st)

	# for CNVs not overlapping with target regions (optional)
	if in_out_cnvname:
		log_print('Reading CNVs outside of exons from provided file...')
		if not os.path.exists(in_out_cnvname):
			log_print('Error: The provided CNV list does not exist!')
			exit(1)
		else:
			in_out_cnv_list_st, in_out_cnv_list_ed, in_out_cn = read_cnv(in_out_cnvname, in_chrs)	
	elif in_tol_cnv_out or in_num_cnv_out or in_cnv_len_file_out:
		log_print('Generating CNVs outside of exons...')
		if not in_ran_m:
			if opt:
				log_print('Exclude missing sequences in the genome...')
			in_ran_m = find_missing(opt, in_chrs, in_seqs)
		in_num_cnv_out_list, tol_out, in_out_cnv_listl = make_num_cnv_list(in_num_cnv_out, \
			in_tol_cnv_out, in_cnv_len_file_out, in_chrs, in_seqs)
		in_out_cnv_list_st, in_out_cnv_list_ed = assign_out_cnv_pos(in_chrs, in_st, in_ed, in_num_cnv_out_list, \
			in_cnv_min_len, in_cnv_max_len, in_seqs, \
			in_cnv_list_st, in_cnv_list_ed, in_method_s, in_method_l, \
			in_out_cnv_listl, in_ran_m, in_flank)
		in_out_cn = assign_copy_numbers(in_chrs, tol_out, in_p_ins, in_min_cn, in_max_cn, \
			in_out_cnv_list_st)

	# write CNV lists into files
	if not in_cnvname:
		log_print('Writing CNVs overlapping with exons to file...')
		write_cnv(in_chrs, in_cnv_list_file, in_cnv_list_st, in_cnv_list_ed, in_cn)
	if in_out_cn and (not in_out_cnvname):
		log_print('Writing CNVs outside of exons to file...')
		write_cnv(in_chrs, in_cnv_out_list_file, in_out_cnv_list_st, in_out_cnv_list_ed, in_out_cn)
	
	# Generate rearranged genome
	log_print('Generating rearranged genome...')
	if in_out_cn:
		for ch in in_chrs:
			in_cnv_list_st[ch] = in_cnv_list_st[ch] + in_out_cnv_list_st[ch]
			in_cnv_list_ed[ch] = in_cnv_list_ed[ch] + in_out_cnv_list_ed[ch]
			in_cn[ch] = in_cn[ch] + in_out_cn[ch]
			list_all = [[]] * len(in_cn[ch])
			for i in range(len(in_cn[ch])):
				list_all[i] = [in_cnv_list_st[ch][i], in_cnv_list_ed[ch][i], in_cn[ch][i]]
			list_all.sort()

			for i in range(len(in_cn[ch])):
				in_cnv_list_st[ch][i] = list_all[i][0]
				in_cnv_list_ed[ch][i] = list_all[i][1]
				in_cn[ch][i] = list_all[i][2]

	st_new, ed_new, seqs_new = gen_rearranged_genome(in_chrs, in_cnv_list_st, in_cnv_list_ed, \
		in_cn, in_st, in_ed, in_seqs)

	# Write rearranged genome and targets
	log_print('Writing rearranged genome and target regions to file...')
	st_newf, ed_newf = write_targets(out_cnv_targets_file, in_chrs, st_new, ed_new, in_inter)
	write_cnv_genome(out_cnv_genome_file, in_chrs, seqs_new)

	if eflag == 0:
		orif_st, orif_ed = write_targets(os.path.join(sim_params['out_dir'],'target_regions.bed'), 
			in_chrs, ori_st, ori_ed, in_inter)
	
	if in_sim_control:
		subprocess.call(['cp', in_genome_file, os.path.join(sim_params['out_dir'],'control.fa')])
		#write_cnv_genome(os.path.join(sim_params['out_dir'],'control.fa'), in_chrs, ori_seqs)
		#shutil.copy2(in_genome_file, os.path.join(sim_params['out_dir'],'control.fa'))
		

	# Simulation with ART
	if in_sim_short_reads:
		log_print('Simulating short reads for rearranged genome...')
		write_exon_fatsta(out_exon_fasta_file, seqs_new, in_chrs, st_newf, ed_newf, in_fl)
		call_art(out_exon_fasta_file, in_cover, rearranged_out_name, in_read_length, \
			in_frag_size, in_stdev, in_paired_end, in_ql, in_qu)
		if in_sim_control:
			log_print('Simulating short reads for control genome...')
			write_exon_fatsta(out_control_fasta_file, in_seqs, in_chrs, orif_st, orif_ed, in_fl)
			call_art(out_control_fasta_file, in_cover, control_out_name, in_read_length, \
				in_frag_size, in_stdev, in_paired_end, in_ql, in_qu)

	# Simulate bam files
	if in_sim_bam:
		ref_genome = os.path.join(sim_params['out_dir'], 'control.fa')
		if not os.path.exists(ref_genome):
			subprocess.call(['cp', in_genome_file, ref_genome])
			#write_cnv_genome(ref_genome, in_chrs, ori_seqs)
			#shutil.copy2(in_genome_file, ref_genome)

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