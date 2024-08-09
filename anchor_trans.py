import sys,os
anchor = f"{sys.argv[1]}.anchor"
with open(anchor,'r') as f:
	dic_ref = {}
	for x in range(14):
		dic_ref[f"Chr{x+1}"]={}
	for l in f:
		l0=l.strip().split()
		if 'refChr' in l:
			#print(l.strip())
			pass
		else:
			if '#' not in l:
				l0=l.strip().split()
				refChr,referenceStart,referenceEnd,queryChr,queryStart,queryEnd,strand,gene,blockIndex,score = l0[0],l0[1],l0[2],l0[3],l0[4],l0[5],l0[6],l0[7],l0[8],l0[9]
				#if refChr != queryChr:
				if int(referenceEnd) - int(referenceStart) < 0:continue
				pos_key = int(referenceStart)
				dic_ref[refChr][pos_key] = f"{referenceStart}\t{referenceEnd}\t{queryStart}\t{queryEnd}\t{strand}\t{refChr}\t{queryChr}\t{gene}\t{blockIndex}\t{score}"
block = 0
block_dic ={}
for x in dic_ref.keys():
	block_dic[x] = {}
	tmp_region_list_ref,tmp_region_list_alt = [],[]
	for num in range(len(list(dic_ref[x]))):
		if num == len(list(dic_ref[x]))-1:
			list_key_sorted_N0 = sorted(list(dic_ref[x]))[num]
			pos_key_N0 = dic_ref[x][list_key_sorted_N0]
			pos_key_list_N0 = pos_key_N0.split()
			pos_key_strand_N0 = pos_key_list_N0[4]
			ref_chr, alt_chr = pos_key_list_N0[5],pos_key_list_N0[6]
			label_Next = 'Next_Chr_break'  #Break
			tmp_region_list_ref='\t'.join(tmp_region_list_ref)
			tmp_region_list_alt='\t'.join(tmp_region_list_alt)
			block +=1
			if pos_key_strand_N0 == "-":
				inversion = 'Inversion'
			else:
				inversion = 'Normal'
			pri=f"Block{block}\t{x}\t{tmp_region_list_ref}\t{alt_chr}\t{tmp_region_list_alt}\t{pos_key_strand_N0}\t{inversion}"
			print(pri)
			block_dic[x][f"Block{block}"] = pri
			#print(f"Block{block}",x,tmp_region_list_ref,tmp_region_list_alt,label_Next,pos_key_strand_N0)
		#	print(num,x,pos_key_N0,label_Next,tmp_region_list_ref,tmp_region_list_alt)
			tmp_region_list_ref,tmp_region_list_alt = [],[]
		else:
			list_key_sorted_N0,list_key_sorted_N1 = sorted(list(dic_ref[x]))[num],sorted(list(dic_ref[x]))[num+1]
			pos_key_N0, pos_key_N1 = dic_ref[x][list_key_sorted_N0], dic_ref[x][list_key_sorted_N1]

			pos_key_list_N0, pos_key_list_N1 = pos_key_N0.split(),pos_key_N1.split()
			pos_key_strand_N0, pos_key_strand_N1 = pos_key_list_N0[4], pos_key_list_N1[4]
			ref_chr, alt_chr = pos_key_list_N0[5],pos_key_list_N0[6]
			if pos_key_strand_N0 == pos_key_strand_N1: #strand consistent or not
				if int(pos_key_list_N0[1])+1 == int(pos_key_list_N1[0]): #ref consistent or not
					if pos_key_strand_N0 == '+':
						#if int(pos_key_list_N0[3])+1 == int(pos_key_list_N1[2]): #query consistent or not
						label_Next = 'Next_continue'
						tmp_region_list_ref.append(pos_key_list_N0[0])
						tmp_region_list_ref.append(pos_key_list_N1[1])
						tmp_region_list_ref = [tmp_region_list_ref[0],tmp_region_list_ref[-1]]
						tmp_region_list_alt.append(pos_key_list_N0[2])
						tmp_region_list_alt.append(pos_key_list_N1[3])	
						tmp_region_list_alt = [tmp_region_list_alt[0],tmp_region_list_alt[-1]]
					#		print(num,x,pos_key_N0,'Next_continue',tmp_region_list_ref,tmp_region_list_alt)
						#else:
						#	label_Next = 'Warning1:Next_continueButQueryBreak'
						#	print(num,x,pos_key_N0,label_Next,tmp_region_list_ref,tmp_region_list_alt)
						#	tmp_region_list_ref,tmp_region_list_alt = [],[]
					elif pos_key_strand_N0 == '-':
						#if int(pos_key_list_N1[3])+1 == int(pos_key_list_N0[2]):
						label_Next = 'Next_continue'
						tmp_region_list_ref.append(pos_key_list_N0[0])
						tmp_region_list_ref.append(pos_key_list_N1[1])
						tmp_region_list_ref = [tmp_region_list_ref[0],tmp_region_list_ref[-1]]
						tmp_region_list_alt.append(pos_key_list_N0[3])
						tmp_region_list_alt.append(pos_key_list_N1[2])	
						tmp_region_list_alt = [tmp_region_list_alt[0],tmp_region_list_alt[-1]]
					#		print(num,x,pos_key_N0,'Next_continue',tmp_region_list_ref,tmp_region_list_alt)
						#else:
						#	label_Next = 'Warning2:Next_continueBurQueryBreak'
							#print(num,x,pos_key_N0,label_Next,tmp_region_list_ref,tmp_region_list_alt)
						#	tmp_region_list_ref,tmp_region_list_alt = [],[]
				else:
					label_Next = 'Next_break'  #Break 
					tmp_region_list_ref='\t'.join(tmp_region_list_ref)
					tmp_region_list_alt='\t'.join(tmp_region_list_alt)
					block += 1
					if pos_key_strand_N0 == "-":
						inversion = 'Inversion'
					else:
						inversion = 'Normal'
					pri=f"Block{block}\t{x}\t{tmp_region_list_ref}\t{alt_chr}\t{tmp_region_list_alt}\t{pos_key_strand_N0}\t{inversion}"
					print(pri)
					block_dic[x][f"Block{block}"] = pri
					#print(f"Block{block}",x,tmp_region_list_ref,tmp_region_list_alt,pos_key_strand_N0,label_Next)
					#print(num,x,pos_key_N0,'Next_break',tmp_region_list_ref,tmp_region_list_alt)
					tmp_region_list_ref,tmp_region_list_alt = [],[]
			elif pos_key_strand_N0 != pos_key_strand_N1:
					label_Next = 'Next_Strand_break'  #Break
					tmp_region_list_ref='\t'.join(tmp_region_list_ref)
					tmp_region_list_alt='\t'.join(tmp_region_list_alt)
					block += 1
					if pos_key_strand_N0 == "-":
						inversion = 'Inversion'
					else:
						inversion = 'Normal'
					pri=f"Block{block}\t{x}\t{tmp_region_list_ref}\t{alt_chr}\t{tmp_region_list_alt}\t{pos_key_strand_N0}\t{inversion}"
					print(pri)
					block_dic[x][f"Block{block}"] = pri
					#print(f"Block{block}",x,tmp_region_list_ref,tmp_region_list_alt,pos_key_strand_N0,label_Next)
					#print(num,x,pos_key_N0,label_Next,tmp_region_list_ref,tmp_region_list_alt)
					tmp_region_list_ref,tmp_region_list_alt = [],[]
dic_ref={}

print('############################################')
for ch in block_dic.keys():
	if list(block_dic[ch]) == []:continue
	block_list = list(block_dic[ch])
	for i in range(len(block_list)):
		if len(block_list) == 1:
			block0 = block_list[i]
			block_info0 = block_dic[ch][block0]
			block_info_list0 = block_info0.split()
			print(block_info0)
		else:
			if i == 0:
				block0 = block_list[i]
				block1 = block_list[i+1]
				block_info0, block_info1 = block_dic[ch][block0], block_dic[ch][block1]
				block_info_list0, block_info_list1 = block_info0.split(), block_info1.split()
				block0_Qsta, block0_Qend, block1_Qsta, block1_Qend=block_info_list0[5],block_info_list0[6],block_info_list1[5],block_info_list1[6]
				block0_strand, block1_strand = block_info_list0[7], block_info_list1[7]
				ref_chr, alt_chr = block_info_list0[1],block_info_list0[4]
				if ref_chr != alt_chr:
					if block0_strand == "+":
						Trans = "TransChr"
						block_info_list0[8] = Trans
						block_info0 = '\t'.join(block_info_list0)
						print(block_info0)
					elif block0_strand == "-":
						Trans = "Inversion/TransChr"
						block_info_list0[8] = Trans
						block_info0 = '\t'.join(block_info_list0)
						print(block_info0)
				else:
					if block0_strand == "+":
						if block1_strand == "+":
							judge = int(block1_Qsta) - int(block0_Qend)
							if judge < 0:
								Trans = "Trans"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
						elif block1_strand == "-":
							judge =  int(block1_Qend) - int(block0_Qend)
							if judge < 0:
								Trans = "Trans"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
					if block0_strand == "-":
						if block1_strand == "+":
							judge =  int(block1_Qsta) - int(block0_Qsta)
							if judge < 0:
								Trans = "Trans/Inversion"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
						elif block1_strand == "-":
							judge =  int(block1_Qend) - int(block0_Qsta)
							if judge < 0:
								Trans = "Trans/Inversion"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
					pass
			else:
				block0 = block_list[i]
				block1 = block_list[i-1]
				block_info0, block_info1 = block_dic[ch][block0], block_dic[ch][block1]
				block_info_list0, block_info_list1 = block_info0.split(), block_info1.split()
				block0_Qsta, block0_Qend, block1_Qsta, block1_Qend=block_info_list0[5],block_info_list0[6],block_info_list1[5],block_info_list1[6]
				block0_strand, block1_strand = block_info_list0[7], block_info_list1[7]
				ref_chr, alt_chr = block_info_list0[1],block_info_list0[4]
				if ref_chr != alt_chr:
					if block0_strand == "+":
						Trans = "TransChr"
						block_info_list0[8] = Trans
						block_info0 = '\t'.join(block_info_list0)
						print(block_info0)
					elif block0_strand == "-":
						Trans = "Inversion/TransChr"
						block_info_list0[8] = Trans
						block_info0 = '\t'.join(block_info_list0)
						print(block_info0)
				else:

					if block0_strand == "+":
						if block1_strand == "+":
							judge = int(block0_Qsta) - int(block1_Qend)
							if judge < 0:
								Trans = "Trans"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
						elif block1_strand == "-":
							judge =  int(block0_Qsta) - int(block1_Qsta)
							if judge < 0:
								Trans = "Trans"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
					if block0_strand == "-":
						if block1_strand == "+":
							judge =  int(block0_Qend) - int(block1_Qend)
							if judge < 0:
								Trans = "Trans/Inversion"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
						elif block1_strand == "-":
							judge =  int(block0_Qend) - int(block1_Qsta)
							if judge < 0:
								Trans = "Trans/Inversion"
								block_info_list0[8] = Trans
								block_info0 = '\t'.join(block_info_list0)
								print(block_info0)
							else:
								print(block_info0)
					



