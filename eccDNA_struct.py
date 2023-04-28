import os
import pandas as pd
import pysam
from Bio import SeqIO
from statistics import median, mean
import random
import numpy as np
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt

def mapppingBamSort(reads, genome_fasta, output = None, 
                    mm2 = 'minimap2', samtools_path="samtools", bamtools_path="bamtools"):
    if not output:
        os.system('{3} -ax map-ont --splice-flank=no -t 150 {0} {1} > {2}.sam'.format(genome_fasta, reads, reads, mm2))
        sam_file = '{0}.sam'.format(reads)
    else:
        os.system('{3} -ax map-ont --splice-flank=no -t 150 {0} {1} > {2}.sam'.format(genome_fasta,  reads, output, mm2))
        sam_file = '{0}.sam'.format(output)

    bam_file = sam_file.rsplit('.', 1)[0] + ".bam"
    sort_bam_file = bam_file.rsplit(r'/', 1)[0] + r"/sorted_" + bam_file.rsplit(r'/', 1)[1]

    ##sam to bam
    samview = '{2} view -Sb {0} > {1}'.format(sam_file, bam_file, samtools_path)
    print(samview)
    os.system(samview)

    ##sort bam
    sortbam = '{2} sort -o {1} -@ 100 {0} '.format(bam_file, sort_bam_file, samtools_path)
    print(sortbam)
    os.system(sortbam)

    ##index
    os.system('{1} index -@ 100 {0}'.format(sort_bam_file, samtools_path))

    return(sort_bam_file)

def run_TH(fastq, outTab, threads=100, min_period=50, min_reps=2):
    print("TideHunter has started ....")
    os.system('{0} -f 1 {1} -t {3} -p {4} -c {5} -l > {2}'.format('TideHunter', fastq, outTab, threads, min_period, min_reps))
    print("TideHunter has finished ..... Generation of the consensus sequence is done")
    return outTab

def truncateIDs_and_selectMonomers(fasta, min_monomer_length = 500, min_rep = 2):
    outname = fasta.rsplit('.',1)[-2] +'.truncated_id.fasta' 
    cnt_all = 0
    cnt_selected = 0
    with open(outname,'w') as out:
        for seq in SeqIO.parse(fasta, 'fasta'):
            sp = seq.id.split("_") # >readName_repN_copyNum readLen_start_end_consLen_aveMatch_fullLen_subPos
            copyNum = sp[6]
            cnt_all += 1
            # if 'cb8fc1ab-6e6d-4414-9dd2-5b1b5ec14e56' in seq.id:
            #     print(seq.id)
            if len(seq.seq) >= min_monomer_length and float(copyNum) >= min_rep:
                cnt_selected += 1
                seq.id = seq.id.split('_rep')[0]
                seq.description = ''
                SeqIO.write(seq,out,'fasta' )
    print(f'{cnt_selected} sequences were kept of {cnt_all}')
    return(outname)
def getReadsFromBam(coordinate, bam_file, out_fasta):
    cigars = []
    ## determine start and end of target TE
    coordinate = coordinate.replace(',','')
    sp = coordinate.split(':')
    chrom = sp[0]
    start, end = [int(i) for i in sp[1].split('..')]
    cigars = []
    fo = open(out_fasta,'w')
    cnt = 0
    for algns in pysam.AlignmentFile(bam_file, 'rb').fetch(chrom, start, end):
        if algns.cigartuples and algns.reference_name != None:
            if algns.reference_name == chrom:
                if not algns.is_supplementary and not algns.is_secondary:
                    st, en =  min(algns.reference_start,  algns.reference_end), max(algns.reference_start,  algns.reference_end)
                    fo.write(f'>{algns.query_name}\n{algns.query_sequence}\n')
                    cnt += 1
    print('Number of sequences selected', cnt)
    fo.close()
    return out_fasta

def getTargetSeq(chrom, start, end, fasta_file, oufile_name):
    with open(oufile_name, 'w') as out:
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            if seq.id == chrom:
                seq_target = seq.seq[start:end]
                out.write(f'>{chrom}_{start}_{end}\n{seq_target}\n')
    return [oufile_name, len(seq_target)]

def getOverlapInfo(previous_coord, current_read_hit, perc_small_overlap = 0.8, over_bp=10):
    current_s, current_e = current_read_hit
    prev_s, prev_e = previous_coord

    len_cur, len_prev = current_e - current_s, prev_e - prev_s
    # junction
    if (abs(current_s - prev_e) <= over_bp) and (current_e > prev_e):
        return True
    ## new hit is insight
    elif ((prev_s - 10) < current_s) and ((prev_e + 10) > current_e):
        return False
    #the end of current hit longer than previous
    elif (abs(prev_s - current_s) <= 10) and (current_e > prev_e):
        return 'replace'
    elif (abs(prev_e - current_e) <= 10) and (current_s < prev_s):
        return 'replace'
    elif (current_s < prev_s)and (current_e > prev_e):
        return 'replace'
    #significant overlap
    elif (prev_s < current_s) and (prev_e < current_e) and (((prev_e-current_s)/len_cur) > perc_small_overlap):
        return False
    else:
        return True
 
  
# def removeRedundancy(list_hit, min_over_to_include_point = 50):    
#     new_coords = []
#     list_hit = sorted(list_hit)
#     for i, points in enumerate(list_hit):
#         if i != 0:
#             if (points - new_coords[-1]) > 10:
#                 new_coords.append(points)
#             else:
#                 #do not remove last element if it is only on the least
#                 if len(new_coords) > 1:
#                     new_coords.pop()
#                 if i == len(list_hit) - 1:
#                     new_coords.append(points)
#         else:
#             new_coords.append(points)
#     return new_coords

def removeRedundancy(list_hit, min_over_to_include_point = 50):
    list_hit = sorted(list_hit, key=lambda x: x[0])
    compressed_list = []
    for i, coord in enumerate(list_hit):
        if i == 0:
            compressed_list.append(coord)
            
        ##non overlapped
        elif coord[0] >= compressed_list[-1][1] + min_over_to_include_point:
            compressed_list.append(coord)
        else:
            merged_list = compressed_list[-1]  + coord
            compressed_list[-1] = [min(merged_list), max(merged_list)]
        
    return unzipLists(compressed_list)

def unzipLists(listOflists):
    toret = []
    for i in listOflists:
        for elem in i:
            toret.append(elem)
    return(sorted(toret))

def getReadLenFromBlasttab(blastTab):
    read_len = {}
    blast_hits = pd.read_csv(blastTab, sep='\t', header=None)
    blast_hits.columns = ['read_id', 'ref', 'read_len', 'ref_len', 'bitscore', 'qstart', 'qend', 'sstart', 'send']
    
    lens = list(blast_hits['read_len'])
    for i, reads in enumerate( list(blast_hits['read_id']) ):
        read_len[reads] = lens[i]
    
    return read_len
        
def isBLASTtabNotEmpty(blastTab):
    cnt = 0
    with open (blastTab) as tab:
        for lines in tab:
            cnt+1
    if cnt==0:
        return False
    else:
        return True
    
def defineBLASTpointsOnTE(blastTab, min_hit_len = 50):
    if not isBLASTtabNotEmpty:
        print('ERROR!!! NO BLAST HITS FOUND!')
        return False
    blast_hits = pd.read_csv(blastTab, sep='\t', header=None)
    blast_hits.columns = ['read_id', 'ref', 'read_len', 'ref_len', 'bitscore', 'qstart', 'qend', 'sstart', 'send']
    # iterate through reads
    per_read_points_ref = {}
    per_read_points_read_itself = {}
    cnt_empty = 0
    for group_name, df_group in blast_hits.sort_values(['read_id','qstart'], ascending=True).groupby('read_id'):
        #print(group_name)
        coords_in_read = []
        coords_in_ref = []
        # iterate through hits for each read
        for row_index, row in df_group.iterrows():
            ref_s, ref_e = min(row['sstart'], row['send']), max(row['sstart'], row['send'])
            if (ref_e - ref_s) >= min_hit_len:
                if coords_in_read:
                    # over will be True if two query hits are junctions
                    for i, prev_coord in enumerate(coords_in_read):
                        over = getOverlapInfo(prev_coord, [row['qstart'], row['qend']], perc_small_overlap = 0.5, over_bp=10)
                        if over == False or over == 'replace':
                            break
                        # if row['read_id'] == '83712fed-4d78-47ec-b395-f12949fda062':
                        #        print(over, prev_coord, [row['qstart'], row['qend']])

                    ## add new intervals
                    if over == True:
                        coords_in_read.append([row['qstart'], row['qend']])
                        coords_in_ref.append([ref_s, ref_e])
                        # if row['read_id'] == '2828ae15-dc3a-49e5-9b46-ec7f0810b3c8':
                        #     print(coords_in_read, coords_in_ref)

                    ## replace one of the intervals
                    elif over == 'replace':
                        coords_in_read[i][0], coords_in_read[i][1] = row['qstart'], row['qend']
                        coords_in_ref[i][0],coords_in_ref[i][1] = row['sstart'], row['send']


                else:
                    coords_in_read.append([row['qstart'], row['qend']])
                    coords_in_ref.append([ref_s, ref_e])
                
        if group_name == '014b1d5d-b775-48b3-8173-f15618eb78eb':
            print(coords_in_read, coords_in_ref)
        #print(group_name)
        # print(coords_in_read)
            print(coords_in_ref)
            print(removeRedundancy(coords_in_ref))
        per_read_points_ref[row['read_id']] = removeRedundancy(coords_in_ref)
        per_read_points_read_itself[row['read_id']] = removeRedundancy(coords_in_read)
        # if len(per_read_points_ref[row['read_id']]) == 1:
        #     print(row['read_id'], per_read_points_ref[row['read_id']])
        #     cnt_empty += 1
    print(cnt_empty)   
    return [per_read_points_ref, per_read_points_read_itself]

def getMatrixAndDrawHist(genome, fastq, TE_coordinate, 
                         prefix, output_folder, min_reps_TH=2, window = 50, runTH=False,
                         min_occurence_for_path_to_include_inHM = 3):
    ##initialize variables
    os.system(f'mkdir {output_folder}')
    collected_reads = f'{output_folder}/collected_TE_reads'
    os.system(f'mkdir {collected_reads}')
    #1. run TideHunter
    print(f'Step 1 of 10. TideHunter selection of concatemer reads with {min_reps_TH} or more monomer repeats....')
    if runTH:
        contigs_all = run_TH(fastq, fastq + ".TH.cons.fasta", threads=100, min_period=500, min_reps=min_reps_TH)
    else:
        print(f'TideHunter is not running. Using preexisted file: {fastq}.TH.cons.fasta')
    
    #2. mapping
    print('Step 2 of 10. Mapping selected reads to the genome....')
    #truncate read ids for minimap2 
    contigs = truncateIDs_and_selectMonomers(fastq+ ".TH.cons.fasta",
                                        min_monomer_length = 500, min_rep = 3)
    ## map ecc_finder contigs 
    bam_file = mapppingBamSort(contigs, genome, f'{output_folder}/TH_selected_{prefix}')
    print(bam_file)
    
    #3. get reads from coordinate regions
    print('Step 3 of 10. Getting reads from TE coordinate....')
    coordinate = TE_coordinate.replace(',','')
    sp = coordinate.split(':')
    chrom = sp[0]
    start, end = [int(i) for i in sp[1].split('..')]

    out_fasta_target = f'{collected_reads}/{TE_coordinate}.fasta'
    out_fasta = out_fasta_target.rsplit('.', 1)[0] + "_reads.fasta"

    selected_reads = getReadsFromBam(coordinate, bam_file, out_fasta)

    target_fasta, target_len = getTargetSeq(chrom, start, end, genome, out_fasta_target)
    print('Target length:', target_len)
    
    #4. BLAST of selected reads vs target sequence
    print('Step 4 of 10. BLAST the target reads vs target sequence....')

    os.system(f'makeblastdb -dbtype nucl -in {target_fasta} -out {target_fasta}')
    os.system(f'blastn -word_size 25 -outfmt "6 qseqid sseqid qlen slen bitscore qstart qend sstart send" -query {selected_reads} -db {target_fasta} -out {selected_reads}.blasttab')
    
    
    #5. Determine non-overlapped TH monomer fragments and the coordinates (0-based) of corresponding ref sequences
    print('Step 5 of 10. Determination of non-overlapped TH monomer fragments and the coordinates (0-based) of corresponding target sequence....')

    per_read_points_ref, per_read_points_itself = defineBLASTpointsOnTE(f'{selected_reads}.blasttab')
    
    
    #####################################################################################
    ########################### PLOTTING ################################################
    #####################################################################################
    
    #6. Make histogram of occurence of individual path
    print('Step 6 of 10. Plotting the histogram')
    cm = 1/2.54  # centimeters in inches
    fig, axs = plt.subplots(2, 2, gridspec_kw={'width_ratios': [1, 1]})
    fig.set_figheight(30*cm)
    fig.set_figwidth(30*cm)
    
    ## Number of ecc dna alignments in ref
    h = [int(len(per_read_points_ref[i])/2) for i in per_read_points_ref]
    numb_var = pd.Series(h)
    
    
    # Plotting the histogram
    sns.histplot(data = numb_var, ax=axs[0,0]).set(title="Number of distinct target fragments for each eccDNA")
    #plt.savefig(f'{output_folder}/{coordinate}_{prefix}_Nalgn_hist.png')
    
    # 7. eccDNA length distribution
    print("Step 7 of 10. Plotting eccDNA length distribution....")
    
    read_len = getReadLenFromBlasttab(f'{selected_reads}.blasttab')
    read_coverage = []
    len_ecc_dna_mapped = {}
    l = []
    for reads in per_read_points_itself:
        len_ecc = 0
        for i, points in enumerate(per_read_points_itself[reads]):
            if (i + 1) % 2 == 0:
                len_ecc += points - per_read_points_itself[reads][i-1]
        len_ecc_dna_mapped[reads] = len_ecc
        l.append(len_ecc)
        cov = (len_ecc*100) / read_len[reads]
        # print(cov)
        # # if cov > 110:
        # #     print(per_read_points_itself[reads], cov, len_ecc, read_len[reads], reads)
        read_coverage.append(cov)
    
    len_var = pd.Series(l)
    #plt.savefig(f'{output_folder}/{coordinate}_{prefix}_eccDNA_length.png')

    # Plotting the histogram
    sns.histplot(data = len_var, ax=axs[0,1], kde=True, color='#ffb347').set(title="eccDNA length based on reference alignments")
    
    sns.histplot(data = pd.Series(read_coverage), bins=20, ax=axs[1,0]).set(title="Monomer coverage by target hits")
    
    #8. Group reads based on the fragments position along target seq
    print("Step 8 of 10. Grouping the reads based on the fragments position along target seq....")

    per_read_points_ref
    ##create windows 
    fragments = [i for i in range(0, target_len, window)]

    #find groups of circles by comparing patterns on chromosomes
    fragment_positions = {}

    for idx, reads in enumerate(per_read_points_ref):
        all_idx_of_read = []
        for point in per_read_points_ref[reads]:
            idx = np.searchsorted(fragments, point)
            all_idx_of_read.append(fragments[idx-1])

        joined_all_idx_of_read = "_".join([str(i) for i in all_idx_of_read])
        if joined_all_idx_of_read not in fragment_positions:
            fragment_positions[joined_all_idx_of_read] = [reads]
        else:
            fragment_positions[joined_all_idx_of_read].append(reads)

    
    #9. PLOT Number of occurence of path
    print("Step 9 of 10. Plotting number of reads....")

    min_occurence = min_occurence_for_path_to_include_inHM
    var_occurence = []
    selected_ecc_variants = {}
    total_reads_cnt = 0
    selected_reads_cnt = 0
    for path in fragment_positions:
        if len(fragment_positions[path]) >= min_occurence:
            selected_ecc_variants[path] = len(fragment_positions[path])
            var_occurence.append(len(fragment_positions[path]))
            selected_reads_cnt+= len(fragment_positions[path])
        total_reads_cnt+=len(fragment_positions[path])
    print('Number of selected variants:', len(selected_ecc_variants))
    print('Number of reads (ALL):', total_reads_cnt)
    print('Number of reads selected:', selected_reads_cnt)

    occurence_var = pd.Series(var_occurence)
    
    # Plotting the histogram
    sns.histplot(data = occurence_var, ax=axs[1,1], bins=20).set(title="Number of occurrence for paths")
    plt.savefig(f'{output_folder}/{coordinate}_{prefix}_path_occurence.png')
    
    #10. Prepare matrix for HeatMap
    print(f"Step 10 of 10. Prepare matrix ({output_folder}/MATRIX_{coordinate}_{prefix}_matrix.tab) for HeatMap....")

    matrix_ = []
    path_names = []
    cnt = min_occurence_for_path_to_include_inHM
    for path in selected_ecc_variants:
        if path:
            #print(path.split('_'))
            zero_list = [0 for i in fragments]
            occurence = selected_ecc_variants[path]
            coords = [int(i) for i in path.split('_')]
            paired_coord = []
            for i, points in enumerate(coords):
                if (i+1)%2 == 0:
                    paired_coord.append([coords[i-1], points])
            #print(path, paired_coord, occurence)
            for idx, window in enumerate(fragments):
                for pairs in paired_coord:
                    if window >= pairs[0] and window < pairs[1]:
                        zero_list[idx] = occurence
            #print(zero_list)
            matrix_.append(zero_list)
            path_names.append(path)
        else:
            cnt+1
    matrix_pd = pd.DataFrame(matrix_)
    matrix_pd.columns = fragments
    matrix_pd.index = path_names

    print("Empty", cnt)
    ## WRITE FINAL MATRIX ##
    matrix_pd.to_csv(f'{output_folder}/MATRIX_{coordinate}_{prefix}_matrix.tab', sep='\t', index=True)
