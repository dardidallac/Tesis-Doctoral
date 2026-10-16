#!/usr/bin/env python3
'''This script generates variant alleles from a VCF file and a reference genome.'''
NAME = "vcf_2_genes.py"
VERSION = "1.1"
REF = 'Not published'
CITE = 'https://github.com/edgardomortiz/vcf2phylip'
GITHUB="https://github.com/maurijlozano"


#Imports
import argparse, os, sys
#from io import StringIO 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

#functions
def parseArgs():
    '''
    Argument parsing is done.
    '''
    parser = argparse.ArgumentParser(description='This script generates variant alleles from a VCF file and a reference genome.')
    # genome files
    parser.add_argument("-i", "--InputGenome",help="Reference genome in fasta format", dest="genome", action='store', required=True)
    parser.add_argument("-v", "--vcfTable",help="variant calling VCF table.", dest="vcf", action='store', required=True)
    parser.add_argument("-c", "--GeneChromosome",help="Gene chromosome location, as in GFF3/VCF file.", dest="chr", action='store', required=False)
    parser.add_argument("-s", "--GeneStart",help="Gene start position, as in GFF3 file.", dest="start", action='store', required=False)
    parser.add_argument("-e", "--GeneEnd",help="Gene end position, as in GFF3 file.", dest="end", action='store', required=False)
    parser.add_argument("--Exon_coordinates",help="Exon start and end coordinates specified as: Exon_1_Start-Exon_1_End,Exon_2_Start-Exon_2_End,...,EnS-EnE. Eg. 1-30,45-60....", dest="exon_coord", action='store', required=False)
    parser.add_argument("--Strand",help="DNA strand, '1' for direct, '0' for reverse.", dest="strand", action='store', required=False)
    parser.add_argument('-m', "--multiple_entries",help="Use this option if you want to extract multiple regions. -m <TSV> Tab separated table with the following fields: chr, start, end, strand, gene structure: Exon_1_Start-Exon_1_End,Exon_2_Start-Exon_2_End,...,EnS-EnE.", dest="table", action='store', required=False)
    parser.add_argument("-t","--Translate",help="Translate DNA to amino acid.", dest="translate", action='store_true')
    parser.add_argument("--CDS",help="Extracts the coding sequence. Requires the --Exon_coordinates argument or the Exon_coordinates as fifth field in the table.", dest="CDS", action='store_true')
    parser.add_argument("--Clean",help="Deletes temporary files.", dest="clean", action='store_true')
    args = parser.parse_args()   
    return(args)

def printSoftName():
    print("\n\n\n")
    print("   ************************************************")
    print("   *****   vcf_2_genes: Get variant alleles from a") 
    print("   *****   VCF file and a reference genome.")
    print("   *****   Version: "+str(VERSION))
    print("   *****   Developed by Mauricio J. Lozano")
    print("   *****   github.com/maurijlozano")
    print("   ************************************************")
    print("   Please cite: "+REF)
    print("   Please cite: "+CITE)
    print("   Downloaded from: "+GITHUB)
    print("\n\n\n")

def get_ref_gene(genome,chr,start,end):
    '''
    Get the ref gene sequence
    '''
    for record in SeqIO.parse(genome,'fasta'):
        if record.id == chr:
            ref_gene_seq = str(record.seq[start:end]).upper()
    return(ref_gene_seq)

def get_sample_names(vcf):
    '''
    Get the sample names and number
    '''
    with open(vcf,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("#CHROM"):
                rec = line.split("\t")
                sample_names = [rec[i].replace("./", "") for i in range(9, len(rec))]
                break
    num_samples = len(sample_names)    
    return(sample_names,num_samples)

def create_seq_dict(ref_gene_seq,sample_names):
    seq_dict = {}
    for sample in sample_names:
        seq_dict[sample] = ref_gene_seq
    return(seq_dict)

#insert gaps...
def align_nt(seq1,seq2):
    seqa = []
    seqb = []
    j=0
    if len(seq1) < len(seq2):
        for i in range(max(len(seq1),len(seq2))):
            f'{i}  {j}'
            if j < len(seq1): 
                if seq1[j] == seq2[i]:
                    seqa.append(seq1[j])
                    seqb.append(seq2[i])
                    j+=1
                else:
                    seqa.append('-')
                    seqb.append(seq2[i])
            else:
                seqa.append('-')
                seqb.append(seq2[i])
        return(''.join(seqa),''.join(seqb))
    elif len(seq2) < len(seq1):
        seq3 = seq2
        seq4 = seq1
        for i in range(max(len(seq3),len(seq4))):
            if j < len(seq3): 
                if seq3[j] == seq4[i]:
                    seqa.append(seq3[j])
                    seqb.append(seq4[i])
                    j+=1
                else:
                    seqa.append('-')
                    seqb.append(seq4[i])
            else:
                seqa.append('-')
                seqb.append(seq4[i])
        return(''.join(seqb),''.join(seqa))
    else:
        return(seq1,seq2)

def gt_to_string(gt,num_to_allele,ref_allele):
    if '/' in gt:
        gt_list = gt.split('/')
    elif '|' in gt:
        gt_list = gt.split('|')
    if len(gt_list) == 2:
        gt_bases = list(set([num_to_allele[n] for n in gt_list]))
        len_ref_allele = len(ref_allele)
        if len_ref_allele > 1:
            #indel. If there are to GTs in the alt, one of them including INDELS, we chose the SNP, unless it is a /3 bases insertion or deletion.
            #3* deletion/insertion. If multiple correct alleles, choose the largest.
            indel = [align_nt(ref_allele,g) for g in gt_bases if ((len(g)-1)%3 == 0 )]
            if len(indel) > 0:
                max_len = 0
                index = 0
                for i, r_a in enumerate(indel):
                    if len(r_a[0]) > max_len:
                        index = i
                        max_len = len(r_a[i])
                return(indel[index])
            else:
                indel = [align_nt(ref_allele,g) for g in gt_bases]
                if len(indel) > 0:
                    max_len = 0
                    index = 0
                    for i, r_a in enumerate(indel):
                        if len(r_a[i]) > max_len:
                            index = i
                            max_len = len(r_a[i])
                    return(indel[index])
                else:
                    return(ref_allele,'N'*len_ref_allele)
        else:    
            if len(gt_bases) > 1:
                if gt_bases[0] != ref_allele:
                    return(ref_allele,gt_bases[0])
                else:
                    return(ref_allele,gt_bases[1])
            else:
                return(ref_allele,gt_bases[0])
    else:
        sys.exit(f'Sorry, vcf_2_genes.py only works with diploid genomes.')

def extract_vcf(table,vcf):
    #split vcf, generates one file per gene
    chrs = list(table.iloc[:,0])
    starts = list(table.iloc[:,1])
    ends = list(table.iloc[:,2])
    strands = list(table.iloc[:,3])
    data = list(zip(chrs,starts,ends,strands))
    try:
        temp_files = []
        temp_file_names = []
        for i,row in table.iterrows():
            temp_file_name = f'{row[0]}_{row[1]}_{row[2]}_{row[3]}.tmp'
            temp_file_names.append(f'{row[0]}_{row[1]}_{row[2]}_{row[3]}.tmp')
            temp_files.append(open(temp_file_name,'w+'))
        #    
        with open(vcf,'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith(f'#'):
                    vcf_rec = line.split('\t')
                    chr = vcf_rec[0]
                    snp_pos = int(vcf_rec[1])
                    for i,d in enumerate(data):
                        if (d[0] == chr) and (snp_pos >= d[1]) and (snp_pos <= d[2]):
                           f = temp_files[i]
                           f.write(line+'\n')
    except:
        print('Error: I/O operation failed.') 
    finally:
        for f in temp_files:
            f.close()
    return(temp_file_names)

def correct_exon_coords(exon_rel_coord_dict,rel_ref_coord,added_gaps):
    for exon in exon_rel_coord_dict:
        if rel_ref_coord < exon_rel_coord_dict[exon][0]:
            exon_rel_coord_dict[exon][0] += added_gaps
            exon_rel_coord_dict[exon][1] += added_gaps
        elif (rel_ref_coord > exon_rel_coord_dict[exon][0]) and (rel_ref_coord < exon_rel_coord_dict[exon][1]):
            exon_rel_coord_dict[exon][1] += added_gaps
    return(exon_rel_coord_dict)

def main(genome, vcf, chr, start, end, strand, ref_gene_seq, sample_names,num_samples, seq_dict,exon_coord,cds):
    #######
    if len(exon_coord) > 0:
        exon_rel_coord_dict = { k:[int(k.split('-')[0]) -1 - start, int(k.split('-')[1]) - start] for k in exon_coord}
        multiple_exons = True
    else:
        multiple_exons = False
    #######
    with open(vcf,'r') as f:
        gap_counter = 0
        for line in f:
            line = line.strip()
            if line.startswith('#') and not line.startswith('##'):
                header = line[1:].split('\t')
            elif line.startswith(f'{chr}\t'):
                vcf_rec = line.split('\t')
                snp_pos = int(vcf_rec[1])
                rel_mut_index = snp_pos -1 - start + gap_counter #adds gaps for position correction!
                if (snp_pos >= start) and ((snp_pos <= end)):
                    if (len(vcf_rec) - 9) != num_samples:
                        print(f'Print skipping {chr} - {snp_pos}, missing data...')
                        continue
                    ref_gt = vcf_rec[3].replace("-","*").upper()
                    ref_alleles = ref_gt.split(',')
                    if len(ref_alleles) > 1:
                        sys.exit('Reference genome must be haploid...')
                    ref_allele = ref_alleles[0]
                    alt_gt = vcf_rec[4].replace("-", "*").replace("<NON_REF>", ref_gt)
                    alleles = alt_gt.split(',')
                    num_to_allele = { str((i+1)):a for i,a in enumerate(alleles)}
                    num_to_allele['0'] = ref_alleles[0]
                    num_to_allele['.'] = 'N'
                    for i in range(9, len(vcf_rec)):
                        gt = vcf_rec[i].split(':')[0]
                        hetero = ''
                        if (len(set(gt.split('/'))) > 1) or (len(set(gt.split('|'))) > 1):
                            hetero = 'H' 
                        ref, alt = gt_to_string(gt,num_to_allele,ref_allele)
                        #update sequences, update numeration for insertions
                        if len(ref) != len(alt):
                            sys.exit(f'Error: reference and alternative alleles should be the same length after their alignment...')
                        else:
                            ref_len = len(ref_allele)
                            ref_aligned_len = len(ref)
                            #verify ref allele:
                            if ref_gene_seq[rel_mut_index:rel_mut_index+ref_len] != ref_allele:
                                print('Error: ref allele doesn\'t match with vcf file.')
                            #add gaps to ref:
                            ref_gene_seq = ref_gene_seq[0:rel_mut_index] + ref + ref_gene_seq[rel_mut_index+ref_len:]
                            #mutate sample seqs
                            seq_dict[sample_names[i-9]] = seq_dict[sample_names[i-9]][0:rel_mut_index] + alt + seq_dict[sample_names[i-9]][rel_mut_index+ref_len:]
                            #Correct exon coordinates with added gaps
                            added_gaps = (ref_aligned_len - ref_len)
                            rel_ref_coord = snp_pos -1 - start 
                            if multiple_exons:
                                exon_rel_coord_dict = correct_exon_coords(exon_rel_coord_dict,rel_ref_coord,added_gaps)
                            #
                            gap_counter += added_gaps
                else:
                    continue
            else:
                continue
    #do for gene in genes
    out_file_name = os.path.splitext(genome)[0] + f'_{chr}_{start}_{end}.fasta'
    with open(out_file_name, 'w+') as f:
        f.write(f'>Ref_Seq_{chr}_{start}_{end}\n{ref_gene_seq}\n')
        for seq in seq_dict:
            f.write(f'>{seq}_{chr}_{start}_{end}\n{seq_dict[seq]}\n')
    if translate:
        from Bio import Seq
        #translate protein to stop
        out_file_name_faa = os.path.splitext(genome)[0] + f'_{chr}_{start}_{end}.faa'
        out_file_name_faa2 = os.path.splitext(genome)[0] + f'_{chr}_{start}_{end}_with_stops.faa'
        if cds:
            out_file_name_fna = os.path.splitext(genome)[0] + f'_{chr}_{start}_{end}.fna'
        with open(out_file_name_faa, 'w+') as f, open(out_file_name_faa2, 'w+') as g, open(out_file_name_fna, 'w+') as h:
            for record in SeqIO.parse(out_file_name, 'fasta'):
                seq = record.seq
                if multiple_exons:
                    coding_seq = Seq.Seq('')
                    for exon in exon_rel_coord_dict:
                        coding_seq = coding_seq + seq[exon_rel_coord_dict[exon][0]:exon_rel_coord_dict[exon][1]]
                else:
                    coding_seq = seq
                if not strand:
                    coding_seq = coding_seq.reverse_complement()
                if cds:
                    h.write(f'>{record.id}\n{coding_seq}\n')
                if len(coding_seq)%3 > 0:
                    add_N = 3 - len(coding_seq)%3 
                    coding_seq = coding_seq + 'N'*add_N
                aa_seq = coding_seq.translate(to_stop=True)
                f.write(f'>{record.id}\n{aa_seq}\n')
                #
                aa_seq2 = coding_seq.translate()
                g.write(f'>{record.id}\n{aa_seq2}\n')
    success = f'The variants present in chromosome: {chr} between start: {start} and end: {end} coordinates were successfully extracted.'
    return(success)


if __name__ == '__main__':
    #parse_vcf_table
    #args
    args = parseArgs()
    #
    genome = args.genome
    vcf = args.vcf
    translate = args.translate
    cds = args.CDS
    clean = args.clean
    #optional
    if args.table:
        table = pd.read_table(args.table, header=None)
        table.iloc[:,0] = table.iloc[:,0].astype('str')
        tmp_file_names = extract_vcf(table,vcf)
        i=0
        for _,row in table.iterrows():
            vcf_file = tmp_file_names[i]
            if os.path.getsize(vcf_file) > 0:
                chr = row.iloc[0]
                start = row.iloc[1]
                end = row.iloc[2]
                strand = row.iloc[3]
                if len(row) > 4:
                    exon_coord = row[5].split(',')
                else:
                    exon_coord = ''
                ref_gene_seq = get_ref_gene(genome,chr,start,end)
                sample_names,num_samples = get_sample_names(vcf)
                seq_dict = create_seq_dict(ref_gene_seq,sample_names)
                to_print = main(genome, vcf, chr, start, end, strand, ref_gene_seq, sample_names,num_samples, seq_dict,exon_coord,cds)
                print(to_print)
            i+=1
        if clean:
            for f in tmp_file_names:
                os.remove(f)
    elif (args.chr) and (args.start) and (args.end) and (args.strand):
        chr = args.chr
        start = int(args.start) - 1
        end = int(args.end)
        strand = int(args.strand)
        if args.exon_coord:
            exon_coord = args.exon_coord.split(',')
        else:
            exon_coord = ''
        ref_gene_seq = get_ref_gene(genome,chr,start,end)
        sample_names,num_samples = get_sample_names(vcf)
        seq_dict = create_seq_dict(ref_gene_seq,sample_names)
        main(genome, vcf, chr, start, end, strand, ref_gene_seq, sample_names,num_samples, seq_dict,exon_coord,cds)
    else:
        sys.exit(f'Error: Either a table containing chr, start, end and strand for the genes to analyze; or a chr, start, end, and strand for a single gene are required.')

