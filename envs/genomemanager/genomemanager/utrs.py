import numpy as np


def is_suspected_to_not_be_5utr(start, end, start_start_codon, end_start_codon, strand):


    if strand == '+':
        if np.isnan(end_start_codon):
            return False
        if start <= end_start_codon:
            if end <= end_start_codon:
                return False
            else:
                return True
        else:
            return True
    else:
        if np.isnan(start_start_codon):
            return False
        if end >= start_start_codon:
            if start >= start_start_codon:
                return False
            else:
                return True
        else:
            return True

def is_suspected_to_not_be_3utr(start, end, start_stop_codon, end_stop_codon, strand):

    if strand == '+':
        if np.isnan(start_stop_codon):
            return False
        if end >= start_stop_codon:
            if start >= start_stop_codon:
                return False
            else:
                return True
        else:
            return True
    else:
        if np.isnan(end_stop_codon):
            return False
        if start <= end_stop_codon:
            if end <= end_stop_codon:
                return False
            else:
                return True
        else:
            return True

def remove_suspicious_5_utrs(whole_annotation, utr_annotation):
    utr_annotation['gene_id'] = utr_annotation['feature'].apply(lambda x: x.split('_utr')[0])
    annot_advanced = whole_annotation.groupby(['gene_id', 'feature', 'strand']).agg({'start': min, 'end': max}).reset_index()
    utr_codon_info = annot_advanced[annot_advanced['feature'] == 'start_codon'][['strand', 'gene_id', 'start', 'end']].merge(utr_annotation, left_on=['gene_id', 'strand'], how='right', right_on=['gene_id', 'strand'], suffixes=('_start_codon', ''))   
    utr_codon_info['is_suspicious'] = utr_codon_info.apply(lambda row: is_suspected_to_not_be_5utr(row['start'], row['end'], row['start_start_codon'], row['end_start_codon'], row['strand']), axis=1)
    suspicious_utrs = utr_codon_info[utr_codon_info['is_suspicious'] == True].index
    return utr_codon_info[utr_codon_info.index.isin(suspicious_utrs) == False]

def remove_suspicious_3_utrs(whole_annotation, utr_annotation):
    utr_annotation['gene_id'] = utr_annotation['feature'].apply(lambda x: x.split('_utr')[0])
    annot_advanced = whole_annotation.groupby(['gene_id', 'feature', 'strand']).agg({'start': min, 'end': max}).reset_index()
    utr_codon_info = annot_advanced[annot_advanced['feature'] == 'stop_codon'][['strand', 'gene_id', 'start', 'end']].merge(utr_annotation, left_on=['gene_id', 'strand'], how='right', right_on=['gene_id', 'strand'], suffixes=('_stop_codon', ''))   
    utr_codon_info['is_suspicious'] = utr_codon_info.apply(lambda row: is_suspected_to_not_be_3utr(row['start'], row['end'], row['start_stop_codon'], row['end_stop_codon'], row['strand']), axis=1)
    suspicious_utrs = utr_codon_info[utr_codon_info['is_suspicious'] == True].index
    
    return utr_codon_info[utr_codon_info.index.isin(suspicious_utrs) == False]