import os
import re
import glob
import csv
from pathlib import Path
import shutil


def update_fasta_file(file, parent_folder):
    lines = []
    with open(file, 'r') as f:
        lines = f.readlines()
    first_line = lines.pop(0)
    length = re.search(r'total_supporting_reads_(\d)', first_line).group(1)
    reads = re.search(r'total_supporting_reads_(\d)', first_line).group(1)
    lines.insert(0, '>' + parent_folder + '\n')
    with open(file, 'w') as f:
        [f.write(line) for line in lines]

    return length, reads

def get_fasta_reads(fasta_file):
    with open(fasta_file, 'r') as f:
        first_line = f.readline()
        reads = re.search(r'total_supporting_reads_(\d)', first_line).group(1)
    return reads

def process_medaka_folder(folder, medaka_id):
    print(f'Processing Medaka folder: {folder}...')
    fasta_file = os.path.join(folder, 'consensus.fasta')
    if os.path.exists(fasta_file):
        p = Path(fasta_file)
        parent_folder = p.parent.parent.name
        parent_folder_name = re.sub(r'^sample_', '', parent_folder)
        reads = get_fasta_reads(fasta_file)
        target_fasta_file = os.path.join(summary_folder, f'{parent_folder_name}-{medaka_id}-RiC{reads}.fasta')
        # shutil.move(fasta_file, target_fasta_file) # TODO
        shutil.copy(fasta_file, target_fasta_file)
        length, reads = update_fasta_file(target_fasta_file, f'{parent_folder_name}-{medaka_id}')

    p = Path(folder)
    medaka_parent_folder = str(p.parent)
    fastq_file_index = 1
    for fastq_file in sorted(glob.glob(f'{medaka_parent_folder}/reads_to_consensus_[0-9]*.fastq')):
        if os.path.exists(fastq_file):
            target_fastq_file = os.path.join(summary_folder, 'FASTQ Files', f'{parent_folder_name}-{fastq_file_index}.fastq')
            # shutil.move(fastq_file, target_fastq_file) # TODO
            shutil.copy(fastq_file, target_fastq_file)
            fastq_file_index += 1

    return f'{parent_folder_name}-{medaka_id}', length, reads

def process_folder(folder):
    medaka_id = 1
    stats = []
    for item in sorted(glob.iglob(folder + '/*', recursive=False)):
        if not os.path.isfile(item):
            p = Path(item)
            folder_name = p.name
            if folder_name.startswith('medaka'):
                filename, length, reads = process_medaka_folder(item, medaka_id)
                stats.append([filename, length, reads, medaka_id])
                medaka_id += 1

    return stats

def merge_fasta_files(folder, output_file):
    with open(os.path.join(folder, output_file), 'w') as f_out:
        for fasta_file in sorted(glob.glob(f'{folder}/*.fasta')):
            with open(fasta_file, 'r') as f:
                f_out.write(f.read())

def generate_summary(stats, writer):

    summary_totals = {
        'Total Unique Samples': 0,
        'Total Consensus Sequences': 0,
        'Total Reads in Consensus Sequences': 0,
    }
    additional_stats = []
    counted_samples = {}

    for row in stats:
        filename = row[0]
        reads = int(row[2])
        if filename.count('ONT') > 1:
            additional_stats.append(row)
        else:
            writer.writerow(row)
            unique_part = re.search(r'^([^-]+)', filename).group(1)
            if unique_part not in counted_samples:
                summary_totals['Total Unique Samples'] += 1
                counted_samples[unique_part] = 1
        summary_totals['Total Consensus Sequences'] += 1
        summary_totals['Total Reads in Consensus Sequences'] += reads

    writer.writerow([])
    writer.writerow(['Total Unique Samples', summary_totals['Total Unique Samples']])
    writer.writerow(['Total Consensus Sequences', summary_totals['Total Consensus Sequences']])
    writer.writerow(['Total Reads in Consensus Sequences', summary_totals['Total Reads in Consensus Sequences']])
    writer.writerow([])
    for row in additional_stats:
        writer.writerow(row)

# source_folder = 'H:\Temp\Test Folder'

source_folder = os.path.dirname(os.path.realpath(__file__)) # TODO
summary_folder = '__Summary__'

# Create Summary folder and subfolder
if not os.path.exists(summary_folder):
    os.mkdir(summary_folder)
if not os.path.exists(os.path.join(summary_folder, 'FASTQ Files')):
    os.mkdir(os.path.join(summary_folder, 'FASTQ Files'))

with open(os.path.join(summary_folder, 'summary.txt'), 'w', encoding='utf-8') as f_out:
    writer = csv.writer(f_out, delimiter='\t', lineterminator='\n')
    writer.writerow(['Filename', 'Length', 'Reads in Consensus', 'Multiple'])

    stats = []
    for item in sorted(glob.iglob(source_folder + '/*', recursive=False)):
        if not os.path.isfile(item):
            folder_stats = process_folder(item)
            stats = stats + folder_stats

    generate_summary(stats, writer)

merge_fasta_files(summary_folder, 'summary.fasta')
print('Done!')