import argparse
import pandas as pd
import os

class Entry:
    def __init__(self, seq_name, cord_start, cord_end, direction):
        self.name   =   seq_name

        if direction == '+': 
            self.start  =   cord_start
            self.end    =   cord_end
        else: #reverse
            self.start  =   cord_end
            self.end    =   cord_start
    
    def retrieve_as_pd_entry(self):
        return {    'name'  :   self.name,
                    'id'    :   self.name,
                    'start' :   self.start,
                    'end'   :   self.end
                }

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    return arg

parser = argparse.ArgumentParser(description= 'Script for generating metadata for zoontic rank from prodigal output file')

parser.add_argument('--input', '-i', metavar='prodigal.output',
        help='Prodigal output file',
        type=lambda x: is_valid_file(parser, x))

parser.add_argument('--output', '-o', metavar='zoontic_rank_metadata.csv',
        help='Output file name where result will be written')

args = parser.parse_args()

df = pd.DataFrame()
with open(args.input) as f:
    for line in f.readlines():
        if line.startswith('# Sequence Data:'):
            seq_name = line.split('seqhdr=')[1].strip().strip('"').split()[0]
        elif line.startswith('# Model Data:'):
            transl_table = line.split('transl_table=')[1].split(';')[0]
        elif line.startswith('>'):
            cds_list = line.strip().split('_')
            cord_start = cds_list[1]
            cord_end = cds_list[2]
            direction = cds_list[3]
            entry = Entry(seq_name, cord_start, cord_end, direction)
            if transl_table == '11':
                df = pd.concat([df, entry.retrieve_as_pd_entry()], ignore_index=True)

df = df.rename(columns={'name':'Name','id':'SequenceID','start':'CodingStart','end':'CodingStop'})
df.to_csv(args.output, index=False, header=True, columns=['Name','SequenceID','CodingStart','CodingStop'])
