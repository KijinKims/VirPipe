import sys
import os
import warnings
import csv
from argparse import Namespace
from typing import Tuple, List, Dict, Any

from virpipe.defined_errors import InputError
from virpipe.arguments_loader import *
from virpipe.config import *

def file_handle_as_csv_dictreader_object(file_handle):
    _, file_extension = os.path.splitext(file_handle.name)

    if file_extension == '.csv':
        reader = csv.DictReader(file_handle)
    elif file_extension == '.tsv':
        reader = csv.DictReader(file_handle, delimiter='\t')
    elif file_extension == '.xlsx' or file_extension == '.xls':
        reader = csv.DictReader(file_handle, dialect='excel')
    else:
        warnings.warn("The file given to --file_input is neither csv, tsv nor xlsx. So it's format is deduced by csv::Sniffer. If not working, please consider changing the file extension in accordance with your file's format.")
        dialect = csv.Sniffer().sniff(file_handle.read(1024))
        file_handle.seek(0)
        reader = csv.DictReader(file_handle, dialect)

    return reader

class CmdGenerator:
    pass

class AnalysisCmdGenerator(CmdGenerator):
    def __init__(self, args):

        argsloader : ArgsLoader
        file_input : str
        
        argsloader, file_input = self.load_args_to_loader(args)

        if file_input != '':
            self.validify_file_input(file_input)
            argsloader += self.convert_file_input_to_argsloader(file_input)
        self.validify_argsloader(argsloader)
        argsloader = self.set_defaults(argsloader)

        self.argsloader = argsloader

    def load_args_to_loader(self, args : Namespace) -> Tuple[ArgsLoader, str]:
        print(args)
        argsloader = ArgsLoader()
        file_input = ''

        for attr in filter(lambda a: not a.startswith('_'), dir(args) ):
            if attr == 'file_input':
                file_input = getattr(args, attr)
            else: # other args shared to all samples
                if type(getattr(args, attr)) == list:
                    if attr =='ref' or attr =='cds':
                        argsloader.add(FileListArg(attr, getattr(args, attr)))
                    elif attr == 'tool':
                        argsloader.add(ValueArg(attr, ' '.join(getattr(args, attr))))
                    else:
                        argsloader.add(ListArg(attr, getattr(args, attr)))
                elif type(getattr(args, attr)) == str:
                    argsloader.add(ValueArg(attr, getattr(args, attr)))
                elif type(getattr(args, attr)) == int:
                    argsloader.add(ValueArg(attr, getattr(args, attr)))
                elif type(getattr(args, attr)) == bool:
                    argsloader.add(BooleanArg(attr, getattr(args, attr)))
                elif type(getattr(args, attr)) == dict:
                    argsloader.add(DictArg(attr, getattr(args, attr)))
                else:
                    print("This case cannot happen.", file=sys.stderr)
                    exit(1)

            if argsloader.has('ref'):
                argsloader.get('ref').resolve_relative_path()

            if argsloader.has('cds'):
                argsloader.get('cds').resolve_relative_path()

        return argsloader, file_input

    def validify_file(self, file_input):
        with open(file_input, newline='') as csvfile:
            reader = file_handle_as_csv_dictreader_object(csvfile)
            header = reader.fieldnames

            try:
                for row in reader:
                    row = dict(filter(lambda item: item[1] is not None, row.items()))

                    # check the length of each row
                    if len(header) != len(row):
                        raise InputError(f"Numbers of columns at header line and line {reader.line_num} in file input are different!")

                    # check if file exists or empty
                    for key in row.keys():
                        if key in file_args:
                            if not os.path.exists(row[key]) or not os.path.getsize(row[key]) > 0:
                                raise InputError(f"File {row[key]} does not exists or is empty.")

            except csv.Error as e:
                sys.exit(f'file {file_input}, line {reader.line_num}: {e}')

    def validify_file_input(self, file_input):
        with open(file_input, newline='') as csvfile:
            reader = file_handle_as_csv_dictreader_object(csvfile)
            header = reader.fieldnames

            for col_name in header:
                if col_name not in sample_specific_args:
                    raise InputError(f"File inputs should contain only sample-specific arguments({'/'.join(sample_specific_args)})! The argument '{col_name}' is not sample-speicific or invalid.")

        self.validify_file(file_input)

    def validify_file_ref_cds(self, file_ref_cds):
        with open(file_ref_cds, newline='') as csvfile:
            reader = file_handle_as_csv_dictreader_object(csvfile)
            header = reader.fieldnames

            for col_name in header:
                if col_name not in ['ref', 'cds']:
                    raise InputError(f"File ref_cds should contain only 'ref' and 'cds' columns! The argument '{col_name}' is neither of them.")

        self.validify_file(file_ref_cds)

    def convert_file_input_to_argsloader(self, file_input):
        argsloader = ArgsLoader()
        with open(file_input, newline='') as csvfile:
            reader = file_handle_as_csv_dictreader_object(csvfile)
            for row in reader:
                for key, value in row.items():
                    argsloader.add(ListArg(key, value))
        
        return argsloader

    def convert_ref_cds_to_argsloader(self, file_ref_cds):
        argsloader = ArgsLoader()
        with open(file_ref_cds, newline='') as csvfile:
            reader = file_handle_as_csv_dictreader_object(csvfile)
            for row in reader:
                for key, value in row.items():
                    argsloader.add(FileListArg(key, value))
        
        return argsloader
    
    def validify_argsloader(self, argsloader : ArgsLoader):

        # validation of parity of sample-specific files
        if argsloader.not_has('prefix'):
            raise InputError("--prefix is required!")

        prefix_num = len(argsloader['prefix'])
        
        sample_num = 0
        if argsloader['task'] != 'report':
            if argsloader.has('x'):
                sample_num = len(argsloader['x'])
            
            if prefix_num != sample_num:
                if sample_num == 0:
                    raise InputError("-x is required!")
                else:
                    raise InputError("--prefix and -x have different lengths!")

            if argsloader['task'] in ['polish', 'post_assembly'] or argsloader['platform'] == 'nanopore':
                if argsloader.has('x2'):
                    warnings.warn("Chosen analysis doesn't need -x2. Given argument will be ignored.")
                if argsloader['task'] == 'polish':
                    if argsloader.not_has('reads'):
                        raise InputError("--reads is missing!")
                    else:
                        if len(argsloader['x']) != len(argsloader['reads']):
                            raise InputError("-x and --reads have different lengths!")
            else:
                if argsloader.not_has('x2'):
                    raise InputError("-x2 is required!")
                else:
                    if len(argsloader['x']) != len(argsloader['x2']):
                        raise InputError("-x and -x2 have different lengths!")
                
                if argsloader['platform'] == 'hybrid':
                    if argsloader.not_has('y'):
                        raise InputError("-y is missing!")
                    else:    
                        if len(argsloader['x']) != len(argsloader['y']):
                            raise InputError("-x and -y have different lengths!")
                else:
                    if argsloader.has('y'):
                        warnings.warn("Chosen analysis doesn't need -y. Given argument will be ignored.")

        if argsloader.has('outdir'):
            if len(argsloader['outdir']) != prefix_num:
                raise InputError("--outdir is in an inappropriate length. It can either be in the same length with --prefix or not given.")

        if argsloader.has('host_genome'):
            if len(argsloader['host_genome']) != 1 and len(argsloader['host_genome']) != prefix_num:
                raise InputError("--host_genome is in an inappropriate length. It can be either in the same length with --prefix, 1 or not given.")

    def set_defaults(self, argsloader : ArgsLoader):

        if argsloader.not_has('outdir'):
            argsloader.add(ListArg('outdir', argsloader['prefix']))

        if argsloader['host_genome'] == 1:
            argsloader.add(ListArg('host_genome', argsloader['host_genome'] * len(argsloader['prefix'])))

        argsloader.add(ListArg('running_report', [ f"{x}/report.html" for x in argsloader['prefix'] ]))
        argsloader.add(ListArg('running_trace', [ f"{x}/trace.txt" for x in argsloader['prefix'] ]))
        argsloader.add(ListArg('running_timeline', [ f"{x}/timeline.html" for x in argsloader['prefix'] ]))

        return argsloader
    
    def get_cmds(self):
        return self.generate_nxf_cmds()

    def generate_nxf_cmds(self) -> str:
        
        nextflow_binary = self.argsloader['nextflow_binary']
        nextflow_modules_dir = self.argsloader['nextflow_modules_dir']
        task = self.argsloader['task']
        subtask = self.argsloader['subtask']

        nxf_cmd = f"{nextflow_binary} run"
        nxf_cmd += " "

        if task in ['post_assembly', 'filter', 'report']:
            nxf_cmd += f"{nextflow_modules_dir}/{task}/{subtask}/{task}_{subtask}.nf"
        else:
            nxf_cmd += f"{nextflow_modules_dir}/{task}/{task}.nf"

        cmds = [nxf_cmd] * len(self.argsloader['prefix'])
        for i, cmd in enumerate(cmds):
            new_cmd = cmd
            for arg_name in self.argsloader:

                if arg_name not in not_forwarded_to_nxf_args:
                    arg = self.argsloader.get(arg_name)

                    if type(arg) == ListArg:
                        new_cmd += ' ' + arg.nxf_cmd(i)
                
                    else:
                        if arg.nxf_cmd():
                            new_cmd += ' ' + arg.nxf_cmd()
            cmds[i] = new_cmd
        return cmds

class BuildDBCmdGenerator(CmdGenerator):
    def __init__(self, args):
        self.argsloader : ArgsLoader = self.load_args_to_loader(args)

    def load_args_to_loader(self, args : Namespace) -> ArgsLoader:
        argsloader = ArgsLoader()

        for attr in filter(lambda a: not a.startswith('_'), dir(args) ):
            if type(getattr(args, attr)) == str:
                argsloader.add(ValueArg(attr, getattr(args, attr)))
            elif type(getattr(args, attr)) == bool:
                argsloader.add(BooleanArg(attr, getattr(args, attr)))
            else:
                print("This case cannot happen.", file=sys.stderr)
                exit(1)

        return argsloader

    def get_cmd(self):

        return self.generate_nxf_cmd()

    def generate_nxf_cmd(self) -> str:
        nextflow_binary = self.argsloader['nextflow_binary']
        nextflow_modules_dir = self.argsloader['nextflow_modules_dir']
        program = self.argsloader['program']

        cmd = f"{nextflow_binary} run {nextflow_modules_dir}/database/build/{program}/database_{program}.nf"

        for arg_name in self.argsloader:

            if arg_name not in not_forwarded_to_nxf_args:
                arg = self.argsloader.get(arg_name)
                cmd += ' ' + arg.nxf_cmd()

        return cmd
