import sys
import os
import warnings
import csv
import copy
from pathlib import Path

from .config import *

class Error(Exception):
    pass

class InputError(Error):
    def __init__(self, message):
        self.message = message

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

def resolve_rpath(rpath):
    return str(Path(rpath).resolve())

class ArgumentsObject():
    def __init__(self, args_dict):

        self.common_dict = args_dict

        nxf_script = f"{args_dict['nextflow_modules_dir']}/{'_'.join([args_dict['task']] + ([args_dict['subtask']] if 'subtask' in args_dict else []))}.nf"
        self.nxf_script = nxf_script

        if 'config' in args_dict:
            self.config = args_dict['config']
        else:
            self.config = f"{args_dict['nextflow_modules_dir']}/{'_'.join([args_dict['task']] + ([args_dict['subtask']] if 'subtask' in args_dict else []))}.config"

        if 'profile' in args_dict:
            self.profile = args_dict['profile']
        else:
            self.profile = None

        if 'resume' in args_dict:
            self.resume = args_dict['resume']
        else:
            self.resume = None

        self.params = {}
        for k, v in args_dict.items():
            if k not in not_params_args:
                self.params[k] = v
        
        self.input_args = {}
        for k, v in args_dict.items():
            if k in run_args:
                self.input_args[k] = v

        self.file_input = args_dict['file_input'] if 'file_input' in args_dict else ''
        
    def add_param_to_params(self, key, value):
        if key in self.params:
            print(f"{key} already exists. It's value has been overwritten.", file=sys.stderr)

        self.params[key] = value
        return self

    def parse_file_input(self, new_obj):
        objs = []

        self.verify_file_input()

        with open(self.file_input, newline='') as csvfile:
            reader = file_handle_as_csv_dictreader_object(csvfile)
            for row in reader:
                for k, v in row.items():
                    new_obj.add_param_to_params(k, v)
            objs.append(new_obj)

        return objs

    def verify_file_input(self):
        with open(self.file_input, newline='') as csvfile:
            reader = file_handle_as_csv_dictreader_object(csvfile)
            header = reader.fieldnames

            for col_name in header:
                if col_name not in run_args:
                    raise InputError(f"File inputs can contain following arguments({'/'.join(run_args)})! The argument '{col_name}' != valid. Please write it in the command.")

            reader = file_handle_as_csv_dictreader_object(csvfile)
            header = reader.fieldnames

            try:
                for row in reader:
                    row = dict(filter(lambda item: item[1] != None, row.items()))

                    # check the length of each row
                    if len(header) != len(row):
                        raise InputError(f"Numbers of columns at header line and line {reader.line_num} in file input are different!")

                    # check if file exists or empty
                    for key in row.keys():
                        if key in ['x','x2']:
                            if not os.path.exists(row[key]) or not os.path.getsize(row[key]) > 0:
                                raise InputError(f"File {row[key]} does not exists or is empty.")

            except csv.Error as e:
                sys.exit(f"file {self.file_input}, line {reader.line_num}: {e}")
    
    def set_params_defaults(self):
        if 'outdir' not in self.params:
            self.params['outdir'] = self.params['prefix']

    def reproduce_run_objs(self):
        run_objs = []

        if self.file_input:
            run_obj = copy.deepcopy(self)
            run_objs += self.parse_file_input(run_obj)

        if not run_objs:
            if 'prefix' not in self.input_args:
                raise InputError("Either --file_input or --prefix is required!")

            n = len(self.input_args['prefix'])

            for k, v in self.input_args.items():
                assert n == len(v), f"Parameter '--{k}' has a different length from that of 'prefix'."

            for i in range(n):
                run_obj = copy.deepcopy(self)
            
                for k, v in self.input_args.items():
                    run_obj.add_param_to_params(k, v[i])

                run_objs.append(run_obj)

        for run_obj in run_objs:
            run_obj.verify_obj()
            run_obj.set_params_defaults()

        return run_objs

    def verify_obj(self):
        if "x" not in self.params:
            raise InputError("--x is required!")

        platform_not_specific = self.common_dict['task'] in ['polish', 'post_assembly'] or (self.common_dict['task'] == 'filter' and self.common_dict['subtask'] in ['map', 'blast', 'contigs'])
        # platform needed
        if not platform_not_specific:
            if "platform" not in self.params:
                raise InputError("--platform is required!")

        # Illumina paired end
        if not platform_not_specific and self.params['platform'] == "illumina":
                if "x2" not in self.params:
                    raise InputError("--x2 is required!")

        if self.common_dict['task'] == 'polish':
            if "reads" not in self.params:
                raise InputError("--reads is missing!")

        if platform_not_specific or self.params['platform'] == "nanopore":
            if "x2" in self.params:
                warnings.warn("Chosen analysis doesn't need --x2. Given argument will be ignored.")

        for k, v in self.params.items():
            # resolve relative path
            if k in path_args:
                if type(v) == list:
                    self.params[k] = [resolve_rpath(x) for x in v]
                else:
                    self.params[k] = resolve_rpath(v)