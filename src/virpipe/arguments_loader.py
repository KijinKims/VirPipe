from typing import List, Dict, Any
import warnings
from pathlib import Path

class Arg:
    def __init__(self, name):
        self.name = name

class ListArg(Arg):
    def __init__(self, name, values):
        super().__init__(name)
        if type(values) == list:
            self.li = values
        else:
            self.li = [values]

    def __len__(self):
        return len(self.li)

    def __add__(self, other_list_arg):
        self.li = self.li + other_list_arg.li
        return self

    def __getitem__(self, idx : int):
        return self.li[idx]

    def set(self, li):
        self.li = li

    def get(self):
        return self.li

    def add(self, value):
        self.li.append(value)

    def nxf_cmd(self, idx):
        d = {   'running_report' : '-with-report',
                'running_trace' : '-with-trace',
                'running_timeline' : '-with-timeline',
        }

        if self.name in d.keys():
            return f"{d[self.name]} {self.li[idx]}"
        else:
            return f"--{self.name} {self.li[idx]}"
            


        

class FileListArg(ListArg):
    def resolve_relative_path(self):
        self.li = [str(Path(p).resolve()) for p in self.li]

    def nxf_cmd(self):
        return f"--{self.name} \"{' '.join(self.li)}\""

class ValueArg(Arg):
    def __init__(self, name, value):
        super().__init__(name)
        self.val = value

    def set(self, value):
        self.val = value

    def get(self):
        return self.val    

    def nxf_cmd(self):
        d = {   'nextflow_config' : '-c',
        }

        if self.name in d.keys():
            return f"{d[self.name]} \"{self.val}\""
        else:
            return f"--{self.name} \"{self.val}\""
            

class BooleanArg(Arg):
    def __init__(self, name, val):
        super().__init__(name)
        self.val = val

    def set(self, value):
        self.val = value

    def get(self):
        return self.val

    def nxf_cmd(self):
        if self.val == True:
            d = {   'resume' : '-resume',
                    'background' : '-bg',
            }
            return d.get(self.name, f"--{self.name}")
        else:
            return None

class DictArg(Arg):
    def __init__(self, name, dict_):
        super().__init__(name)
        self._dict = dict_

    def set(self, dict_):
        self._dict = dict_

    def get(self):
        return self._dict

class ArgsLoader:
    def __init__(self):
        self.container : Dict[str, Arg] = {}

    def add(self, arg : Arg):
        if arg.name in self.container:
            if issubclass(type(arg),ListArg) :
                self.container[arg.name] += arg
            else:
                warnings.warn("Non-list type argument cannot be added. Nothing changes.")
        else:
            self.container[arg.name] = arg

    def __add__(self, other_args_loader):
        for name in other_args_loader.names():
            self.add(other_args_loader.get(name))

        return self

    def names(self) -> List[str]:
        return list(self.container.keys())

    def delete(self, arg_name : str):
        if arg_name in self.container:
            del self.container[arg_name]

    def get(self, arg_name : str) -> Arg:
        if arg_name in self.container:
            return self.container[arg_name]
        else:
            return None

    def __getitem__(self, arg_name : str):
        if arg_name in self.container:
            return self.container[arg_name].get()
        else:
            return None

    def __iter__(self):
        return iter(self.container)

    def has(self, arg_name : str) -> bool :
        if arg_name in self.names():
            return True
        else:
            return False
    
    def not_has(self, arg_name : str) -> bool :
        if arg_name not in self.names():
            return True
        else:
            return False