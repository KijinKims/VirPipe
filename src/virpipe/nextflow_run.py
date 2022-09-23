import subprocess
from .arguments_object import ArgumentsObject

class Pipeline():
    def __init__(self, run_obj : ArgumentsObject):
        self.run_obj = run_obj

    def generate_command(self):
        cmd = f"nextflow {self.run_obj.nxf_script} "

        cmd += f"-c {self.run_obj.config} "

        if self.run_obj.profile:
            if type(self.run_obj.profile) == list:
                cmd += f"-profile {','.join(self.run_obj.profile)} "
            else:
                cmd += f"-profile {self.run_obj.profile}" 

        for param, value in self.run_obj.params.items():
            if type(value) == list:
                cmd += f"--{param} {' '.join(value)} "
            elif type(value) == bool:
                if value:
                    cmd += f"--{param} "
            else:
                cmd += f"--{param} {value} "

        return cmd


    def run(self):
        command_string = self.generate_command()
        print(command_string)
        #subprocess.run(command_string, shell=True, check=True)