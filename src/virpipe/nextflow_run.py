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
                cmd += f"-profile {self.run_obj.profile} "

        if self.run_obj.resume:
            cmd += f"-resume "

        for param, value in self.run_obj.params.items():
            if param in ["with_report", "with_trace", "with_timeline"]:
                cmd += f"-{param.replace('_','-')} {self.run_obj.params['outdir']}/{value} "
            elif type(value) == list:
                cmd += f"--{param} {' '.join(value)} "
            elif type(value) == bool:
                if value:
                    cmd += f"--{param} "
            else:
                cmd += f"--{param} {value} "

        return cmd

    def generate_summary_command(self):
        cmd = f"Rscript {self.run_obj.summary_script} --markdown {self.run_obj.summary_markdown} --outdir {self.run_obj.params['outdir']} --prefix {self.run_obj.params['prefix']}"
        return cmd

    def run(self):
        command_string = self.generate_command()
        #print(command_string)
        subprocess.run(command_string, shell=True, check=True)

    def summary_run(self):
        command_string = self.generate_summary_command()
        subprocess.run(command_string, shell=True, check=True)
