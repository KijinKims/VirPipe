import sys
from .arguments.parser import Parser
from .arguments.validation import Validator
from .command.nextflow import CommandConverter
from .command.run_docker import DockerRunner

def main():
    p = Parser(sys.argv[1:])
    
    args = p.get_args()

    args = Validator(args).validate()

    nextflow_command = CommandConverter(args).generate_nextflow_script_command()

    DockerRunner(args).run(nextflow_command)