import sys
from .arguments.parser import Parser
from .arguments.validation import Validator
from .command.nextflow import CommandConverter
from .command.run_docker import DockerRunner

if __name__ == '__main__':
    args = Parser(sys.argv[1:]).get_args()

    args = Validator(args).validate()

    nextflow_command = CommandConverter(args).generate_nextflow_script_command()

    DockerRunner(args).run(nextflow_command)