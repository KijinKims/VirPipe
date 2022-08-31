import sys
import subprocess
from virpipe.arguments_parser import Parser
from virpipe.cmd_generator import AnalysisCmdGenerator

if __name__ == '__main__':
    parser = Parser(sys.argv[1:])
    args = parser.get_args()
    

    cmd_generator = AnalysisCmdGenerator(args)

    for cmd in cmd_generator.get_cmds():
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)