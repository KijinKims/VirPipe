import sys
import subprocess
from virpipe.arguments_parser import Parser
from virpipe.cmd_generator import AnalysisCmdGenerator, BuildDBCmdGenerator
from virpipe.donwloader import DBDownloader

if __name__ == '__main__':
    parser = Parser(sys.argv[1:])
    args = parser.get_args()
    
    if args.task != 'database':
        cmd_generator = AnalysisCmdGenerator(args)
        
        for cmd in cmd_generator.get_cmds():
            print(cmd)
            subprocess.run(cmd, shell=True, check=True)
    else:
        if args.subtask == 'build':
            cmd_generator = BuildDBCmdGenerator(args)
            print(cmd_generator.get_cmd())
            subprocess.run(cmd_generator.get_cmd(), shell=True, check=True)
        else:
            downloader = DBDownloader(args)
            downloader.download()