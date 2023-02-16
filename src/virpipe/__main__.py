import sys

from .arguments_parser import Parser
from .arguments_object import ArgumentsObject
from .nextflow_run import Pipeline

if __name__ == '__main__':
    parser = Parser(sys.argv[1:])
    args = parser.get_args()

    arg_objs = ArgumentsObject(args)
    run_objs = arg_objs.reproduce_run_objs()

    for run_obj in run_objs:
        pipeline = Pipeline(run_obj)
        pipeline.run()