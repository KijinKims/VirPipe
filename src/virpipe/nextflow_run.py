from virpipe.arguments_object import ArgumentsObject

class Pipeline():
    def __init__(self, run_obj : ArgumentsObject):
        self.run_obj = run_obj

    def generate_command(self):
        for param in self.run_obj.params:
            if type(param) == list:
                pass
            elif type(param) == bool:
                pass
            else:
                pass

    def run():
        pass