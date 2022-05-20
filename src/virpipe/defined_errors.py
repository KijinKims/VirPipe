class Error(Exception):
    pass

class InputError(Error):
    def __init__(self, message):
        self.message = message