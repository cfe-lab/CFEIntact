
class UserError(ValueError):
    def __init__(self, message: str):
        self.message = message
