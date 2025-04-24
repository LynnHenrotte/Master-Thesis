class PharmVarApiException(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class InvalidArgumentError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class NoDataFoundError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class ValidationError(Exception):
    """Base class for validation errors"""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)