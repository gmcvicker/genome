
import sys

class DBAdaptor(object):
    """Base class for database adaptors."""
    
    def __init__(self, track):
        """Initializes a database adaptor with an open database track"""
        self.track = track
