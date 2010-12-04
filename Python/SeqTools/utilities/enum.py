def Enum(*sequential, **named):
    """Simple enum implementation from stackoverflow.
    
    See post by Alec Thomas here:
    http://stackoverflow.com/questions/36932
    
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)
