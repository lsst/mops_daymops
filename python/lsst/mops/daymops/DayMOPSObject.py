"""
Convenience super-class that automatically adds getters and setters according to
LSST desiderata.
"""
class DayMOPSObject(object):
    """
    Auto-create getters and setters.
    """
    def __getattr__(self, name):
        if name[:3] == "get":
            def getter(obj, key='_%s%s' %(name[3].lower(), name[4:])):
                return getattr(obj, key)
            setattr(self.__class__, name, getter)
            return(getattr(self, name))
        elif name[:3] == "set":
            def setter(obj, value, key='_%s%s' %(name[3].lower(), name[4:])):
                return(setattr(obj, key, value))
            setattr(self.__class__, name, setter)
            return(getattr(self, name))
        raise(AttributeError(name))
