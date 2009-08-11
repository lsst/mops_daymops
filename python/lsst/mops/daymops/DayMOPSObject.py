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
    
    # Implement simple default getId() method.
    def getId(self):
        # Try and understand which attribute is an id. Usually it is an 
        # attribute of the form <className>[0].lower() + <className>[1:] + 'Id'
        # so the method we want is 'get' + <className> + 'Id'()
        methodName = 'get%sId' %(self.__class__.__name__)
        if(hasattr(self, methodName)):
            return(getattr(self, methodName)())
        return(None)
    
    # Default comparison is by id.
    def __lt__(self, other):
        if(other == None):
            return(False)
        return(self.getId() < other.getId())
    
    def __le__(self, other):
        if(other == None):
            return(False)
        return(self.getId() <= other.getId())
    
    def __eq__(self, other):
        if(other == None):
            return(False)
        return(self.getId() == other.getId())
    
    def __ne__(self, other):
        if(other == None):
            return(True)
        return(self.getId() != other.getId())

    def __gt__(self, other):
        if(other == None):
            return(False)
        return(self.getId() > other.getId())
    
    def __ge__(self, other):
        if(other == None):
            return(False)
        return(self.getId() >= other.getId())
