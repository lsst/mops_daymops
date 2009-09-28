"""
DayMOPSStage

Abstract class that provides some basic functionality shared by all DayMOPS 
Stages.
"""
import lsst.pex.harness.Stage as Stage
import lsst.pex.policy as policy
import lsst.pex.logging as logging



class DayMOPSStage(object, Stage.Stage):
    def __init__(self, stageId=-1, policy=None):
        """
        Standard Stage initializer.
        """
        # First init the parent class.
        Stage.Stage.__init__(self, stageId, policy)
        
        # Make sure that we load all the nested policy files in out Stage 
        # policy.
        self._policy.loadPolicyFiles()
        
        # If we have a parameter which is itself a policy file, merge it with 
        # self._policy fr added convenience. We do not do this recursevly!!!!
        for policyName in self._policy.policyNames():
            p = self._policy.get(policyName)
            for key in p.names():
                self._policy.set(key, p.get(key))
        
        
        # Then setup logging.
        self._log = logging.Log(logging.Log.getDefaultLog(), 
                                self.__class__.__name__)
        if(isinstance(self._log, logging.ScreenLog)):
            self._log.setScreenVerbose(True)
        
        # Read some basic verbose levels from policy. Supported verbosityLevels
        # are: DEBUG, INFO, WARN, FATAL. Default is INFO
        userLevel = self._policy.get('verbosityLevel')
        if(userLevel and hasattr(logging.Log, userLevel)):
            self.verbosityLevel = getattr(logging.Log, userLevel)
            msg = 'Setting log verbosity to %s' %(userLevel)
        elif(not userLevel):
            self.verbosityLevel = logging.Log.DEBUG
            msg = 'Verbosity level not specified. Set to DEBUG.'
        else:
            self.verbosityLevel = logging.Log.DEBUG
            msg = 'Verbosity level "%s" not supported. Set to DEBUG.' \
                  %(userLevel)
        
        logging.Trace_setVerbosity('lsst.daymops', self.verbosityLevel)
        self._log.setThreshold(self.verbosityLevel)
        self.logIt('INFO', msg)
        return
    
    def getValueFromPolicy(self, paramName, defaultValue=None):
        """
        Fetch the given parameter from policy and assign it to the given default
        value. If no default value is not provided and the parameter is not 
        defined in the policy, then raise an exception.
        """
        val = self._policy.get(paramName)
        if(val == None and defaultValue != None):
            msg = '%s not found in the policy file. Using default value (%s).' \
                  %(paramName, str(defaultValue))
            self.logIt('INFO',  msg)
            return(defaultValue)
        elif(val == None):
            raise(Exception('%s not found in policy and no defaut specified.' \
                            %(paramName)))
        
        self.logIt('DEBUG', 'Policy: %s=%s' %(paramName, str(val)))
        return(val)
    
    def logIt(self, level, logString):
        """
        Write logString to self._log using the given verbosity level.
        
        @param logString: the message to write in the logs.
        @param level: verbosity level. Accepted values DEBUG, INFO, WARN, FATAL.
        """
        return(logging.Rec(self._log, getattr(logging.Log, level)) << logString << logging.endr)





