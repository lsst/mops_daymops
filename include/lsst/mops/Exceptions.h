// -*- LSST-C++ -*-


/* jmyers 1/15/09
 *
 * This is a really basic set of exceptions to be used until pex_exceptions
 * are sufficiently stable; no real attempt is made to build these with 
 * forward-compatible call semantics or anything, so use caution.
 *
 */


#ifndef LSST_PANSTARRS_INTERIM_EXCEPTION
#define LSST_PANSTARRS_INTERIM_EXCEPTION




#ifdef PANSTARRS
// use the really naive exceptions, not the LSST pex exceptions
#include <exception>
#include <string>
#include <iostream>

       

#define LSST_EXCEPT(class, str)                 \
        class(str);
        
#define LSST_EXCEPTION_TYPE(t, b, c)                    \
        class t : public b {                            \
        public:                                         \
        t() : b () { };                                 \
        t(std::string s) { std::cerr << s; exit(1);  }; \
        };
        
        LSST_EXCEPTION_TYPE(ProgrammerErrorException, std::exception, 
                            lsst::mops::ProgrammerErrorException)
        
        LSST_EXCEPTION_TYPE(BadParameterException, ProgrammerErrorException, 
                            lsst::mops::BadParameterException)
        
        LSST_EXCEPTION_TYPE(UserErrorException, std::exception,
                            lsst::mops::UserErrorException)
        
        LSST_EXCEPTION_TYPE(InputFileFormatErrorException, UserErrorException, 
                            lsst::mops::InputFileFormatErrorException)
        
        LSST_EXCEPTION_TYPE(CommandlineParseErrorException, UserErrorException, 
                            lsst::mops::CommandlineParseErrorException)
        
        LSST_EXCEPTION_TYPE(GSLException, std::exception, 
                            lsst::mops::GSLException)
        
        LSST_EXCEPTION_TYPE(FileException, std::exception,
                            lsst::mops::FileException)

        LSST_EXCEPTION_TYPE(MemoryException, std::exception,
                            lsst::mops::MemoryException)

        LSST_EXCEPTION_TYPE(UninitializedException, std::exception,
                            lsst::mops::UninitializedException)

        LSST_EXCEPTION_TYPE(BadIndexException, std::exception,
                            lsst::mops::BadIndexException)
        
        LSST_EXCEPTION_TYPE(KnownShortcomingException, std::exception,
                            lsst::mops::KnownShortcomingException)
    
    


#else
// this is not for PANSTARRS - use real LSST pex exceptions


#include <string>
#include "lsst/pex/exceptions/Exception.h"
#include "lsst/pex/exceptions/Runtime.h"


LSST_EXCEPTION_TYPE(ProgrammerErrorException, lsst::pex::exceptions::Exception, 
                    lsst::mops::ProgrammerErrorException)

LSST_EXCEPTION_TYPE(BadParameterException, ProgrammerErrorException, 
                    lsst::mops::BadParameterException)

LSST_EXCEPTION_TYPE(UserErrorException, pexExcept::Exception,
                    lsst::mops::UserErrorException)

LSST_EXCEPTION_TYPE(InputFileFormatErrorException, UserErrorException, 
                    lsst::mops::InputFileFormatErrorException)

LSST_EXCEPTION_TYPE(CommandlineParseErrorException, UserErrorException, 
                    lsst::mops::CommandlineParseErrorException)

LSST_EXCEPTION_TYPE(GSLException, lsst::pex::exceptions::Exception, 
                    lsst::mops::GSLException)

LSST_EXCEPTION_TYPE(FileException, lsst::pex::exceptions::Exception,
                    lsst::mops::FileException)

LSST_EXCEPTION_TYPE(UninitializedException, lsst::pex::exceptions::Exception,
                    lsst::mops::UninitializedException)

LSST_EXCEPTION_TYPE(BadIndexException, lsst::pex::exceptions::Exception,
                    lsst::mops::BadIndexException)

LSST_EXCEPTION_TYPE(KnownShortcomingException, lsst::pex::exceptions::Exception,
                    lsst::mops::KnownShortcomingException)

LSST_EXCEPTION_TYPE(MemoryException, lsst::pex::exceptions::Exception,
                    lsst::mops::MemoryException)




#endif


#endif
