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

namespace collapseTracklets {
    namespace exceptions {
       

#define LSST_EXCEPT(class, str)                 \
        class(str);
        
#define LSST_EXCEPTION_TYPE(t, b, c)                    \
        class t : public b {                            \
        public:                                         \
        t() : b () { };                                 \
        t(std::string s) { std::cerr << s; exit(1);  }; \
        };
        
        LSST_EXCEPTION_TYPE(ProgrammerErrorException, std::exception, 
                            collapseTracklets::exceptions::ProgrammerErrorException)
        
        LSST_EXCEPTION_TYPE(BadParameterException, ProgrammerErrorException, 
                            collapseTracklets::exceptions::BadParameterException)
        
        LSST_EXCEPTION_TYPE(UserErrorException, std::exception,
                            collapseTracklets::exceptions::UserErrorException)
        
        LSST_EXCEPTION_TYPE(InputFileFormatErrorException, UserErrorException, 
                            collapseTracklets::exceptions::InputFileFormatErrorException)
        
        LSST_EXCEPTION_TYPE(CommandlineParseErrorException, UserErrorException, 
                            collapseTracklets::exceptions::CommandlineParseErrorException)
        
        LSST_EXCEPTION_TYPE(GSLException, std::exception, 
                            collapseTracklets::exceptions::GSLException)
        
        LSST_EXCEPTION_TYPE(FileException, std::exception,
                            collapseTracklets::exceptions::FileException)

        LSST_EXCEPTION_TYPE(MemoryException, std::exception,
                            collapseTracklets::exceptions::MemoryException)

        LSST_EXCEPTION_TYPE(UninitializedException, std::exception,
                            collapseTracklets::exceptions::UninitializedException)
        
    
    }
}

#else
// this is not for PANSTARRS - use real LSST pex exceptions


#include <string>
#include "lsst/pex/exceptions/Exception.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace collapseTracklets {
namespace exceptions {

LSST_EXCEPTION_TYPE(ProgrammerErrorException, lsst::pex::exceptions::Exception, 
                    collapseTracklets::exceptions::ProgrammerErrorException)

LSST_EXCEPTION_TYPE(BadParameterException, ProgrammerErrorException, 
                    collapseTracklets::exceptions::BadParameterException)

LSST_EXCEPTION_TYPE(UserErrorException, pexExcept::Exception,
                    collapseTracklets::exceptions::UserErrorException)

LSST_EXCEPTION_TYPE(InputFileFormatErrorException, UserErrorException, 
                    collapseTracklets::exceptions::InputFileFormatErrorException)

LSST_EXCEPTION_TYPE(CommandlineParseErrorException, UserErrorException, 
                    collapseTracklets::exceptions::CommandlineParseErrorException)

LSST_EXCEPTION_TYPE(GSLException, lsst::pex::exceptions::Exception, 
                    collapseTracklets::exceptions::GSLException)

LSST_EXCEPTION_TYPE(FileException, lsst::pex::exceptions::Exception,
                    collapseTracklets::exceptions::FileException)

LSST_EXCEPTION_TYPE(UninitializedException, lsst::pex::exceptions::Exception,
                    collapseTracklets::exceptions::UninitializedException)


#define CALL_AND_CATCH_EXCEPTIONS(INSTRUCTIONS)                         \
    try {                                                               \
        INSTRUCTIONS;                                                   \
    }                                                                   \
    catch (collapseTracklets::exceptions::BadParameterException bpe) {  \
        std::cerr << "Encountered a programmer error "                  \
                  << "due to bad parameters to a function. "            \
                  << std::endl;                                         \
        std::cerr << "Message: " << std::endl << std::endl;             \
        std::cerr << bpe.what() << std::endl;                           \
        exit(1);                                                        \
    }                                                                   \
    catch (collapseTracklets::exceptions::ProgrammerErrorException      \
           const& pe) {                                                 \
        std::cerr << "Encountered a programmer error "                  \
                  << std::endl;                                         \
        std::cerr << "Message: " << std::endl << std::endl;             \
        std::cerr << pe.what() << std::endl;                            \
        exit(2);                                                        \
    }                                                                   \
    catch (collapseTracklets::exceptions::UserErrorException            \
           const& ue) {                                                 \
        std::cerr << ue.what() << std::endl;                            \
        exit(3);                                                        \
    }                                                                   \
    catch (collapseTracklets::exceptions::GSLException const& gsle) {   \
        std::cerr << "Encountered a problem with GSL." << std::endl;    \
        std::cerr << gsle.what() << std::endl;                          \
        exit(4);                                                        \
    }                                                                   \
    catch (collapseTracklets::exceptions::FileException const& fe) {    \
        std::cerr << fe.what() << std::endl;                            \
        exit(5);                                                        \
    }                                                                   \
    catch (lsst::pex::exceptions::Exception const& e) {                 \
        std::cerr << "Unexpected pex_exception:\n\n";                  \
    std::cerr << e.what() << std::endl;                                 \
    exit(6);                                                            \
    }                                                                   \


}}

#endif
#endif
