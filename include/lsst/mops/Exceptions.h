// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 


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





#endif


#endif
