//
// Created by vsevolod on 19/12/22.
//

#ifndef SRC_LOGGER_H
#define SRC_LOGGER_H

//
// Created by vsevolod on 9/29/21.
//

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)

#define LOG_SILENT 0
#define LOG_ERR 1
#define LOG_WARN 2
#define LOG_INFO 3
#define LOG_TIME 4
#define LOG_DEBUG 5
#define LOG_DEFAULT 4

//
//
//#define LOG_INIT_COUT() logger Log(std::cout, __PRETTY_FUNCTION__)
//#define LOG_INIT_CERR() logger Log(std::cerr, __PRETTY_FUNCTION__)
//#define LOG_INIT_CLOG() logger Log(std::clog, __PRETTY_FUNCTION__)
//#define LOG_INIT_CUSTOM(X) logger Log((X), __PRETTY_FUNCTION__)
//
//#ifdef FORMAT_LOG_NO_COLORS
//
//#define FORMAT_LOG_TIME "[ TIME    ]"
//#define FORMAT_LOG_DEBUG "[ DEBUG   ]"
//#define FORMAT_LOG_ERROR "[ ERROR   ]"
//#define FORMAT_LOG_WARNING "[ WARNING ]"
//#define FORMAT_LOG_INFO "[ INFO    ]"
//
//#else
//
#define FORMAT_LOG_TIME    "\033[0;35m[ TIME    ]\033[0;0m"
#define FORMAT_LOG_DEBUG   "[ DEBUG   ]"
#define FORMAT_LOG_ERROR   "\033[0;31m[ ERROR   ]\033[0;0m"
#define FORMAT_LOG_WARNING "\033[0;33m[ WARNING ]\033[0;0m"
#define FORMAT_LOG_INFO    "\033[0;34m[ INFO    ]\033[0;0m"
//
//#endif

int CurrLogLevel = LOG_DEBUG;


class logger{
private:
    unsigned _message_level;
    unsigned _loglevel = CurrLogLevel;
    std::ostream& _fac;
    std::ostream& _err;
    std::string _name;

    static std::vector<std::string> split_location( const char *location ) {
        std::string tmp ( location );
        std::replace(tmp.begin(), tmp.end(), '/', ' ');
        std::vector<std::string> elements;
        std::stringstream ss(tmp);
        std::string temp;
        while (ss >> temp)
            elements.push_back(temp);
        return elements;
    }

    static std::string getFname( const char *location ){
        auto fpath = split_location(location);
        return fpath[fpath.size()-1];
    }

public:

    inline logger(std::ostream& f, std::ostream& ferr, unsigned ll, std::string n)
            : _message_level(LOG_SILENT), _fac(f), _err(ferr), _name( std::move(n) ) {
//        time(&_now);
//        time(&_start);

        _loglevel = ll;
    }

    inline logger(std::ostream& f, std::ostream& ferr, std::string n)
            : _message_level(LOG_SILENT), _fac(f), _err(ferr), _name( std::move(n) ) {
//        time(&_now);
//        time(&_start);
    }

//    inline void set_log_level( unsigned ll ) { _loglevel = ll; }

//    static unsigned& _loglevel() {
//        static unsigned _ll_internal = LOG_DEFAULT;
//        return _ll_internal;
//    };

    template <typename T>
    friend logger& operator<<(logger& l, const T& s );

    inline logger& operator()(unsigned ll, const char *location ){
        _message_level = ll;
//        if ( _loglevel == LOG_DEBUG ){
//            if (_message_level == LOG_ERR)
//                _err << " " << location << "\n" << "            ";
//            else
//                _fac << " " << location << "\n" << "            ";
//        } else {
//            if (_message_level <= _loglevel) {
//                if (_message_level == LOG_INFO) {
//                    _fac << FORMAT_LOG_INFO << " : ";
//                    _fac << "[ " << getFname(location) << " ] : ";
//                }
//                if (_message_level == LOG_WARN) {
//                    _fac << FORMAT_LOG_WARNING << " : ";
//                    _fac << "[ " << getFname(location) << " ] : ";
//                }
//                if (_message_level == LOG_ERR) {
//                    _err << FORMAT_LOG_ERROR << " : ";
//                    _err << "[ " << getFname(location) << " ] : ";
//                }
//
//                ": "; // prep_level(*this) << prep_time(*this) << prep_name(*this) <<
//            }
//        }

        if (_message_level <= _loglevel) {
            if (_message_level == LOG_INFO) {
                _fac << FORMAT_LOG_INFO << " : ";
                _fac << "[ " << getFname(location) << " ] : ";
            }
            if (_message_level == LOG_WARN) {
                _fac << FORMAT_LOG_WARNING << " : ";
                _fac << "[ " << getFname(location) << " ] : ";
            }
            if (_message_level == LOG_ERR) {
                _err << FORMAT_LOG_ERROR << " : ";
                _err << "[ " << getFname(location) << " ] : ";
            }

            ": "; // prep_level(*this) << prep_time(*this) << prep_name(*this) <<
        }

        return *this;
    }

    inline unsigned getLogLevel() const {return _loglevel;}
};
// FRIEND
template <typename T>
logger& operator<<(logger& l, const T& s ) {
    if ( l._message_level <= l._loglevel ) {
        l._fac << s;
        return l;
    } else {
        return l;
    }
}

#endif //SRC_LOGGER_H
