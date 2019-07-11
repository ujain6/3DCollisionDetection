#ifndef _DEBUG_H_
#define _DEBUG_H_

#include <iostream>
#include <cstdlib>
#include <cstdarg>

template <typename T>
void _debugexp(const char* expr_str, T expr, const char* file, int line){
#ifdef DEBUG
    std::cout << "[" << file << ", line" << line << "] "; // << std::endl;
    std::cout << "DebugExpr: " << expr_str << " = " << expr << std::endl;
#endif
}

#define debugexp(Expr)     _debugexp(#Expr, Expr, __FUNCTION__, __LINE__)

#define  info(      Msg, ...)  _info(             __FUNCTION__, __LINE__, Msg, ##__VA_ARGS__)

#define debug(      Msg, ...) _debug(             __FUNCTION__, __LINE__, Msg, ##__VA_ARGS__)

#define  warn(Expr, Msg, ...)  _warn(#Expr, Expr, __FUNCTION__, __LINE__, Msg, ##__VA_ARGS__)

void _warn(const char* expr_str, bool expr, const char* file, int line, const char* msg, ...);

void _guard(const char* expr_str, bool expr, const char* file, int line, const char* msg, ...);

void _debug(const char* file, int line, const char* msg, ...);

void  _info(const char* file, int line, const char* msg, ...);

void die();

template <int N, typename T> 
void usage(int argc, T str)
{ 
    if (argc != N){
        std::cout << "Usage: " << str << std::endl; 
        die();
    }
    
}

void _warn(const char* expr_str, bool expr, const char* file, int line, const char* msg, ...)
{
#ifdef DEBUG
    if (!expr) { 
        /* Construct buffer from vargs */ 
        const int MAXBUF = 1024;
        char buffer[MAXBUF]; va_list args; va_start (args, msg); vsnprintf (buffer, MAXBUF-1, msg, args); va_end (args);
        
        std::cout << "[" << file << ", line" << line << "] " ; // << std::endl;
        std::cout << "Warn: " << expr_str << std::endl;
        std::cout << buffer << std::endl;
    }
#endif
}

void _debug(const char* file, int line, const char* msg, ...)
{
#ifdef DEBUG
    /* Construct buffer from vargs */ 
    const int MAXBUF = 1024;
    char buffer[MAXBUF]; va_list args; va_start (args, msg); vsnprintf (buffer, MAXBUF-1, msg, args); va_end (args);

    std::cout << "[" << file << ", line" << line << "] "; // << std::endl;
    std::cout << "Debug: " << buffer << std::endl;
#endif
}

void _info(const char* file, int line, const char* msg, ...)
{
    /* Construct buffer from vargs */ 
    const int MAXBUF = 1024;
    char buffer[MAXBUF]; va_list args; va_start (args, msg); vsnprintf (buffer, MAXBUF-1, msg, args); va_end (args);

    std::cout << "[" << file << ", line" << line << "] "; // << std::endl;
    std::cout << "Info: " << buffer << std::endl;
}


void die() 
{
    exit(EXIT_FAILURE);
}

#endif // _DEBUG_H_