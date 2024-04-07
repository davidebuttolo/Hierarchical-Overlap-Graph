#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <string>

#ifdef _WIN32
    #include "windows.h"
    #include "psapi.h"
#endif // _WIN32

void memory_usage_MB(std::string message);
bool ends_with(std::string const &a, std::string const &b);
void to_upper(std::string &s);
std::string trim(const std::string &in, const std::string &search_value);

#endif // UTIL_HPP