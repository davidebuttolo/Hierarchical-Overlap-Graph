#include "util.hpp"

void memory_usage_MB(std::string message)
{
    #ifdef _WIN32   

    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS *)&pmc, sizeof(pmc));
    SIZE_T virtualMemUsedByMe = pmc.PrivateUsage;
    SIZE_T physMemUsedByMe = pmc.WorkingSetSize;

    std::cout << message << physMemUsedByMe / (1024 * 1024) << " MB" << std::endl;

    #endif // _WIN32
}

bool ends_with(std::string const &a, std::string const &b)
{
    if (b.size() > a.size())
        return false;
    return std::equal(b.rbegin(), b.rend(), a.rbegin());
}

void to_upper(std::string &s)
{
    size_t length = s.size();
    for (size_t i = 0; i < length; i++)
    {
        s[i] -= 32 * (s[i] >= 'a' && s[i] <= 'z');
    }
}

std::string trim(const std::string &in, const std::string &search_value)
{
    size_t begin = in.find_first_not_of(search_value);
    size_t end = in.find_last_not_of(search_value);

    return in.substr(begin, end - begin + 1);
}