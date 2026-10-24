#include "TextUtils.h"
#include <filesystem>
#include <stdexcept>
#include <vector>
#include <windows.h>

namespace text_utils
{
// Vstup a výstup sú v UTF-8; normalizácia vo Windows pracuje s UTF-16.
std::string withoutDiacritics(const std::string &text)
{
    if (text.empty())
        return text;
    const auto wideName = std::filesystem::u8path(text).wstring();
    int capacity = NormalizeString(NormalizationD, wideName.c_str(), -1, nullptr, 0);
    if (capacity <= 0)
        throw std::runtime_error("Cannot normalize UTF-8 text");
    std::wstring normalized;
    int length;
    for (;;)
    {
        normalized.resize(capacity);
        length = NormalizeString(NormalizationD, wideName.c_str(), -1, normalized.data(), capacity);
        if (length > 0)
            break;
        if (length < 0 && GetLastError() == ERROR_INSUFFICIENT_BUFFER)
            capacity = -length;
        else
            throw std::runtime_error("Cannot normalize UTF-8 text");
    }
    std::vector<WORD> types(length);
    if (!GetStringTypeW(CT_CTYPE3, normalized.data(), length, types.data()))
        throw std::runtime_error("Cannot classify characters in UTF-8 text");
    std::wstring plainName;
    for (int i = 0; i < length - 1; ++i) // Vynechá ukončovací nulový znak.
        if (!(types[i] & C3_NONSPACING))
            plainName += normalized[i];
    return std::filesystem::path(plainName).u8string();
}

}

