#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

const char        TAB               {  '\t' };

const std::string SEQUENCE_LINE     { "@SQ" };
const std::string SEQUENCE_NAME     { "SN:" };
const std::string SIZE_NAME         { "LN:" };
const char        BARCODE_FLAG[]    {  "BX" };
const size_t      BARCODE_FLAG_SIZE {    6  };
const size_t      LINE_LEN          {   60  };

#endif
