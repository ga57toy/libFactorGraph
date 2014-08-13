#ifndef FGCOMMON_H
#define FGCOMMON_H

#define _PI 3.14159265359

namespace FG
{
    enum NodeType{Variable, Factor};
    enum pdfType {Single, Gaussian, Binary, Beta}; //pdfType Single is actually an "Generate Distribution", see: http://en.wikipedia.org/wiki/Degenerate_distribution
}

#endif // FGCOMMON_H
