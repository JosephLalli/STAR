#ifndef H_PersonalizedInput_helpers
#define H_PersonalizedInput_helpers
#include "IncludeDefine.h"
#include <cstring>
#include <map>
#include <zlib.h>
#include <stdexcept>
namespace PersonalizedInputHelpers {
vector<string> splitTabs(const string &line);
string personalizedOutputPath(const string &prefix, const string &sample, const string &kind);
void mergeCounts(map<string,uint64> &destination, const map<string,uint64> &source);
// Preserve row bytes, including CR, while accepting plain files, gzip/BGZF and FIFOs.
class LineReader {
    gzFile file;
    vector<char> buffer;
    size_t begin=0, end=0;
public:
    explicit LineReader(const string &path) : file(gzopen(path.c_str(),"rb")), buffer(1<<20) {
        if (file==NULL) throw runtime_error("cannot open input: " + path);
        gzbuffer(file,1<<20);
    }
    ~LineReader() { gzclose(file); }
    LineReader(const LineReader &) = delete;
    LineReader &operator=(const LineReader &) = delete;
    bool next(string &line) {
        line.clear();
        for (;;) {
            if (begin==end) {
                const int length=gzread(file,buffer.data(),buffer.size());
                int code;
                const char *message=gzerror(file,&code);
                if (length<0 || (length==0 && code!=Z_OK && code!=Z_STREAM_END))
                    throw runtime_error(string("input decompression/read failed: ")+message);
                if (length==0) return !line.empty();
                begin=0;
                end=length;
            };
            const char *first=buffer.data()+begin;
            const char *newline=static_cast<const char*>(memchr(first,'\n',end-begin));
            const size_t length=newline ? newline-first : end-begin;
            line.append(first,length);
            begin+=length;
            if (newline!=NULL) { ++begin; return true; };
        }
    }
};
}
#endif
