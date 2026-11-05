#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <fstream>
#include <iostream>
#include <memory>

// Fixture conversion is deliberately independent of GWLD's input adapter.
int main(int argc,char* argv[]) {
    if (argc!=4) return 1;
    if (std::string(argv[3])=="wz") {
        std::ifstream input(argv[1],std::ios::binary);
        BGZF* output=bgzf_open(argv[2],"w");
        if (!input || !output) return 1;
        char bytes[65536]; bool okay=true;
        while (input.read(bytes,sizeof(bytes)) || input.gcount())
            if (bgzf_write(output,bytes,input.gcount())!=input.gcount()) {okay=false; break;}
        if (bgzf_close(output)<0) okay=false;
        return okay ? 0 : 1;
    }
    htsFile* input=hts_open(argv[1],"r");
    if (!input) return 1;
    bcf_hdr_t* header=bcf_hdr_read(input);
    htsFile* output=hts_open(argv[2],argv[3]);
    bcf1_t* record=bcf_init();
    bool okay=header && output && record;
    if (okay && bcf_hdr_write(output,header)<0) okay=false;
    int status=-1;
    if (okay) while ((status=bcf_read(input,header,record))==0)
        if (record->errcode || bcf_write(output,header,record)<0) {okay=false; break;}
    if (status< -1) okay=false;
    if (record) bcf_destroy(record);
    if (header) bcf_hdr_destroy(header);
    if (hts_close(input)<0) okay=false;
    if (output && hts_close(output)<0) okay=false;
    return okay ? 0 : 1;
}
