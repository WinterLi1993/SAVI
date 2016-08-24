CC=g++
CFLAGS=-O2

all: pileup2multiallele_vcf add_to_info make_qvt savi_poster savi_conf savi_comp savi_poster_accum savi_poster_merge savi_poster_prt savi_unif_prior savi_txt2prior filter_pileup

pileup2multiallele_vcf: pileup2multiallele_vcf.cpp
	$(CC) $(CFLAGS) pileup2multiallele_vcf.cpp -o pileup2multiallele_vcf

add_to_info: add_to_info.cpp
	$(CC) $(CFLAGS) add_to_info.cpp -o add_to_info

# for up to date versions of gcc, use: -std=c++11
# needed to use: http://www.cplusplus.com/reference/string/stoi/
make_qvt: make_qvt.cpp
	$(CC) $(CFLAGS) -std=c++0x make_qvt.cpp -o make_qvt

filter_pileup: filter_pileup.cpp
	$(CC) $(CFLAGS) -std=c++0x filter_pileup.cpp -o filter_pileup 
	
savi_poster: savi.cpp fasd/dcdflib.cpp fasd/dcdflib.h
	$(CC) $(CFLAGS) -Dmposter savi.cpp fasd/dcdflib.cpp -o savi_poster

savi_conf: savi.cpp fasd/dcdflib.cpp fasd/dcdflib.h
	$(CC) $(CFLAGS) -Dmconf savi.cpp fasd/dcdflib.cpp -o savi_conf

savi_comp: savi.cpp fasd/dcdflib.cpp fasd/dcdflib.h
	$(CC) $(CFLAGS) -Dmcomp savi.cpp fasd/dcdflib.cpp -o savi_comp

savi_poster_accum: savi.cpp fasd/dcdflib.cpp fasd/dcdflib.h
	$(CC) $(CFLAGS) -Dmposter_accum savi.cpp fasd/dcdflib.cpp -o savi_poster_accum

savi_poster_merge: savi.cpp fasd/dcdflib.cpp fasd/dcdflib.h
	$(CC) $(CFLAGS) -Dmposter_merge savi.cpp fasd/dcdflib.cpp -o savi_poster_merge

savi_poster_prt: savi.cpp fasd/dcdflib.cpp fasd/dcdflib.h
	$(CC) $(CFLAGS) -Dmposter_prt savi.cpp fasd/dcdflib.cpp -o savi_poster_prt 

savi_unif_prior: savi.cpp fasd/dcdflib.cpp fasd/dcdflib.h
	$(CC) $(CFLAGS) -Dmunif savi.cpp fasd/dcdflib.cpp -o savi_unif_prior

savi_txt2prior: savi.cpp fasd/dcdflib.cpp fasd/dcdflib.h
	$(CC) $(CFLAGS) -Dmtxt2prior savi.cpp fasd/dcdflib.cpp -o savi_txt2prior

clean:
	rm -f pileup2multiallele_vcf add_to_info make_qvt savi_poster savi_conf savi_comp savi_poster_accum savi_poster_merge savi_poster_prt savi_unif_prior savi_txt2prior filter_pileup
