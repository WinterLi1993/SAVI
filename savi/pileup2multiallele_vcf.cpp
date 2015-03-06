#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cctype>

using namespace std;

// global variables
bool verbose = false;		// bool verbose
bool debug = false;		// bool debug
bool allvar = false;		// bool - if true, print all lines not just variants
bool bool_pval = false;		// bool - strand bias p val
int qual_offset = 33; 		// quality offset

// statistical functions taken verbatim from here:
// http://genome.sph.umich.edu/w/images/d/d8/Bios615-fa12-lec03-handout.pdf
// http://genome.sph.umich.edu/wiki/Biostatistics_615/815
/*
int fac(int n)  // calculates factorial
{
	int ret;
	for(ret=1; n > 0; --n) { ret *= n; }
	return ret;
}

double hypergeometricProb(int a, int b, int c, int d) 
{
	int num = fac(a+b) * fac(c+d) * fac(a+c) * fac(b+d);
	int den = fac(a) * fac(b) * fac(c) * fac(d) * fac(a+b+c+d);
	return (double)num/(double)den;
}
*/

double fac(int n) // main() function remains the same
{
	double ret; // use double instead of int
	for(ret=1.; n > 0; --n) { ret *= n; }
	return ret;
}

double hypergeometricProb(int a, int b, int c, int d) 
{
	double num = fac(a+b) * fac(c+d) * fac(a+c) * fac(b+d);
	double den = fac(a) * fac(b) * fac(c) * fac(d) * fac(a+b+c+d);
	return num/den; // use double to calculate factorials
}

double logFac(int n) 
{
	double ret;
	for(ret=0.; n > 0; --n) { ret += log((double)n); }
	return ret;
}

double logHypergeometricProb(int a, int b, int c, int d) 
{
	return logFac(a+b) + logFac(c+d) + logFac(a+c) + logFac(b+d) - logFac(a) - logFac(b) - logFac(c) - logFac(d) - logFac(a+b+c+d);
}

void initLogFacs(double* logFacs, int n) 
{
	logFacs[0] = 0;
	for(int i=1; i < n+1; ++i) 
	{
		logFacs[i] = logFacs[i-1] + log((double)i); // only n times of log() calls
	}
}

double logHypergeometricProb(double* logFacs, int a, int b, int c, int d) 
{
	return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d] - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}

double fastFishersExactTest(int a, int b, int c, int d) 
{
	// int a = atoi(argv[1]), b = atoi(argv[2]), c = atoi(argv[3]), d = atoi(argv[4]);
	int n = a + b + c + d;
	double* logFacs = new double[n+1]; // *** dynamically allocate memory logFacs[0..n] ***
	initLogFacs(logFacs, n); // *** initialize logFacs array ***
	double logpCutoff = logHypergeometricProb(logFacs,a,b,c,d); // *** logFacs added
	double pFraction = 0;
	for(int x=0; x <= n; ++x)
	{
		if ( a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0 ) 
		{
			double l = logHypergeometricProb(logFacs,x,a+b-x,a+c-x,d-a+x);
			if ( l <= logpCutoff ) pFraction += exp(l - logpCutoff);
		}
	}
	double logpValue = logpCutoff + log(pFraction);
	// std::cout << "Two-sided log10-p-value is " << logpValue/log(10.) << std::endl;
	// std::cout << "Two-sided p-value is " << exp(logpValue) << std::endl;
	delete [] logFacs;
	return exp(logpValue);
	// return 0;
}
// http://genome.sph.umich.edu/w/images/d/d8/Bios615-fa12-lec03-handout.pdf

// a triple for forward counts, reverse counts, quality
struct triple 
{
	int f;		// forward 
	int r;		// reverse
	int q;		// qual 
	// bool isdel;	// bool for "is deletion" 

	void setzero()
	{
		f = 0;
		r = 0;
		q = 0; 
	//	isdel = 0; 
	}
};

// sample ref counts and quals 
struct counts
{
	// ref counts
	int sdp; // raw depths
	int rdf; // reverse reads
	int rdr; // forward reads
	int rbq; // read base qual

	void setzero()
	{
		sdp = 0;
		rdf = 0;
		rdr = 0;
		rbq = 0;
	}
};

// get qual score of ASCII char give the quality score offset
int ascii2qual(char x, int qual_offset)
{
	// see http://en.wikipedia.org/wiki/Pileup_format:
	// "This is an optional column. If present, the ASCII value of the character minus 33 gives the mapping Phred quality of each of the bases in the previous column 5. This is similar to quality encoding in the FASTQ format"
	// qual_offset = 33;
	return int(x) - qual_offset;
}

// http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
	std::stringstream ss(s);
	std::string item;

	while (std::getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	return elems;
}

// http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
std::vector<std::string> split(const std::string &s, char delim)
{
	// Example usage:
	/*
	string x = "awef:dfwe:adsf:";
	vector<string> vec_format_field;
	char delim = ':';
	vec_format_field = split(x,delim);
	for (int i = 0; i < vec_format_field.size(); i++)	
	{	
		cout << vec_format_field[i] << endl;
	}
	*/

	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

// read stuff from std:in and re-echo it
void read_from_std()
{	
	// http://stackoverflow.com/questions/201992/how-to-read-until-eof-from-cin-in-c
	std::string line;

	while (std::getline(std::cin, line))
	{
		std::cout << line << std::endl;
	}
}

// convert integer to string
string convertInt(int number)
{
	// from:
	// http://www.cplusplus.com/forum/beginner/7777/

	stringstream ss;	//create a stringstream
	ss << number;		//add number to the stream
	return ss.str();	//return a string with the contents of the stream
}

// convert double to string
string convertDouble(double number)
{
	// from:
	// http://www.cplusplus.com/forum/beginner/7777/

	stringstream ss;	//create a stringstream
	ss << number;		//add number to the stream
	return ss.str();	//return a string with the contents of the stream
}

// replace a char in a string
string replaceChar(string str, char ch1, char ch2) 
{
	// from:
	// http://www.cplusplus.com/forum/beginner/33835/

	for (int i = 0; i < str.length(); ++i) 
	{
		if (str[i] == ch1)
		str[i] = ch2;
	}
	return str;
}

// make a vector of triples equal to the number of samples
// note: no need to return anything because this is passed by reference 
void intialize_vec(std::vector<triple> &vec_triple, int number_samples)
{
	// make a triple for each sample
	for (int j = 0; j < number_samples; j++)
	{
		triple mytriple;
		// set all counts to zero
		mytriple.setzero();
		vec_triple.push_back(mytriple);
	}
}

// increment the sampleindex-th sample in the map of indels or mismatch to vector of triples
// (this was originally called increment_indel bc it was just for indels but now it's for indels and mismatches, too)
void increment_indel(std::map< string,vector<triple> > &map_indel, string myindel, int sampleindex, int prevqual, bool isrev)
{
	// forward
	if (not isrev)
	{
		// increment forward count of i_th sample 
		map_indel[myindel][sampleindex].f += 1;
	}
	// reverse
	else
	{
		map_indel[myindel][sampleindex].r += 1;
	}

	map_indel[myindel][sampleindex].q += prevqual;
}

// print vcf for mismatches
void print_vcf(string prependstr, vector<counts> &vec_count, std::map< string,vector<triple> > &map_mis, string ref)
{	
	// loop over mismatches
	for (std::map< string,vector<triple> >::iterator iter = map_mis.begin(); iter != map_mis.end(); ++iter)
	{	
		// key
		string mykey = iter->first;
		// value
		vector<triple> myval = iter->second;

		string alt_format = "";

		// loop thro samples to accumulate format string 
		for(int i = 0; i < vec_count.size(); i++)	
		{		
			string rbq = ".";
			if (vec_count[i].rdf + vec_count[i].rdr > 0)
			{
				rbq = convertInt(int(vec_count[i].rbq/(vec_count[i].rdf + vec_count[i].rdr)));
			}

			string gt = "0/0";
			string abq = ".";
			string p_val = ".";

			if (myval[i].f + myval[i].r > 0) 
			{
				gt = "0/1";
				abq = convertInt(int(myval[i].q/(myval[i].f + myval[i].r)));

				if (bool_pval)
				{
					p_val = convertDouble(fastFishersExactTest(myval[i].f, myval[i].r, vec_count[i].rdf, vec_count[i].rdr));
				}
			}

			// GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
			alt_format = alt_format + "\t" + gt + ":.:" + convertInt(vec_count[i].sdp) + ":.:" + convertInt(vec_count[i].rdf + vec_count[i].rdr) + ":" + convertInt(myval[i].f + myval[i].r) + ":.:" + p_val + ":" + rbq + ":" + abq +":" + convertInt(vec_count[i].rdf) + ":" + convertInt(vec_count[i].rdr)  + ":" + convertInt(myval[i].f) + ":" + convertInt(myval[i].r);
		}

		// only print line if have a variant
		//#CHROM  POS     ID      REF	ALT     QUAL    FILTER  INFO    FORMAT
		cout << prependstr << "\t" << mykey << "\t" << "." << "\t" << "." << "\t" << "." << "\t" << "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR" << alt_format << endl; 
	}

	// if allvar flag hi AND map empty, print also non-variants
	if (allvar && map_mis.empty())
	{	
		// format field
		string alt_ref_format = "";

		// loop thro samples to accumulate format string 
		for(int i = 0; i < vec_count.size(); i++)	
		{		
			string rbq = ".";
			if (vec_count[i].rdf + vec_count[i].rdr > 0)
			{
				rbq = convertInt(int(vec_count[i].rbq/(vec_count[i].rdf + vec_count[i].rdr)));
			}

			alt_ref_format = alt_ref_format + "\t" + "0/0" + ":.:" + convertInt(vec_count[i].sdp) + ":.:" + convertInt(vec_count[i].rdf + vec_count[i].rdr) + ":" + "0" + ":.:.:" + rbq  + ":" + "0" + ":" + convertInt(vec_count[i].rdf) + ":" + convertInt(vec_count[i].rdr)  + ":" + "0" + ":" + "0";
		}

		// print
		cout << prependstr << "\t" << ref << "\t" << "." << "\t" << "." << "\t" << "." << "\t" << "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR" << alt_ref_format << endl;
	}

}

// print vcf (for indels)
void print_indel_vcf(string prependstr, vector<counts> &vec_count, std::map< string,vector<triple> > &map_indel, bool isdel, string ref)
{	
	// loop over indels
	for (std::map< string,vector<triple> >::iterator iter = map_indel.begin(); iter != map_indel.end(); ++iter)
	{	
		// key
		string mykey = iter->first;
		// value
		vector<triple> myval = iter->second;

		string alt_format = "";

		// loop thro samples to accumulate format string 
		for(int i = 0; i < vec_count.size(); i++)	
		{		
			string rbq = ".";
			if (vec_count[i].rdf + vec_count[i].rdr > 0)
			{
				rbq = convertInt(int(vec_count[i].rbq/(vec_count[i].rdf + vec_count[i].rdr)));
			}

			string gt = "0/0";
			string abq = ".";
			if (myval[i].f + myval[i].r > 0) 
			{
				gt = "0/1";
				abq = convertInt(int(myval[i].q/(myval[i].f + myval[i].r)));
			}

			// GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
			alt_format = alt_format + "\t" + gt + ":.:" + convertInt(vec_count[i].sdp) + ":.:" + convertInt(vec_count[i].rdf + vec_count[i].rdr) + ":" + convertInt(myval[i].f + myval[i].r) + ":.:.:" + rbq + ":" + abq +":" + convertInt(vec_count[i].rdf) + ":" + convertInt(vec_count[i].rdr)  + ":" + convertInt(myval[i].f) + ":" + convertInt(myval[i].r);
		}

		// if deletion
		if (isdel) 
		{
			cout << prependstr << mykey << "\t" << ref << "\t" << "." << "\t" << "." << "\t" << "." << "\t" << "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR" << alt_format << endl; 
		}
		// if insertion
		else
		{
			cout << prependstr << "\t" << ref << mykey << "\t" << "." << "\t" << "." << "\t" << "." << "\t" << "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR" << alt_format << endl; 
		}
	}
}

// read mpileup
void read_pileup(char *s, string sample_names)
{	
	string line;				// line from file
	string chr, pos, ref;			// fields from vcf:
	int number_samples;			// number of samples

	vector<string> vec_line;		// line as a vector 
	vector<counts> vec_count;		// vector of sample counts

	// maps for insertions and deletions
	std::map< string,vector<triple> > map_in;
	std::map< string,vector<triple> > map_del;

	// map for mismatches
	std::map< string,vector<triple> > map_mis;
	
	bool indelisrev;			// indel is in reverse direction
	bool isfirstiteration = 1;		// boolean for the first iteration (in which we print the header then turn this off)

	// header
	string header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

	// loop over std:in
	while ( std::getline(std::cin, line) )
	{
		if (debug)
		{
			cout << "Full line:" << endl;
			cout << line << endl;
		}

		if (isfirstiteration) cout << header; 

		// split line on tabs and get fields
		vec_line = split(line, '\t');
		chr = vec_line[0];
		pos = vec_line[1];
		ref = vec_line[2];
	
		if (vec_line.size() % 3 == 2)
		{
			/*
			This is a hack! It addresses the following problem.
			The output of samtools mpileup can look like this:
			chr1  2275948  G  8  .....,,,  >A>DCDD?  8  ,,,.,.,,   :HJDDIBD  0  *  *
			chr1  2275949  C  8  GGGGGggg  C>CDCDDA  8  gggGgGgg   :HJDDJBD  0
			chr1  2275950  T  8  .....,,,  9?9C9DDA  8  ,,,.,.,,   :BJ9DHDD  1  ,  2
			Note that sometimes pileup uses asterisks, and sometimes it uses nothing, to represent empty fields.
			However, this trips up the split function, which treats
			0 \t * \t * \n
			0 \t \t \n
			differently, causing problems.
			Splitting on tab, it says the former has 3 elts, while the latter has 2 elts
			(python, in contrast, does NOT have this problem - it's split function says 3 for both)
			To prevent this problem, I add an empty element to the line vector if vector size mod 3 == 2
			thus making it mod 3 == 0
			*/
			vec_line.push_back("");
			// DEBUG
			// cout << "DEBUG " << vec_line.size() << " " << vec_line.size() % 3 << endl;
		}

		// get number of samples
		number_samples = (vec_line.size() - 3)/3;

		// define string for first part of vcf row, which we'll need later
		string prependstr = chr + "\t" + pos + "\t" + "." + "\t" + ref;

		// if user has given sample names as input, print them
		if (isfirstiteration && sample_names.length() > 0)
		{
			// replace comma with tab
			cout << "\t" + replaceChar(sample_names,',','	');
		}

		// loop thro samples (i is the sample index)
		for (int i = 0; i < number_samples; i++)
		{
			// print the rest of the header, which depends on the number of samples (if user hasnt supplied names)
			if (isfirstiteration && sample_names.length() == 0)
			{
				cout << "\tsample_" + convertInt(i+1); 
			}

			int depth = atoi(vec_line[3+3*i].c_str());	// read depth
			string bases = vec_line[4+3*i];			// bases
			string quals = vec_line[5+3*i];			// qualities

			// make a counts struct per sample
			counts mycounts;
			// set all counts to zero
			mycounts.setzero();
			// push counts for every sample
			vec_count.push_back(mycounts);

			// skip sample if raw depth == 0
			if ( depth == 0 )
			{
				continue;
			}
			
			// set raw sample depth 
			vec_count[i].sdp += depth;

			// more variables
			int flag_indel = 0;		// flag for indel 
			bool bool_skip = 0;		// bool for skipping the line
			 				// (1 -> hit a "+" char, -1 -> hit a "-" char, 2 -> inside insertion, -2 -> inside deletion)
			string str_indel = "";		// string for indel
			int indel_counter = 0;		// number for counting bases after indel
			// totnumbase = len(bases) - 1	// tot number of bases
			int realbasecount = 0;		// this counter tracks the position of the current real base to get the corresponding qual score
			 				// indels and read segment start stops are not counted since they dont have a qual ASCII
			int prevqual = 0;		// prevqual 
							// this is hack-y - this stores the last-seen qual score, and uses it for indels

			// loop thro read char string 
			for (int k = 0; k < bases.size(); k++)
			{
				// cout << bases[k] << endl;

				// if skip, skip and turn off skip bool
				if (bool_skip)
				{
					bool_skip = 0;
					continue;
				}

				bool bool_regular_base = 0;	// bool for regular base (AaCcGgTt)

				// one after "+" or "-" char, so get length of indel
				if ( abs(flag_indel) == 1 )
				{
					// set flag for tracking is the indel on forward or rev strand?
					indelisrev = 1;

					// if it's a digit, get it and skip to next iteration of loop; OTHERWISE change flag and stay in loop
					// 48 to 57 are ASCII vals for 0 to 9
					if (int(bases[k]) >= 48 && int(bases[k]) <= 57)
					{
						// get number from ASCII, set counter equal to length of indel
						indel_counter = indel_counter*10 + int(bases[k]) - int('0');
						continue;
					}
					// if not a number, change flag and proceed
					else
					{
						// 2 value denotes we're inside an indel - set this in preparation for next iteration
						flag_indel = 2 * flag_indel;
					}
				}

				// if in indel, create str_indel
				if (indel_counter > 0 && abs(flag_indel) == 2)
				{
					// check if elements of indel are upper case - if so, call it forward
					// if (int(bases[k]) >= 65 && int(bases[k]) <= 90)
					if (isupper(bases[k]))
					{
						indelisrev = 0;
					}

					bases[k]=toupper(bases[k]);
					str_indel = str_indel + bases[k];
					indel_counter = indel_counter - 1;
					continue;
				}

				// if counter = 0 but flag = 2 this means we just finished accumulating an indel
				// so make new entry for it in the dict if it does not exist, 
				// or increment it if it already exists
				// and, importantly, proceed to parse next character as ordinary
				if (abs(flag_indel) == 2 and indel_counter == 0)
				{
					// insertion
					if (flag_indel > 0)
					{
						// initialize if not exist
						if (map_in.count(str_indel) == 0)
						{
							// vector of triples
							vector<triple> vec_triple; 

							// make an empty triple for each sample
							intialize_vec(vec_triple, number_samples);

							// map insertion to vector of triples
							map_in[str_indel] = vec_triple;
						}

						// increment indel count for i_th sample
						increment_indel(map_in, str_indel, i, prevqual, indelisrev);
					}
					// deletion
					else
					{
						// initialize if not exist
						if (map_del.count(str_indel) == 0)
						{
							// vector of triples
							vector<triple> vec_triple; 

							// make an empty triple for each sample
							intialize_vec(vec_triple, number_samples);

							// map insertion to vector of triples
							map_del[str_indel] = vec_triple;
						}

						// increment indel count for i_th sample
						increment_indel(map_del, str_indel, i, prevqual, indelisrev);
					}

					// reset counter 
					flag_indel = 0;
					// reset indel str
					str_indel = "";
				}

				if (bases[k] == '*')
				{
					// do nothing
				}
				else if (bases[k] == '$')
				{
					// read end - do nothing
				}
				else if (bases[k] == '^')
				{
					// if read start == ^ skip this char as well as the nxt char
					bool_skip = 1;
				}
				// match ref
				else if (bases[k] == '.')
				{
					// increment ref matching base count
					vec_count[i].rdf += 1;
					// increment qual score
					vec_count[i].rbq += ascii2qual(quals[realbasecount], qual_offset);
					// set regular base = true
					bool_regular_base = 1;
				}
				else if (bases[k] == ',')
				{
					vec_count[i].rdr += 1;
					vec_count[i].rbq += ascii2qual(quals[realbasecount], qual_offset);
					bool_regular_base = 1;
				}
				// indel
				else if (bases[k] == '+')
				{
					flag_indel = 1;
				}
				else if (bases[k] == '-')
				{
					flag_indel = -1;
				}
				// mismatch
				else if (toupper(bases[k]) == 'A' || toupper(bases[k]) == 'C' || toupper(bases[k]) == 'T' || toupper(bases[k]) == 'G')
				{
					// base on reverse strand
					bool is_reverse = true; 
					if (toupper(bases[k]) == bases[k]) is_reverse = false;

					// bases[k] is a char - recast to str
					string elt;
					stringstream ss;
					bases[k] = toupper(bases[k]);
					ss << bases[k];
					ss >> elt;

					// initialize if not exist
					if (map_mis.count(elt) == 0)
					{
						// vector of triples
						vector<triple> vec_triple; 

						// make an empty triple for each sample
						intialize_vec(vec_triple, number_samples);

						// map insertion to vector of triples
						map_mis[elt] = vec_triple;
					}

					// increment mismatch count for i_th sample
					increment_indel(map_mis, elt, i, ascii2qual(quals[realbasecount], qual_offset), is_reverse);

					bool_regular_base = 1;
				}

				// if regular base
				if (bool_regular_base)
				{
					// store current qual score in prevqual
					prevqual = ascii2qual(quals[realbasecount], qual_offset);
					// increment number of real bases
					realbasecount += 1;
				}
			}
			// loop thro bases

			// must take care of case where flag == 2 and we're at the end of the loop
			// if this is the case, we havent yet updated the map so do that now
			if (abs(flag_indel) == 2 and indel_counter == 0)
			{
				// insertion
				if (flag_indel > 0)
				{
					// initialize if not exist
					if (map_in.count(str_indel) == 0)
					{
						// vector of triples
						vector<triple> vec_triple; 

						// make an empty triple for each sample
						intialize_vec(vec_triple, number_samples);

						// map insertion to vector of triples
						map_in[str_indel] = vec_triple;
					}

					// increment indel count for i_th sample
					increment_indel(map_in, str_indel, i, prevqual, indelisrev);
				}
				// deletion
				else
				{
					// initialize if not exist
					if (map_del.count(str_indel) == 0)
					{
						// vector of triples
						vector<triple> vec_triple; 

						// make an empty triple for each sample
						intialize_vec(vec_triple, number_samples);

						// map insertion to vector of triples
						map_del[str_indel] = vec_triple;
					}

					// increment indel count for i_th sample
					increment_indel(map_del, str_indel, i, prevqual, indelisrev);
				}

				// reset counter 
				flag_indel = 0;
				// reset indel str
				str_indel = "";
			}

		}
		// loop thro samples

		// carriage return after header if 1st iteration
		if (isfirstiteration) 
		{
			cout << endl;
			isfirstiteration = 0;
		}

		print_vcf(prependstr, vec_count, map_mis, ref);			// print mismatches
		print_indel_vcf(prependstr, vec_count, map_in, 0, ref);		// print insertions
		print_indel_vcf(prependstr, vec_count, map_del, 1, ref);	// print deletions
		
		// clear variables
		vec_count.clear();
		map_mis.clear();
		map_in.clear();
		map_del.clear();
	}
	// loop thro line
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void pileup2vcf(int argc, char **argv) 
{    
	bool pipein = true;		// bool true if pipe input (default)
	bool bool_header = false;	// bool header
	string sample_names = "";	// string for the sample names
	char* myfile;			// vcf file (if not pipe)
	string myhelp = "Convert pileup format to vcf format, treating multi-allelic variants properly\n\nFlags:\n -h,--help: help\n -v,--verbose: verbose mode\n -q,--quality: quality offset (default: 33)\n --header: print vcf-double # lines\n -s, --samples: a comma-delimited list of sample names\n -a, --all: print non-variant lines as well as variants\n -p, --pval: print the p-value for strand bias\n\nUsage:\n To use, simply pipe your pileup file into this program. E.g.,\n cat mpileup.txt | pileup2multiallele_vcf";

	/*
	if (argc == 1)
	{
		cout << myhelp << endl;
		exit(0);
	}
	*/

	// http://stackoverflow.com/questions/4220124/c-how-to-switch-through-arguments
	// "The switch statement can't use strings. You'll need to replace it with a string of if-else if statements instead."
	// "C and C++ switch statements only operate on integral types."
	for (int i = 1; i < argc;)
	{	
		if (strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
		{
			// help
			cout << myhelp << endl;
			exit(0);
		}
		else if (strcmp(argv[i],"-a") == 0 || strcmp(argv[i],"--all") == 0)
		{
			// debug 
			allvar = true;
			i++;
		}
		else if (strcmp(argv[i],"-d") == 0 || strcmp(argv[i],"--debug") == 0)
		{
			// debug 
			debug = true;
			i++;
		}
		else if (strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0)
		{
			// verbose
			verbose = true;
			i++;
		}
		else if (strcmp(argv[i],"-q") == 0 || strcmp(argv[i],"--quality") == 0)
		{
			// quality
			qual_offset = atoi(argv[i+1]);
			i=i+2;
		}
		else if (strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--samples") == 0)
		{
			// quality
			sample_names = argv[i+1];
			i=i+2;
		}
		else if (strcmp(argv[i],"--header") == 0)
		{
			// header 
			bool_header = true;
			i++;
		}
		else if (strcmp(argv[i],"-p") == 0 || strcmp(argv[i],"--pval") == 0)
		{
			// strand bias p value 
			bool_pval = true;
			i++;
		}
		else
		{
			i++;
		}
	}

	// vcf header describe FORMAT fields
	string header="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Raw Read Depth as reported by SAMtools\">\n##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">\n##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Quality Read Depth of bases with Phred score >= 0\">\n##FORMAT=<ID=FREQ,Number=1,Type=Float,Description=\"Variant allele frequency\">\n##FORMAT=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">\n##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description=\"Average quality of reference-supporting bases (qual1)\">\n##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description=\"Average quality of variant-supporting bases (qual2)\">\n##FORMAT=<ID=RDF,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on forward strand (reads1plus)\">\n##FORMAT=<ID=RDR,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on reverse strand (reads1minus)\">\n##FORMAT=<ID=ADF,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on forward strand (reads2plus)\">\n##FORMAT=<ID=ADR,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on reverse strand (reads2minus)\">";

	if (bool_header)
	{
		cout << header << endl;
	}

	// read vcf file
	read_pileup(myfile, sample_names);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int main(int argc, char **argv) 
{    
	// to compile: g++ -O2 myprog.cpp -o out.exe
	
	pileup2vcf(argc, argv);
	return 0;
}
