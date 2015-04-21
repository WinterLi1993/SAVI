#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cstdio>
#include <limits>

using namespace std;

/*
Name			Description								Size*		Range*
char			Character or small integer.						1byte		signed: -128 to 127 unsigned: 0 to 255
short int 		(short)	Short Integer.							bytes		signed: -32768 to 32767 unsigned: 0 to 65535
int			Integer.								bytes		signed: -2147483648 to 2147483647 unsigned: 0 to 4294967295
long int 		(long)	Long integer.							bytes		signed: -2147483648 to 2147483647 unsigned: 0 to 4294967295
bool			Boolean value. It can take one of two values: true or false.		1byte		true or false
float			Floating point number.							4bytes		+/- 3.4e +/- 38 (~7 digits)
double			Double precision floating point number.					8bytes		+/- 1.7e +/- 308 (~15 digits)
long double		Long double precision floating point number.				8bytes		+/- 1.7e +/- 308 (~15 digits)
wchar_t			Wide character.								2 or 4 bytes	1 wide character
*/

/*
preprocessor directives
#define #error #import	#undef #elif #if #include #using #else #ifdef #line #endif #ifndef #pragma
*/

// global variables
bool rdplusad = false;		// use RD+AD as total depth, not SDP
bool hybrid = false;		// use RD+AD as total depth for sample 1, SDP as total depth for sample 2

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

// convert integer to string
string convertInt(int number)
{
	// from:
	// http://www.cplusplus.com/forum/beginner/7777/

	stringstream ss;	//create a stringstream
	ss << number;		//add number to the stream
	return ss.str();	//return a string with the contents of the stream
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//**************************************//
// file reading functions		//
//**************************************//

// read vcf file
/*
http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
There are 8 fixed, mandatory columns:
#CHROM
POS
ID
REF
ALT
QUAL
FILTER
INFO

followed by:
FORMAT
SAMPLE1
SAMPLE2
...
*/

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

// read vcf from argument s OR std:in if pipein
void read_vcf(char *s, bool leadingone, int samplenum_1, string &samplestr, bool pipein)
{	
	// fields from vcf:
	string chrom, pos, id, ref, alt, qual, filter, info, format;

	string samp_1;				// sample 1 if the 2samples flag
	string samp_2;				// sample 2 if the 2samples flag

	string line;				// line from file
	ifstream myfile (s);			// ifstream from file
	vector<string> vec_line;		// line as a vector 
	vector<string> vec_format_field;	// vector of substrings from FORMAT field
	vector<string> vec_sample_indices;	// vector of sample indices from the 2samples flag
	vector<string> vec_s1_field;		// vector of substrings from first sample field
	vector<string> vec_s2_field;		// vector of substrings from second sample field

	int samplenum_2 = 1;			// index of second sample for the 2samples flag

	char delim = ':';			// the delimiter for the FORMAT field
	char delim_tab = '	';		// the delimiter for the line
	char delim_com = ',';			// the delimiter for the 2samples flag
	bool first_read_line = true;		// a boolean which is true for the first line
	bool wellformed_line;			// a boolean which is true if the line is well-formed

	int index_ref_qual = 0;			// index of reference qual
	int index_var_qual = 0;			// index of variant qual
	int index_ref = 0; 			// index of ref depth
	int index_var = 0; 			// index of variant depth
	int index_tot = 0; 			// index of tot depth

	// if 2samples flag not empty
	if ( samplestr.compare("") != 0 )
	{
		vec_sample_indices = split(samplestr, delim_com);
		samplenum_1 = atoi(vec_sample_indices[0].c_str());
		samplenum_2 = atoi(vec_sample_indices[1].c_str());
	}

	// loop if file open OR pipe from std:in
	if (pipein or myfile.is_open())
	{
		while ( (pipein and std::getline(std::cin, line)) or (!pipein and !myfile.eof() ) )
		{
			if (not pipein) getline(myfile,line);

			// if line begins with # character, skip this iteration
			if ( line.substr(0,1).compare("#") == 0 )
			{
				continue;
			}

			// split line on tabs and get fields
			vec_line = split(line, delim_tab);
			chrom = vec_line[0];
			pos = vec_line[1];
			id = vec_line[2];
			ref = vec_line[3];
			alt = vec_line[4];
			qual = vec_line[5];
			filter = vec_line[6];
			info = vec_line[7];
			format = vec_line[8];
			
			// if 2samples flag not empty
			if ( samplestr.compare("") != 0 and (samplenum_1 > 0 && 8+samplenum_1 <= vec_line.size()) and (samplenum_2 > 0 && 8+samplenum_2 <= vec_line.size()) )
			{
				samp_1 = vec_line[8+samplenum_1];	
				samp_2 = vec_line[8+samplenum_2];	
			}
			// else if only one sample: check if sample in range and, if so, get sample field
			else if (samplenum_1 > 0 && 8+samplenum_1 <= vec_line.size())
			{
				samp_1 = vec_line[8+samplenum_1];	
			}
			else
			{
				cerr << "sample index out of range" << endl;
				exit(1);
			}

			// get format field
			vec_format_field = split(format,delim);

			// get sample 1
			vec_s1_field = split(samp_1,delim);

			// get sample 2 if requested
			if ( samplestr.compare("") != 0 )
			{
				vec_s2_field = split(samp_2,delim);
			}

			// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
			// ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
			// ##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Raw Read Depth as reported by SAMtools">
			// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 15">
			// ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
			// ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
			// ##FORMAT=<ID=FREQ,Number=1,Type=Float,Description="Variant allele frequency">
			// ##FORMAT=<ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
			// ##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">
			// ##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
			// ##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
			// ##FORMAT=<ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
			// ##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
			// ##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">

			// if first line, get indicies of desired fields
			// just need to do this once
			if (first_read_line)
			{
				for (int i = 0; i < vec_format_field.size(); i++)	
				{	
					if (vec_format_field[i].compare("RBQ") == 0)
					{
						index_ref_qual = i;
					}
					else if (vec_format_field[i].compare("ABQ") == 0)
					{
						index_var_qual = i;
					}
					else if (vec_format_field[i].compare("SDP") == 0)
					{
						index_tot = i;
					}
					else if (vec_format_field[i].compare("RD") == 0)
					{
						index_ref = i;
					}
					else if (vec_format_field[i].compare("AD") == 0)
					{
						index_var = i;
					}
				}

				first_read_line = false;
			}

			wellformed_line = true;
			if ( index_ref_qual <= vec_s1_field.size() && index_var_qual <= vec_s1_field.size() && index_ref <= vec_s1_field.size() && index_var <= vec_s1_field.size() )
			{
				// if field has . character, not well formed
				if ( vec_s1_field[index_ref].substr(0,1).compare(".") == 0 or vec_s1_field[index_var].substr(0,1).compare(".") == 0 )
				{
					wellformed_line = false;
				}
			}
			else
			{
				wellformed_line = false;
			}

			// if fields exist, print q v t
			if ( wellformed_line )
			{
				int myqual = ( atoi(vec_s1_field[index_ref_qual].c_str()) + atoi(vec_s1_field[index_var_qual].c_str()) )/2;
				// total depth at the position
				int mytot;

				if ( rdplusad || hybrid ) 
				{
					// if rdplusad, use mytot = RD + AD
					mytot = atoi(vec_s1_field[index_ref].c_str()) + atoi(vec_s1_field[index_var].c_str());
				}
				else
				{
					// else use mytot = SDP
					mytot = stoi(vec_s1_field[index_tot]);
				}

				// 1 q v d 
				if ( leadingone )
				{
					cout << "1\t" << myqual << "\t" << vec_s1_field[index_var] << "\t" << mytot;
				}
				// q v d 
				else
				{
					cout << myqual << "\t" << vec_s1_field[index_var] << "\t" << mytot;
				}
			}
			// else print 0's
			else
			{
				// 1 q v d 
				if ( leadingone )
				{
					cout << "1\t0\t0\t0";
				}
				// q v d 
				else
				{
					cout << "0\t0\t0";
				}
			}

			// if second sample specified, append to line
			if ( samplestr.compare("") != 0 )
			{
				wellformed_line = true;
				if ( index_ref_qual <= vec_s2_field.size() && index_var_qual <= vec_s2_field.size() && index_ref <= vec_s2_field.size() && index_var <= vec_s2_field.size() )
				{
					// if field has . character, not well formed
					if ( vec_s2_field[index_ref].substr(0,1).compare(".") == 0 or vec_s2_field[index_var].substr(0,1).compare(".") == 0 )
					{
						wellformed_line = false;
					}
				}
				else
				{
					wellformed_line = false;
				}

				if ( wellformed_line )
				{
					int mytot;

					if ( rdplusad ) 
					{
						// if rdplusad, use mytot = RD + AD
						mytot = atoi(vec_s2_field[index_ref].c_str()) + atoi(vec_s2_field[index_var].c_str());
					}
					else
					{
						// else use mytot = SDP
						mytot = stoi(vec_s2_field[index_tot]);
					}

					int myqual = ( atoi(vec_s2_field[index_ref_qual].c_str()) + atoi(vec_s2_field[index_var_qual].c_str()) )/2;
					cout << "\t" << myqual << "\t" << vec_s2_field[index_var] << "\t" << mytot << endl;
				}
				// else print 0's
				else
				{
					cout << "\t0\t0\t0" << endl;
				}
			}
			else
			{
				cout << endl;
			}

			/*
			for (int k = 0; k < format.size(); k++)
			{
				cout << format[k] << endl;
			}
			*/
		}
		if (not pipein) myfile.close();
	}
	else cerr << "Unable to open file";    		
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void make1p(int argc, char **argv) 
{    
	bool leadingone = false;	// bool print leading one
	bool verbose = false;		// bool verbose
	bool pipein = true;		// bool true if pipe input (default)
	char* myfile;			// vcf file (if not pipe)
	int samplenum = 1;		// index of the sample for which to produce q v t 
	string samplestr = "";		// comma-delimited list of sample indices (if more than one)

	string myhelp = "Make the prior files for savi\nTake a vcf input and produce a file with columns: quality, variant depth, total depth\nNOTE: the sub-fields of the FORMAT field of the vcf file MUST BE THE SAME throughout the file\n\nFlags:\n -1: print a column of leading ones\n -f,--vcf: vcf file (if you leave this flag off, the program expects input piped from std:in)\n -s,--sample: sample index (e.g., 1)\n -2s,--2sample: indices of two samples in a comma-delimited list (e.g., 2,3)\n -h,--help: help\n -rdplusad,--rdplusad: use reference agreeing reads plus variant calling reads (RD+AD) for the total depth, rather than SDP\n -hybrid,--hybrid: if 2sample flag, use reference agreeing reads plus variant calling reads (RD+AD) for the total depth in FIRST sample and use SDP as total depth in SECOND sample; if not, this flag is identical to --rdplusad\n -v,--verbose: verbose mode\n\nUsage: To use, simply pipe your vcf file into this program. E.g.,\n cat variants.vcf | make_qvt -1";

	if (argc == 1)
	{
		cout << myhelp << endl;
		exit(0);
	}

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
		else if (strcmp(argv[i],"-1") == 0)
		{
			// print leading one
			i++;
			leadingone = true;
		}
		else if (strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--vcf") == 0)
		{
			// vcf file
			i++;
			// use file instead of std:in
			pipein = false;
			myfile = argv[i];
		}
		else if (strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--sample") == 0)
		{
			// sample index
			i++;
			samplenum = atoi(argv[i]);
		}
		else if (strcmp(argv[i],"-2s") == 0 || strcmp(argv[i],"--2sample") == 0)
		{
			// sample indices
			i++;
			samplestr = argv[i];
		}
		else if (strcmp(argv[i],"-rdplusad") == 0 || strcmp(argv[i],"--rdplusad") == 0)
		{
			// use RD + AD, not SDP, as total depth 
			i++;
			rdplusad = true;
		}
		else if (strcmp(argv[i],"-hybrid") == 0 || strcmp(argv[i],"--hybrid") == 0)
		{
			// use RD + AD as tot depth for first sample , and SDP as total depth for second
			i++;
			hybrid = true;
		}
		else if (strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0)
		{
			// verbose
			i++;
			verbose = true;
		}
		else
		{
			i++;
		}
	}

	/*
	cout << leadingone << endl;
	cout << myfile << endl;
	cout << samplenum << endl;
	*/

	// read vcf file
	read_vcf(myfile, leadingone, samplenum, samplestr, pipein);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int main(int argc, char **argv) 
{    
	// to compile: g++ -O2 myprog.cpp -o out.exe
	
	make1p(argc, argv);
	// read_from_std();
	return 0;
}
