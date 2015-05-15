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
void read_vcf(char *s, char *t, bool header)
{	
	// fields from vcf:
	string chrom, pos, id, ref, alt, qual, filter, info, format, restofline;

	string line;				// line from std:in
	string fileline;			// line from file
	ifstream myfile (s);			// ifstream from file
	ofstream myoutfile;			// ofstream for output file

	vector<string> vec_line;		// line as a vector 

	char delim_tab = '	';		// the delimiter for the line

	if (header)
	{
		myoutfile.open(t);
	}

	// loop if file open OR pipe from std:in
	while ( std::getline(std::cin, line) )
	{
		// if line begins with # character, skip this iteration
		if ( line.substr(0,1).compare("#") == 0 )
		{
			if (header)
			{
				myoutfile << line << endl;
			}
			continue;
		}

		// split line from std:in on tabs and get fields
		vec_line = split(line, delim_tab);

		// get line from file
       		if (myfile.is_open())
        	{
			getline(myfile,fileline);
		}

		// chrom = vec_line[0];
		// pos = vec_line[1];
		// id = vec_line[2];
		// ref = vec_line[3];
		// alt = vec_line[4];
		// qual = vec_line[5];
		// filter = vec_line[6];
		// info = vec_line[7];
		// format = vec_line[8];

		for (int i = 0; i <= 6; i++)	
		{	
			cout << vec_line[i] << "\t";
		}

		// print INFO field with stuff from the file stuck on
		cout << vec_line[7] << ";" << fileline << "\t";

		for (int i = 8; i < vec_line.size(); i++)	
		{	
			if (i < vec_line.size() - 1)
			{
				cout << vec_line[i] << "\t";
			}
			else
			{
				cout << vec_line[i] << endl;
			}
		}
	}

	myfile.close();

	if (header)
	{
		myoutfile.close();
	}
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void addtoinfo(int argc, char **argv) 
{    
	bool verbose = false;		// bool verbose
	bool header = false;		// bool true if save header to file 
	char* myfile;			// textfile
	char* myoutfile;		// file to write the output header to 

	string myhelp = "Paste the contents of a text file onto the INFO field of a vcf piped in from std:in.\nNOTE: the text file must have the same number of lines as the non-header portion of the vcf.\nIt will be attached to the preexisting INFO field with a semi-colon delimiter.\n -f,--file: text file\n -o,--header: save the vcf header as the argument of this flag\n -h,--help: help\n -v,--verbose: verbose mode";

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
		else if (strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--header") == 0)
		{
			// file
			i++;
			myoutfile = argv[i];
			header = true;
		}
		else if (strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--file") == 0)
		{
			// file
			i++;
			myfile = argv[i];
		}
		else if (strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0)
		{
			// verbose
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
	read_vcf(myfile, myoutfile, header);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int main(int argc, char **argv) 
{    
	// to compile: g++ -O2 myprog.cpp -o out.exe
	
	addtoinfo(argc, argv);
	// read_from_std();
	return 0;
}
