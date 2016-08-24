#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>

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
// minimum depth sought
int mindepth = 10;
// minimum depth of variant reads sought
int minaltdepth = 2;
// require the first sample (assumed to be normal) to have min read depth 
bool normaldepth = false;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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

// parse read string to count variants, etc
int parseRead(const std::string &s)
{

	// skip a character
	bool skip = false;
	// variant count
	int count = 0;

	// parse read
	for (int i=0; i< s.length(); i++) 
	{
		if (skip)
		{
			skip = false;
		}
		else
		{
			// read start => ignore next character
			if (s[i] == '^')
			{
				skip = true;
			}
			// indel (this could be worked out - 2 or 3 digit indels are a pain to treat)
			// for now, indels will be counted as mismatches
			else if (s[i] == '+' or s[i] == '-')
			{
				;
			}
			// mismatch
			else if (s[i] == 'A' or s[i] == 'a' or s[i] == 'C' or s[i] == 'c' or s[i] == 'T' or s[i] == 't' or s[i] == 'G' or s[i] == 'g' or s[i] == 'N' or s[i] == 'n')
			{
				count++;
			}
		}
	}
	return count;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

// filter pileup file
void filterpileup() 
{    
	// line
	std::string line;

	while (std::getline(std::cin, line))
	{
		// at least one sample has min depth bool
		bool is_gt_mindepth = false;
		// at least one sample has min alt depth bool
		bool is_gt_minaltdepth = false;
		// vector of fields comprising the line
		vector<string> x;
		// character to split line on: tab
		char delim = '	';

		x = split(line, delim);
		for (int i = 0; i < x.size(); i++)	
		{	
			// ignore first three columns: chr pos depth
			// depth field
			if (i > 2 && i % 3 == 0)
			{
				// if any sample meets min depth requirements, turn this to true
				if (atoi(x[i].c_str()) >= mindepth)
				{
					is_gt_mindepth = true;
				}
			}
			// read field
			else if (i > 2 && i % 3 == 1)
			{
				// if any sample meets min alt depth requirements, turn this to true
				if (parseRead(x[i]) >= minaltdepth)
				{
					is_gt_minaltdepth = true;
				}
			}
			/*
			cout << x[i] << endl;
			cout << "is min depth(";
			cout << mindepth;
			cout << "): ";
			cout << is_gt_mindepth << endl;
			cout << "is min alt depth(";
			cout << minaltdepth;
			cout << "): ";
			cout << is_gt_minaltdepth << endl;
			*/
		}

		// print line if satisfies requirements of mindepth and min alt depth
		if (is_gt_mindepth && is_gt_minaltdepth)
		{
			std::cout << line << std::endl;
		}
	}
}    

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

// parse arguments
void argparse(int argc, char **argv) 
{    
	string myhelp = "Filter mpileup file based on coverage depths and number of variants\n\nFlags:\n -h,--help: help\n -m,--mindepth [INT]: if at least one sample satisfies min read depth, print line (default: 10)\n -a,--minaltdepth [INT]: if at least one sample satisfies min alt read depth (default: 2)\n -n,--normaldepth: require the first sample (assumed to be normal) to have min read depth (default: off)\n\nUsage:\n To use, simply pipe your pileup file into this program. E.g.,\n cat mpileup.txt | filter_pileup";

	/*
	if (argc == 1)
	{
		cout << myhelp << endl;
		exit(0);
	}
	*/

	for (int i = 1; i < argc;)
	{	
		if (strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
		{
			// help
			cout << myhelp << endl;
			exit(0);
		}
		else if (strcmp(argv[i],"-m") == 0 || strcmp(argv[i],"--mindepth") == 0)
		{
			// min alt read depth to print line
			mindepth = atoi(argv[i+1]);
			i=i+2;
		}
		else if (strcmp(argv[i],"-a") == 0 || strcmp(argv[i],"--minaltdepth") == 0)
		{
			// min read depth to print line
			minaltdepth = atoi(argv[i+1]);
			i=i+2;
		}
		else if (strcmp(argv[i],"-n") == 0 || strcmp(argv[i],"--normaldepth") == 0)
		{
			// require the first sample (assumed to be normal) to have min read depth 
			normaldepth = true;
			i=i+1;
		}
		else
		{
			i++;
		}
	}

	// filter pileup file
	filterpileup();
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int main(int argc, char **argv) 
{    
	// to compile: g++ -O2 myprog.cpp -o out
	
	argparse(argc, argv);
	return 0;
}
