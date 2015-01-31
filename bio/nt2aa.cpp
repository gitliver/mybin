#include <iostream>
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
#include <iomanip>

using namespace std;

/*
About: 
Translate nucleotides to amino acids using codon table
This script takes one argument - the codon amino acid 
translations in the file: humangencode.txt
As input, it takes an uppercase nucleotide sequence piped in.

Example:
If we enter:
$ echo -e "acggtcctt\nAGT"
acggtcctt
AGT

Useage:
$ echo -e "acggtcctt\nAGT" | tr [:lower:] [:upper:] | nt2aa humangencode.txt
TVL
S

Notes:
Anything that's not a proper codon will be rendered as an "X". For example,
$ echo AGANN | nt2aa humangencode.txt
RX

Notes:
Dont forget to compile: g++ -O2 nt2aa.cpp -o nt2aa
*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//**************************************//
// functions				//
//**************************************//

// read humangencode.txt (codon	\t AA)
void read_gene_code(const char *f, map<string, string> &x) 
{
	ifstream is(f);
	string s1, s2;
	
	x.clear();
	while (1) 
	{
		is >> s1;
		if (is.eof()) break;
		is >> s2;

		x[s1] = s2; 
	}
}

void nt2aa(string &s, map<string, string> &g, string &r) 
{
	// nucleotide to amino acid

	string e;
	
	r = "";
	for(unsigned i = 0; i < s.length(); i += 3) {
		e = s.substr(i, 3);
		if (g.find(e) != g.end()) r = r + g[e];
		else r = r + 'X';
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void main_nt2aa(char **argv) 
{
	/*
	arguments:
	
	$1 humangencode.txt (file: codon \t AA)
	
	to compile: g++ -O2 nt2aa.cpp -o nt2aa
	
	to run, e.g.:
	cat nt_file | nt2aa humangencode.txt	
	*/
			
	string ntseq;	// a line piped in
	string aaseq;	// a line piped out
	
	// istringstream iss;	// input stream

	map<string, string> g;	// codon to amino acid mapping
	
	// cerr << "reading gene code" << endl;
	read_gene_code(argv[1], g);	// read humangencode.txt	
	
	while (true)
	{
		// read line
		std::getline(std::cin, ntseq);

		// read file piped in until hit an empty line				
		if (ntseq.empty())
		{
			break;
		}
		else 
		{
			nt2aa(ntseq, g, aaseq);
			cout << aaseq << endl;
		}
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    
int main(int argc, char **argv) 
{    
	// main_revc(argv);
	main_nt2aa(argv);
	return 0;
}
