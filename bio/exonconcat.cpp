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
Concatenate Exons

Example:
Suppose we have a list of exon positions for some given genes on a particular chromosome, 
and we want to stitch these together to get the sequences. 
For example, say we have ccds exon positions for the genes PAX7 and RIMS3:

chr1	CCDS186.1	PAX7	+	18830685	18830769	
chr1	CCDS186.1	PAX7	+	18833384	18833619	
chr1	CCDS186.1	PAX7	+	18834192	18834321	
chr1	CCDS186.1	PAX7	+	18835318	18835452	
chr1	CCDS186.1	PAX7	+	18890835	18891034	
chr1	CCDS186.1	PAX7	+	18899734	18899899	
chr1	CCDS186.1	PAX7	+	18902175	18902377	
chr1	CCDS186.1	PAX7	+	18934713	18935120
	
chr1	CCDS30687.1	RIMS3	-	40864776	40864988	
chr1	CCDS30687.1	RIMS3	-	40867069	40867208	
chr1	CCDS30687.1	RIMS3	-	40867530	40867631	
chr1	CCDS30687.1	RIMS3	-	40871328	40871440	
chr1	CCDS30687.1	RIMS3	-	40874175	40874316	
chr1	CCDS30687.1	RIMS3	-	40879968	40880184	         

We have to give exonconcat two inputs: one is simply the whole chromosome 1 sequence; 
the other is a file in the format start_position, end_position, id, sense:

18830685	18830769	1	1
18833384	18833619	1	1
18834192	18834321	1	1
18835318	18835452	1	1
18890835	18891034	1	1
18899734	18899899	1	1
18902175	18902377	1	1
18934713	18935120	1	1
40864776	40864988	2	0
40867069	40867208	2	0
40867530	40867631	2	0
40871328	40871440	2	0
40874175	40874316	2	0
40879968	40880184	2	0

where IMPORTANTLY the file must be sorted by (unique) id then start_position, as it is here. 
For the sense, 1 --> + and 0 --> -. A reverse complement is performed for sequences not in the + sense. 
Let's call the above file "test", and run it through exonconcat:

Useage:
$ exonconcat /path/chr/1/seq test | tr '[:lower:]' '[:upper:]'

For each id, a row with the stitched together sequence has been outputted. 
Is it right? We can check the CCDS Database (http://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi) 
to make sure.

Notes:
Dont forget to compile: g++ -O2 exonconcat.cpp -o exonconcat		
*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

//**********************************************//
// Vlad's Templates and Functions		//
//**********************************************//

// a vector
template <typename T>
struct vec 
{
	typedef vector<T> type;
};

// a matrix
template <typename T>
struct mat 
{
	typedef vector<vector<T> > type;
};	

template<typename T>
unsigned length(T &v) 
{
	return v.size();
}
	
template<typename T>
unsigned nrows(T &m) 
{
	return m.size();
};

template<typename T>
unsigned ncols(T &m) 
{
	return m[0].size();
};

template<typename T1, typename T2>
void const_vec(T1 &v, unsigned l, T2 z) 
{
	v.clear();
	v.resize(l, z);
};

template<typename T1, typename T2>
void const_mat(T1 &m, unsigned r, unsigned c, T2 z) 
{
	m.clear();
	m.resize(r, typename vec<T2>::type(c, z));
};

// a function for reading files
template<typename T>
void read_text_mat(istream &is, T &v) 
{
	typedef typename T::value_type::value_type value_type;
	
	typename vec<value_type>::type t;
	istringstream iss;
	string s;
	value_type d;
	bool f;
	unsigned i;
	
	f = true;
	i = 0;
	do {
		getline(is, s);
		if (is.eof()) break;
		iss.clear();
		iss.str(s);

		if (f) 
		{
			while (!iss.eof()) 
			{
				iss >> d >> ws;
				t.push_back(d);
			}
			f = false;
		}
		else
			for(unsigned i = 0; i < t.size(); i++)
				iss >> t[i];
		
		v.push_back(t);
	} while (1);
};

template<typename T>
void read_text_mat(char *s, T &m) 
{
	ifstream is(s);
	read_text_mat(is, m);
};

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

bool comp(const vec<unsigned>::type &v1, const vec<unsigned>::type &v2) 
{
	return v1[3] < v2[3];
}

char ntcompl(char c) 
{
	static char x[] = "ACGTRYPQHSJKLMN";
	static char y[] = "TGCAYRPQSHMLKJN";
	char *p;
	
	p = strchr(x, toupper(c));
	if (p) return y[p - x]; else return 'N';
}

void nt2aa(string &s, map<string, string> &g, string &r) 
{
	string e;
	
	r = "";
	for(unsigned i = 0; i < s.length(); i += 3) {
		e = s.substr(i, 3);
		if (g.find(e) != g.end()) r = r + g[e];
		else r = r + 'X';
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void reverse_complement(string &str) 
{
	// output the rev comp of the string
	
	string letter;
  	string::reverse_iterator rit;
  	for (rit=str.rbegin(); rit < str.rend(); rit++)
  	{
  		letter=*rit;
	
  		if (letter.compare("A") == 0)
	  		cout << "T";
	  	else if (letter.compare("T") == 0)
	  		cout << "A";
	  	else if (letter.compare("C") == 0)
	  		cout << "G";
	  	else if (letter.compare("G") == 0)
	  		cout << "C";
	  	else if (letter.compare("a") == 0)
	  		cout << "t";
	  	else if (letter.compare("t") == 0)
	  		cout << "a";
	  	else if (letter.compare("c") == 0)
	  		cout << "g";
	  	else if (letter.compare("g") == 0)
	  		cout << "c";		  		  				
  	}
};

string reverse_complement_seq(const string &str) 
{
	// return the rev comp of the string, but preserve dashes
	  	
  	string s;
  	
  	s.resize(str.size());
  	
  	int r;
  	
  	for (int i = str.size() - 1; i >= 0; --i)
  	{
  		r = str.size() - i - 1;
  		
  		switch(str[i]) 
  		{
  			case 'A':
  				s[r] = 'T';
  				break;
  			case 'T':
  				s[r] = 'A';
  				break;
  			case 'C':
  				s[r] = 'G';
  				break;
  			case 'G':
  				s[r] = 'C';
  				break;
  			case '-':
  				s[r] = '-';
  				break;  				  				
  			case 'a':
  				s[r] = 't';
  				break;
  			case 't':
  				s[r] = 'a';
  				break;
  			case 'c':
  				s[r] = 'g';
  				break;
  			case 'g':
  				s[r] = 'c';
  				break;  			
  		}
  	}
  	
  	return s;

	/*
	string letter;
	string::reverse_iterator rit;
	stringstream s;
	for (rit=str.rbegin(); rit < str.rend(); rit++)
	{
  		letter=*rit;
	
  		if (letter == "A")
	  		s << "T";
	  	else if (letter == "T")
	  		s << "A";
	  	else if (letter == "C")
	  		s << "G";
	  	else if (letter == "G")
	  		s << "C";
	  	else if (letter == "a")
	  		s << "t";
	  	else if (letter == "t")
	  		s << "a";
	  	else if (letter == "c")
	  		s << "g";
	  	else if (letter == "g")
	  		s << "c";		  		  				
  	}
  	return s.str();
  	*/
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void concat_exons(char **argv) 
{
	
	/*
	This program stitches together the exons of a ccds
	
	It takes two arguments:
	$1: a file containing the full nucleotide sequence for a chr on a single line
	$2: a tab delimited file of exons of the form: begin (1-based counting), end, ccdsid, sense --- sorted by ccdsid then begin position
	
	It outputs a file of the form: ccdsid, ccds_seq
	*/

	unsigned i;		// loop index
	unsigned ccdsid;	// ccdsid
	int sense;		// the sense
	string line;		// the full sequence for a chr read in from $1
	string seq;		// the seq of the given ccds (subset of line)
	mat<unsigned>::type m1; // 4x4 matrix read in from $2 - ASSUME SORTED BY CCDS THEN POSITION
  	
  	//cerr << "reading 1" << endl; // read the sequence  	
  	ifstream myfile (argv[1]);
  	if (myfile.is_open())
  	{
		getline (myfile,line);
		myfile.close();
  	}

	//cerr << "reading 2" << endl; // read the 16_exon file
	read_text_mat(argv[2], m1);
	
	seq.clear(); // clear the ccds seq string
    
	ccdsid = m1[0][2]; // the first ccds - this is a variable that changes only when the ccds changes

	// cout << ccdsid << "\t";	// tab
	
	for(i = 0; i < nrows(m1); i++)
	{
		if (m1[i][2]==ccdsid) // if ccds stays the same, grow seq
		{
			seq += line.substr(m1[i][0]-1, m1[i][1]-m1[i][0]+1);	
		}
		else // if ccds changes, print previous seq
		{
			if (m1[i-1][3]==1) // prev sense = 1 
				cout << seq << endl;     			
			else if (m1[i-1][3]==0) // prev sense = 0
			{
				reverse_complement(seq);
				cout << endl; 
			}

			seq.clear(); 
			seq += line.substr(m1[i][0]-1, m1[i][1]-m1[i][0]+1);

			ccdsid = m1[i][2];
			//cout << ccdsid << "\t"; 
		}
	}
		
	// print final one
	if (m1[i-1][3]==1) // sense = 1 
		cout << seq << endl;     			
	else if (m1[i-1][3]==0) // sense = 0
	{
		reverse_complement(seq);
		cout << endl; 
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    
int main(int argc, char **argv) 
{    
	concat_exons(argv);		
	return 0;
}
