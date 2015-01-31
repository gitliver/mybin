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
Reverse Complement a Sequence

Example:
Suppose
$ cat example/seq.txt 
TGCGGG
AAA
TTATG
C
G

Useage:
$ cat example/seq.txt | revc 
CCCGCAT
TTT
CATAA
G
C

Notes:
Dont forget to compile: g++ -O2 revc.cpp -o revc 
*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// functions				    //
//******************************************//

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
};


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void main_revc(char **argv) 
{	
	// a line piped in
	string line;
	// input stream
	// istringstream iss;
	
	while (true)
	{
		// read line
		std::getline(std::cin, line);

		// read file piped in until hit an empty line				
		if (line.empty())
		{
			break;
		}
		else 
		{
			cout << reverse_complement_seq(line) << endl;
		}
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    
int main(int argc, char **argv) 
{    
	main_revc(argv);

  	return 0;
}
