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
Get substrings from a sequence 

Example:
Suppose
$ cat example/test_seq.txt
AAAAGCCGAGAGACCTTTTTTTTTT

$ cat example/positions.txt 
1	5
6	6
10	11

Useage:
$ cat example/positions.txt | ./get_nt example/test_seq.txt
AAAAG
C
GA

Notes:
Dont forget to compile: g++ -O2 get_nt.cpp -o get_nt
*/

// Get Substring
void get_nucleotide(string &full_seq, unsigned vstart, unsigned vfinish)
{
	string actual_ref = full_seq.substr(vstart - 1, vfinish - vstart +1);	
	cout << actual_ref << endl;
}

//
void get_nt(char **argv) 
{	
	// a line piped in
	string line;
	istringstream iss;

	// position
	unsigned start, finish;

	// sequence
	string full_seq;
	
	ifstream myfile (argv[1]);

	if (myfile.is_open())
	{
		getline (myfile,full_seq);
		myfile.close();
	}
	else cout << "Unable to open file";  	 	

	while (true)
	{
		// read line, get current vid		
		std::getline(std::cin, line);

		if (line.empty()) 
			break;
		else 
		{
			iss.clear();
			iss.str(line);
			iss >> start >> ws >> finish;
			get_nucleotide(full_seq, start, finish);
		}
	}

}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    
int main(int argc, char **argv) 
{    
	get_nt(argv);
	return 0;
}
