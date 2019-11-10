#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _WIN32
#include "include/getopt.h"
#else
#include <unistd.h>
#endif

using namespace std;

vector<vector<string>> table(1);

int main(int argc, char *argv[]) {
	cin.sync_with_stdio(0);
	cin.tie(0);

	if (argc < 3) {
		cout << "Usage: txtFilter -i <input file> -o <output file> -m <mode> -r <raw file> -p <pep threshold>\n";
		exit(1);
	}

	string input_filename, output_filename, raw_filename;
	int mode = 0;
	double pep_threshold = -1.0;

	int c;
	while ((c = getopt(argc, argv, "i:o:m:r::p::")) != -1) {
		switch (c) {
		case 'i':
			input_filename = string(optarg);
			break;
		case 'o':
			output_filename = string(optarg);
			break;
		case 'm':
			mode = stoi(optarg);
			break;
		case 'r':
			raw_filename = string(optarg);
			break;
		case 'p':
			pep_threshold = stod(optarg);
			break;
		default:
			cout << "Usage: txtFilter -i <input file> -o <output file> -m <mode> -r <raw file> -p <pep threshold>\n";
			exit(1);
		}
	}

	if (mode < 0 || mode > 2) {
		cout << "Error: <mode> must be 0, 1, or 2\n";
		exit(1);
	}

	if (mode == 1 && raw_filename == "") {
		cout << "Error: <raw file> cannot be empty if <mode> is 1\n";
		exit(1);
	}

	if (mode == 2 && pep_threshold < 0) {
		cout << "Error: <pep threshold> must be positive\n";
		exit(1);
	}

	fstream infile(input_filename.c_str(), ios::in);
	if (!infile.is_open()) {
		cout << "Error: unable to open input file\n";
		exit(2);
	}

	fstream outfile(output_filename.c_str(), ios::out | ios::trunc);
	if (!outfile.is_open()) {
		cout << "Error: unable to create output file\n";
		exit(2);
	}

	cout << "Processing input...";

	int is_csv = (input_filename.substr(input_filename.size() - 3) == "csv") ? 1 : 0;

	unordered_map<string, int> headers;
	int idx = 0;

	string input;
	getline(infile, input);

	stringstream ss(input);
	while (ss.good()) {
		if (is_csv) getline(ss, input, ',');
		else getline(ss, input, '\t');
		table[0].push_back(input);
		headers.insert({ input, idx++ });
	}

	while (1) {
		getline(infile, input);
		if (input == "") break;
		vector<string> line;

		stringstream ss(input);
		while (ss.good()) {
			if (is_csv) getline(ss, input, ',');
			else getline(ss, input, '\t');
			line.push_back(input);
		}

		if (mode == 0) {
			table.push_back(line);
		}
		else if (mode == 1) {
			string rf = line[headers.find("Raw file")->second];
			if (rf == raw_filename)
				table.push_back(line);
		}
		else if (mode == 2) {
			string pep_str = line[headers.find("PEP")->second];
			if (pep_str != "NaN" && pep_str != "" && stod(pep_str) <= pep_threshold)
				table.push_back(line);
		}
	}

	for (int i = 0; i < table[0].size(); i++) {
		outfile << table[0][i];
		if (i < table[0].size() - 1)
			outfile << ',';
	}
	outfile << '\n';

	for (int i = 1; i < table.size(); i++) {
		for (int j = 0; j < table[0].size(); j++) {
			outfile << table[i][j];
			if (j < table[0].size() - 1)
				outfile << ',';
		}
		outfile << '\n';
	}

	cout << "Done\n";
	infile.close();
	outfile.close();

	return 0;
}
