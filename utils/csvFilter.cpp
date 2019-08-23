#include <algorithm>
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

#define Q_THRESHOLD 0.01
typedef pair<double, vector<string>*> pair_dvs;

using namespace std;

vector<vector<string>> table(1);

int main(int argc, char *argv[]) {
	// Speed up input
	cin.sync_with_stdio(0);
	cin.tie(0);

	if (argc < 4) {
		cout << "Usage: csvFilter -i <input csv> -o <output csv> -m <mode> [-a <min rt> -b <max rt>]\n";
		exit(1);
	}

	string input_filename, output_filename;
	int mode;
	double min_rt, max_rt;

	int c;
	while ((c = getopt(argc, argv, "i:o:m:a:b:")) != -1) {
		switch (c) {
		case 'i':
			input_filename = string(optarg);
			break;
		case 'o':
			output_filename = string(optarg);
			break;
		case 'm':
			mode = stoi(optarg);
			if (mode < 0 || mode > 2) {
				cout << "Error: mode must be 0 (q filtering), 1 (rt/mz/int projection), or 2 (int sorting)\n";
			}
			break;
		case 'a':
			min_rt = stod(optarg);
			break;
		case 'b':
			max_rt = stod(optarg);
			break;
		default:
			cout << "Usage: csvFilter -i <input csv> -o <output csv> -m <mode> [-a <min rt> -b <max rt>]\n";
			exit(1);
		}
	}

	// Redirect standard IO to the correct csv files
	freopen(input_filename, "r", stdin);
	freopen(output_filename, "w", stdout);

	string input;
	cin >> input;

	// Get the headers
	unordered_map<string, string> headers;
	int idx = 0;

	stringstream ss1(input);
	while (ss1.good()) {
		getline(ss1, input, ',');
		table[0].push_back(input);
		headers.insert({ input, idx++ });
	}

	// Projection only requires these three columns
	if (mode == 1) {
		vector<string> header { "RT", "m/z", "Intensity" };
		table.push_back(header);
	}

	// Simultaneously read and filter the input
	for (int i = 0; cin >> input; i++) {
		vector<string> line;
		stringstream ss2(input);
		while (ss2.good()) {
			getline(ss2, input, ',');
			line.push_back(input);
		}

		if (mode == 0) {
			double rt = stod(line[headers.find("RT")->second]);
			double q = stod(line[headers.find("q")->second]);

			if (rt >= min_rt && rt <= max_rt && q >= Q_THRESHOLD) {
				table.push_back(line);
			}
		}
		else if (mode == 1) {
			string rt = line[headers.find("RT")->second];
			string mz = line[headers.find("m/z")->second];
			string intensity = line[headers.find("Intensity")->second];

			vector<string> projection { rt, mz, intensity };
			table.push_back(projection);
		}
		else {
			table.push_back(line);
		}
	}

	idx = mode == 1 ? 1 : 0;
	for (int i = 0; i < table[idx].size(); i++) {
		cout << table[idx][i];
		if (i < table[idx].size() - 1) {
			cout << ',';
		}
		cout << '\n';
	}

	if (mode == 2) {
		vector<pair_dvs> sorted;
		for (int i = 1; i < table.size(); i++) {
			double intensity = stod(line[headers.find("Intensity")->second]);
			sorted.push_back({ intensity, &table[i] });
		}

		sort(sorted.begin(), sorted.end(), [](const pair_dvs &a, const pair_dvs &b) {
			return a.first >= b.first;
		});

		for (int i = 0; i < sorted.size(); i++) {
			vector<string> *temp = sorted[i].second;
			for (int j = 0; j < temp->size(); j++) {
				cout << (*temp)[j];
				if (j < temp->size() - 1) {
					cout << ',';
				}
			}
			cout << '\n';
		}
	}
	else {
		for (int i = idx + 1; i < table.size(); i++) {
			for (int j = 0; j < table[i].size(); j++) {
				cout << table[i][j];
				if (j < table[i].size() - 1) {
					cout << ',';
				}
			}
			cout << '\n';
		}
	}

	return 0;
}
