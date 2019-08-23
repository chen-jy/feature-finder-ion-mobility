#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

vector<vector<string>> table(1), filtered;
vector<string> empty_vec;

int main(int argc, char *argv[]) {
	cin.sync_with_stdio(0);
	cin.tie(0);

	// Temporary arg processing until Win 10 gets its full Linux kernel
	if (argc < 3) {
		cout << "Usage: csvFilter <csv filename> <mode (0 or 1)>\n";
		exit(1);
	}

	int mode = atoi(argv[2]);
	if (mode != 0 && mode != 1) {
		cout << "Error: mode must be 0 (q filtering) or 1 (rt/mz/int filtering)\n";
		exit(1);
	}

	freopen(argv[1], "r", stdin);

	string line, entry;
	cin >> line;

	// Get headers
	stringstream ss1(line);
	while (ss1.good()) {
		getline(ss1, line, ',');
		table[0].push_back(line);
	}

	// Read in the entire csv first
	// TODO: can read and filter in one step
	int i = 1;
	while (cin >> line) {
		table.push_back(empty_vec);
		stringstream ss(line);
		while (ss.good()) {
			getline(ss, line, ',');
			table[i].push_back(line);
		}
		i++;
	}

	//cout << "Rows: " << table.size() << "\n";

	// Target column indices
	int idx_rt = distance(table[0].begin(), find(table[0].begin(), table[0].end(), "RT"));
	int idx_mz = distance(table[0].begin(), find(table[0].begin(), table[0].end(), "m/z"));
	int idx_int = distance(table[0].begin(), find(table[0].begin(), table[0].end(), "Intensity"));
	int idx_q = distance(table[0].begin(), find(table[0].begin(), table[0].end(), "initialPeakQuality"));

	for (i = 1; i < table.size(); i++) {
		double rt = atof(table[i][idx_rt].c_str());
		double mz = atof(table[i][idx_mz].c_str());
		double intensity = atof(table[i][idx_int].c_str());
		double q = atof(table[i][idx_q].c_str());

		// The RT range must be changed for new data
		if (mode == 0 && 3000 <= rt && rt <= 3100 && q >= 0.01) {
			filtered.push_back(table[i]);
		}
		else if (mode == 1) {
			vector<string> row;
			row.push_back(to_string(rt));
			row.push_back(to_string(mz));
			row.push_back(to_string(intensity));
			filtered.push_back(row);
		}
	}

	//cout << "Results: " << filtered.size() << "\n";

	// Write results out to a new csv file
	freopen("filtered.csv", "w", stdout);

	if (mode == 0) {
		for (i = 0; i < table[0].size(); i++) {
			cout << table[0][i];
			if (i < table[0].size() - 1)
				cout << ",";
		}
		cout << "\n";
	}
	else if (mode == 1) {
		cout << "RT,m/z,Intensity\n";
	}

	for (i = 0; i < filtered.size(); i++) {
		for (int j = 0; j < filtered[i].size(); j++) {
			cout << filtered[i][j];
			if (j < filtered[i].size() - 1)
				cout << ",";
		}
		cout << "\n";
	}

	return 0;
}
