import argparse
import csv
import matplotlib.pyplot as plt
import statistics as stat

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Histogram')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--bins', action='store', required=True, type=int)
    args = parser.parse_args()
    
    data = []
    
    with open(args.input, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pep_str = row['PEP']
            if pep_str == '' or pep_str == 'NaN':
                pep = float(2)
            else:
                pep = float(pep_str)
                
            if pep < 1:
                data.append(pep)
            
    
    n, bins, patches = plt.hist(data, args.bins)
    plt.show()
    
    print("min:", min(data))
    print("max:", max(data))
    print("mean:", stat.mean(data))
    print("mode:", stat.mode(data))
    