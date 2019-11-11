import argparse
import csv

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Projection')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)
    args = parser.parse_args()
    
    data = []
    
    with open(args.input, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            rt = row['Retention time']
            mz = row['m/z']
            intensity = row['Intensity']
            data.append([rt, mz, intensity])
            
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    