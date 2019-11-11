import argparse
import csv
from operator import itemgetter

def checkFloat(val):
    try:
        return float(val)
    except ValueError:
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Projection')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)
    parser.add_argument('--start', action='store', required=True, type=float)
    parser.add_argument('--stop', action='store', required=True, type=float)
    args = parser.parse_args()
    
    data = []
    
    with open(args.input, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            rt = checkFloat(row['Retention time'])
            mz = checkFloat(row['m/z'])
            intensity = checkFloat(row['Intensity'])
            
            if rt is False or mz is False or intensity is False:
                continue
            if mz <= 10 or intensity <= 10:
                continue
            if args.start <= rt <= args.stop:
                data.append([rt, mz, intensity])
                
    data = sorted(data, key=itemgetter(2))
            
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    