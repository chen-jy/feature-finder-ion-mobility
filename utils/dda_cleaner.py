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
    parser.add_argument('--file', action='store', required=True, type=str)
    parser.add_argument('--start', action='store', required=True, type=float)
    parser.add_argument('--stop', action='store', required=True, type=float)
    parser.add_argument('--pep', action='store', required=False, type=float)
    
    args = parser.parse_args()
    PEP = args.pep if args.pep is not None else float('inf')
    
    data = []
    errorFalse, errorSmall = 0, 0
    cleanRows = 0
    
    with open(args.input, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            rf = row['Raw file']
            rt = checkFloat(row['Retention time'])
            mz = checkFloat(row['m/z'])
            intensity = checkFloat(row['Intensity'])
            pep = checkFloat(row['PEP'])
            
            if rf != args.file:
                continue
            if rt is False or mz is False or intensity is False:
                errorFalse += 1
                continue
            if mz <= 10 or intensity <= 10:
                errorSmall += 1
                continue
            if args.start <= rt <= args.stop and pep < PEP:
                data.append([rt, mz, intensity])
                cleanRows += 1
                
    data = sorted(data, key=itemgetter(2))
            
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    
    print("'False' errors:", errorFalse)
    print("'Small' errors:", errorSmall)
    print("Total rows:", cleanRows)
    