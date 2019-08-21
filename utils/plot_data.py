import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Reads an input file where each line has three floats, representing a point in three-
# dimensional space, and plots all of the points in a scatter plot. The output can be
# saved if required, and the axes renamed.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot 3D scatter data.')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=False, type=str)
    parser.add_argument('--xlabel', action='store', required=False, type=str)
    parser.add_argument('--ylabel', action='store', required=False, type=str)
    parser.add_argument('--zlabel', action='store', required=False, type=str)

    args = parser.parse_args()

    with open(args.input) as f:
        input = f.readlines()
    input = [x.strip() for x in input]

    data = []
    for line in input:
        x = line.split()
        data.append([float(y) for y in x])

    transpose = list(zip(*data))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(list(transpose[0]), list(transpose[1]), list(transpose[2]))

    xlabel = 'RT threshold' if args.xlabel is None else args.xlabel
    ylabel = 'm/z threshold' if args.ylabel is None else args.ylabel
    zlabel = '# of features' if args.zlabel is None else args.zlabel

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    if args.output is None:
        plt.show()
    else:
        plt.savefig(args.output)
