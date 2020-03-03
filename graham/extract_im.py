import csv
import diapysef.timsdata as tims
import numpy as np
import pyopenms as ms

td = tims.TimsData('raw-files/20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_2_A1_01_2768.d')
conn = td.conn

q = conn.execute('SELECT COUNT(*) FROM Frames')
row = q.fetchone()
N = row[0]

scan_number_axis = np.arange(N, dtype=np.float64)
ook0_axis = td.scanNumToOneOverK0(1, scan_number_axis)

lines = []
with open('temp/evidence-projected.csv', 'r') as file:
    reader = csv.reader(file)
    for line in reader:
        lines.append(line)

for i in range(len(lines)):
    lines[i][2] = ook0_axis[int(lines[i][2])]

with open('temp/evidence-projected-im.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(lines)
