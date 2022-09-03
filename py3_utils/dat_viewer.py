import sys
import math
import pandas as pd
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
# group  = parser.add_mutually_exclusive_group()
parser.add_argument ("--path", default="/home/wdklotz/pyOrbit/examples+/SNS_Linac/pyorbit_linac_model", help="path to data source")
parser.add_argument ("--file", default="pyorbit_twiss_sizes_ekin.dat", help="data source file")
parser.add_argument ("--x", default="position", help="Abscissa (def = position")
parser.add_argument ("--y1", default="betaX", help="Ordinate (def = betaX")
parser.add_argument ("--y2", default="", help="Ordinate (def = '')")
parser.add_argument ("--nh", action="store_true", help="no header (def = False)")
args = vars(parser.parse_args())

file_dir   = args['path']
file_name  = args['file']
qfile_name = f'{file_dir}/{file_name}'
print(f'data source: {qfile_name}')

header = not args["nh"]

if header:
    df = pd.read_csv(qfile_name, sep="\s+")
    print(df.info())
    x_column  = args['x']
    y1_column = args['y1']
    y2_column = args['y2']
    y_columns = [y1_column]
    if y2_column != "": y_columns.append(y2_column)
    print('====================================================================')
    print(f'plotting ordinates "{y_columns}" against abscissa "{x_column}"')
    print('====================================================================')
else:
    df = pd.read_csv(qfile_name, header=None, sep="\s+")
    print(df.info())
    x_column  = int(args['x'])
    y1_column = int(args['y1'])
    y2_column = args['y2']
    print(y2_column)
    y_columns = [y1_column]
    if y2_column != "": y_columns.append(int(y2_column))
    print('====================================================================')
    print(f'plotting ordinates "{y_columns}" against abscissa "{x_column}"')
    print('====================================================================')

df.plot(x_column,y_columns,linestyle='-',linewidth="0.6", title=qfile_name)
plt.show()
