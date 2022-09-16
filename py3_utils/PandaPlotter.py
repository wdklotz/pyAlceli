import sys
import math
import pandas as pd
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse

class PandaPlotter(object):
    def __init__(self,args):
        self.dir        = args.get('path')
        self.file       = args.get('file')
        self.header     = args.get('nh')
        self.abscissa   = args.get('x')
        self.ordinate1  = args.get('y1')
        self.ordinate2  = args.get('y2')
        self.ordinate3  = args.get('y3')
        self.ordinate4  = args.get('y4')
    def do_plot(self):
        qfile_name = f'{self.dir}/{self.file}'
        print(f'data source: {qfile_name}')

        header    = not self.header
        if header:
            df = pd.read_csv(qfile_name, sep="\s+")
        else:
            df = pd.read_csv(qfile_name, header=None, sep="\s+")
        df.info()

        x_column  = self.abscissa
        y1_column = self.ordinate1
        y2_column = self.ordinate2
        y3_column = self.ordinate3
        y4_column = self.ordinate4
        y_columns = [y1_column]
        if y2_column != "": y_columns.append(y2_column)
        if y3_column != "": y_columns.append(y3_column)
        if y4_column != "": y_columns.append(y4_column)
        print('====================================================================')
        print(f'plotting ordinates "{y_columns}" against abscissa "{x_column}"')
        print('====================================================================')

        df.plot(x_column,y_columns,linestyle='-',linewidth="0.6", title=qfile_name)
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # group  = parser.add_mutually_exclusive_group()
    parser.add_argument ("--path", default="/home/wdklotz/PyOrbit/examples+/SNS_Linac/pyorbit_linac_model", help="path to data source")
    parser.add_argument ("--file", default="pyorbit_twiss_sizes_ekin.dat", help="data source file")
    parser.add_argument ("--x",    default="position", help="Abscissa (def = position)")
    parser.add_argument ("--y1",   default="betaX", help="Ordinate (def = betaX)")
    parser.add_argument ("--y2",   default="", help="2nd Ordinate")
    parser.add_argument ("--y3",   default="", help="3rd Ordinate")
    parser.add_argument ("--y4",   default="", help="4th Ordinate")
    parser.add_argument ("--nh",   action="store_true", help="no header (def = False)")
    args = vars(parser.parse_args())
    # print(args)

    plotter = PandaPlotter(args)
    plotter.do_plot()

