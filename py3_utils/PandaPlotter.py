import sys, os, math, argparse
import pandas as pd
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

class PandaPlotter(object):
    def __init__(self,args):
        # print(args)
        self.dir        = args.get('px')
        self.file       = args.get('file')
        self.header     = args.get('nh')
        self.abscissa   = None
        self.ordinate1  = None
        self.ordinate2  = None
        self.ordinate3  = None
        self.ordinate4  = None
        qfile_name = f'{self.dir}/{self.file}' 
        self.qfile_name = os.path.normpath(qfile_name)
        # print(f'data source: {qfile_name}')
        header    = not self.header
        if header:
            self.df = pd.read_csv(qfile_name, sep="\s+")
        else:
            self.df = pd.read_csv(qfile_name, header=None, sep="\s+")

    def do_plot(self,args):
        self.abscissa   = args.get('x')
        self.ordinate1  = args.get('y')
        self.ordinate2  = args.get('y2')
        self.ordinate3  = args.get('y3')
        self.ordinate4  = args.get('y4')

        x_column        = self.abscissa
        y1_column       = self.ordinate1
        y2_column       = self.ordinate2
        y3_column       = self.ordinate3
        y4_column       = self.ordinate4
        y_columns       = [y1_column]
        if y2_column != None: y_columns.append(y2_column)
        if y3_column != None: y_columns.append(y3_column)
        if y4_column != None: y_columns.append(y4_column)
        print('================================================================================================================')
        print(f'plotting ordinates "{y_columns}" against abscissa "{x_column}"')
        print('================================================================================================================')

        self.df.plot(x_column,y_columns,linestyle='-',linewidth="0.6", title=self.qfile_name)
        plt.show()
    def data_info(self):
        self.df.info()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument ("--px",   default="",          help="prefix to data source")
    parser.add_argument ("file",   default="data.txt",  help="data source file")
    parser.add_argument ("x",      default="position",  help="Abscissa")
    parser.add_argument ("y",      default="betaX",     help="Ordinate")
    parser.add_argument ("--y2",   default=None,        help="2nd Ordinate")
    parser.add_argument ("--y3",   default=None,        help="3rd Ordinate")
    parser.add_argument ("--y4",   default=None,        help="4th Ordinate")
    parser.add_argument ("--nh",   action="store_true", help="[False] no header")
    parser.add_argument ("--ex",   action="store_true", help="[False] cli command example")
    args = vars(parser.parse_args())
    if args['ex']:
        print('CLI example => python PandaPlotter.py --px ".." pyorbit_twiss_sizes_ekin.dat position emittX --y2 emittZ')
    else:
        PandaPlotter(args).do_plot(args)

