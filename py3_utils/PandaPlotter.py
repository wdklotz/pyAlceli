import sys, os, math, argparse
import pandas as pd
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

class PandaPlotter(object):
    def __init__(self,args):
        self.dir        = args.get('px')
        self.file       = args.get('file')
        self.header     = args.get('nh')
        self.abscissa   = args.get('x')
        self.ordinate1  = args.get('y')
        self.ordinate2  = args.get('y2')
        self.ordinate3  = args.get('y3')
        self.ordinate4  = args.get('y4')
    def do_plot(self):
        qfile_name = f'{self.dir}/{self.file}' if self.dir != "" else f'{self.file}'
        qfile_name = os.path.normpath(qfile_name)
        print(f'data source: {qfile_name}')

        header    = not self.header
        if header:
            df = pd.read_csv(qfile_name, sep="\s+")
        else:
            df = pd.read_csv(qfile_name, header=None, sep="\s+")

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

        try:
            df.plot(x_column,y_columns,linestyle='-',linewidth="0.6", title=qfile_name)
            plt.show()
        finally:
            df.info()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument ("--px",   default="",          help="prefix to data source")
    parser.add_argument ("file",   default="data.txt",  help="data source file")
    parser.add_argument ("x",      default="position",  help="Abscissa")
    parser.add_argument ("y",      default="betaX",     help="Ordinate")
    parser.add_argument ("--y2",   default="",          help="2nd Ordinate")
    parser.add_argument ("--y3",   default="",          help="3rd Ordinate")
    parser.add_argument ("--y4",   default="",          help="4th Ordinate")
    parser.add_argument ("--nh",   action="store_true", help="no header (def = False)")
    args = vars(parser.parse_args())
    # print(args)
    PandaPlotter(args).do_plot()

