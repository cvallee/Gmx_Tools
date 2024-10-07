'''
GROMACS .xvg files tools
'''
import numpy as np
import os

class XVG:
    '''Convert an .xvg file into a python object'''
    def __init__(self, file, parse=True):
        self.filepath = os.path.abspath(file)
        if '\\' in self.filepath:
            self.filename = self.filepath.split('\\')[-1]
        else:
            self.filename = self.filepath.split('/')[-1]
        self.title = ''
        self.x_column = []
        self.x_label = ''
        self.y_columns = []
        self.y_label = ''
        self.legends = []
        if parse:
            self.parse_file()
    
    def parse_file(self):
        if self.filename.split('.')[-1] != 'xvg':
            raise ValueError('ERROR! File provided is not an XVG file')
        
        with open(self.filepath, 'r') as f:
            data = {}
            for line in f:
                # Ignore comments
                if line[0] == '#':
                    continue
                # Store title and x/y labels
                if line[0] == '@':
                    if line.split()[1] == 'title':
                        self.title = line.split('"')[-2]
                    if line.split()[1] == 'xaxis':
                        self.x_label = line.split('"')[-2]
                    if line.split()[1] == 'yaxis':
                        self.y_label = line.split('"')[-2]
                    if len(line.split()) > 2:
                        if line.split()[2] == 'legend':
                            self.legends.append(line.split('"')[-2])
                    continue
                
                l = line.split()
                for i, value in enumerate(l):
                    if i in data.keys():
                        data[i].append(float(value))
                    else:
                        data[i] = [float(value)]
                # Extract x and y
                self.x_column = data[0]
                self.y_columns = [data[i] for i in range(1,len(data))]
                    
    def get3Dcoord(self, np_array=False):
        '''Return a list of N 3D coordinates for each time step'''
        if not self.x_column and not self.y_columns:
            self.parse_file()
        assert len(self.y_columns) % 3 == 0
        coordinates = []
        # Loop over each time step
        for i, time in enumerate(self.x_column):
            coord = []
            n = 0
            # Loop over all X, Y and Z
            for j in range(int(len(self.y_columns)/3)):
                # Get the X, Y, Z  from the three following column
                x = self.y_columns[n][i]
                n += 1
                y = self.y_columns[n][i]
                n += 1
                z = self.y_columns[n][i]
                n += 1
                # Append (X, Y, Z) to coord list
                coord.append((float(x),float(y),float(z)))
            # Append list of N coord for each time step
            if np_array:
                coordinates.append(np.asarray(coord))
            else:
                coordinates.append(coord)
            
        return coordinates

    def plot(self, title='', x_label='', y_label='', legends=[], show=False, save=False, **kwargs):
        if not self.x_column and not self.y_columns:
            self.parse_file()
        import matplotlib.pyplot as plt

        # Give the opportunity for the user to overwrite the title, x and y labels, and the legends
        if title:
            plot_title = title
        else:
            plot_title = self.title
        if x_label:
            plot_x_label = x_label
        else:
            plot_x_label = self.x_label
        if y_label:
            plot_y_label = y_label
        else:
            plot_y_label = self.y_label
        if legends:
            plot_legends = legends
        else:
            plot_legends = self.legends

        assert len(self.y_columns) == len(self.legends)

        fig, ax = plt.subplots()
        for y, leg in zip(self.y_columns, plot_legends):
            ax.plot(self.x_column, y, label=leg, **kwargs)
            ax.title.set_text(plot_title)
            ax.set_xlabel(plot_x_label)
            ax.set_ylabel(plot_y_label)
            ax.legend()
        if show:
            plt.show()
        if save:
            if bool(save):
                plt.savefig(self.filepath.replace('.xvg','.png'))
            if str(save):
                plt.savefig(save)
        
        return fig, ax


        