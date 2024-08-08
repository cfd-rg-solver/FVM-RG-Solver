#from fileinput import filename
import matplotlib.pyplot as plt
# import numpy as np
import pandas as pd
import os

# Matplotlib Style:

plt.style.use('style.mpltstyle')

def get_file_extension(file_path):
    _, file_extension = os.path.splitext(file_path)
    return file_extension.lower()

class DataPlotter:
    def __init__(self, data_file, path_out, file_name):
        self.data_file = data_file
        self.path_out = path_out
        self.file_name = file_name

    def load_data(self, str_sep, legend_names=None, num_skip_row=0, num_read_rows=None, \
                  sheet_excel=0):
        try:
            type = get_file_extension( self.data_file)
            if type in ('.txt', '.csv'):
                self.data = pd.read_csv(self.data_file, sep = str_sep, header=0, names = legend_names, \
                                        index_col=0, skiprows = num_skip_row, nrows = num_read_rows) # header=0
            elif type in ('.xlsx', '.xls'):
                self.data = pd.read_excel(self.data_file, sheet_name=sheet_excel, index_col=0, header=0, \
                                          skiprows = num_skip_row, nrows = num_read_rows)
            else:
                print(f"File type '{type}' not found.")
            self.df = pd.DataFrame(self.data)
        except FileNotFoundError:
            print(f"File '{self.data_file}' not found.")
            self.data = None

    def plot_data(self, ylabel):
        if self.data is not None:
            plt.figure()

            self.df.plot(legend=None) # meanwhile legend is not required
            # plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            
            plt.savefig(self.path_out + self.file_name)
            #plt.show()
        else:
            print("Data not loaded. Please load data first.")

if __name__ == "__main__":

    current_directory = os.getcwd()
    parent_directory = os.path.dirname(current_directory)
    path_data = parent_directory + '\\toPlot\\'
    path_plot = parent_directory + '\\fig\\'

    files = os.listdir(path_data) 

    # legend_names = []
    plot_names = ['Density, kg/m^3', 'Sp. heat ratio', 'Pressure, Pa', 'Temperature, K', 'Temperature-vibr, K', 'Velocity, m/s', 'Velocity normal, m/s', 'Velocity tangent, m/s']

    for i in range(len(files)):
        
        plotter = DataPlotter(path_data + files[i], path_plot, files[i].replace('.txt',''))
        plotter.load_data(' ')   
        plotter.plot_data(plot_names[i])

    

