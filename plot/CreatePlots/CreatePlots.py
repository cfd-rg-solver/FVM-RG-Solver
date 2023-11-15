#from fileinput import filename
import matplotlib.pyplot as plt
# import numpy as np
import pandas as pd
import os

# Matplotlib Style:

plt.style.use('style.mpltstyle')


class DataPlotter:
    def __init__(self, data_file, path, file_name):
        self.data_file = data_file
        self.path = path
        self.file_name = file_name

    def load_data(self):
        try:
            self.data = pd.read_csv(self.data_file, sep=' ', header=None)  
            self.df = pd.DataFrame(self.data)
        except FileNotFoundError:
            print(f"File '{self.data_file}' not found.")
            self.data = None

    def plot_data(self):
        if self.data is not None:
            plt.figure()

            x = self.df[0]
            y = self.df[1]

            #df.plot()
            plt.plot(x,y)
            file_name = self.file_name.replace('.txt','')
            plt.ylabel(file_name)
            plt.xlabel('y')
            
            plt.savefig(self.path + file_name)
            #plt.show()
        else:
            print("Data not loaded. Please load data first.")

if __name__ == "__main__":

    current_directory = os.getcwd()
    parent_directory = os.path.dirname(current_directory)
    path_data = parent_directory + '\\toPlot\\'
    path_plot = parent_directory + '\\fig\\'

    plotName = os.listdir(path_data) 

    for name in plotName:
        
        plotter = DataPlotter(path_data + name, path_plot, name)
        plotter.load_data()
        plotter.plot_data()
    

