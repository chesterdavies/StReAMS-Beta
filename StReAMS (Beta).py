"""

@author: Chester Davies

Title: StReAMS (Statistical Reservoir and Model Significance)

"""

####################################################################################################################
##Pip Installs##
####################################################################################################################
# !pip install tkfilebrowser
# !pip install numpy
# !pip install astroML
# !pip install astropy
# !pip install scipy
# !pip install scikit-gstat
# !pip install seaborn
# !pip install pandas

####################################################################################################################
##Package Imports##
####################################################################################################################
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
import scipy.stats as st

from tkinter import*
from statistics import mean
from scipy import signal
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
from astropy.timeseries import LombScargle
from skgstat.models import spherical
from scipy.ndimage import gaussian_filter1d
from astroML.time_series import lomb_scargle, search_frequencies

# Importing tkinter packages for use with the UI
try:
    import tkinter as tk
    import tkinter.ttk as ttk
    from tkinter import filedialog
    from tkinter import messagebox
except ImportError:
    import ttk
    import tkFileDialog as filedialog

####################################################################################################################
##UI Code##
####################################################################################################################

window = Tk()

window.title('Reservoir Modelling Statistical Tool')
window.geometry("550x300+10+20")

root = tk.Tk()
root.overrideredirect(True)
root.geometry('0x0+0+0')
root.focus_force()
ttl = 'Select File'

#Defining the variable to hold the array of all target files to be read in from the target folder         
targetFiles = []
global targetFolder
global destinationFolder

#Module to clear all of the tickboxes in the UI form, resetting the inputs
def clearForm():
    v5.set(0)
    v6.set(0)
    
#Checks the tickbox values and determines whether the input parameters have been satisfied
#If at least 1 statistical method is selected, and a target and destination folder is selected, then allows for programme to run
#If not, then an error is thrown and the user is prompted to rectify the issue
def startProcess():
    if (((v5.get() == 0) or (v6.get() == 0))):
        if (v5.get() == 0):
            messagebox.showinfo("Error","There is no target folder selected")
        if (v6.get() == 0):
            messagebox.showinfo("Error","There is no destination folder selected")
    else:
        #Placeholder
        messagebox.showinfo("Information","Brill")
        fileReader(targetFolder, targetFiles, destinationFolder)
        dataSorter()
        graphingFunction()

####################################################################################################################
##File Reader Code##
####################################################################################################################        

#Reads in the full extent of the files held within the targetFiles array, trims any excess lines, compiles and formats the data into a single file
def fileReader(targetFolder, targetFiles, destinationFolder):
    
    data_list = []
    modelSet = 1
    
    modelRun = []
    
    #Opens each file incrementally, and replaces the entire data list with only those that relate to a facies code (integer)
    for targetFile in targetFiles[0]:
        
        #Setting the file being used to the file path and file name
        print(targetFile)
        file = targetFolder + '/' + str(targetFile)
        print(file)
        
        #Opens the file and removes all whitespace
        readFile = pd.read_csv(file, delim_whitespace = 'true')
        
        print(readFile)
        
        #Saves only the rows with facies codes to the original variable
        readFile = readFile[((readFile['Statistics']== '0') | (readFile['Statistics']== '1') | (readFile['Statistics'] == '3') | (readFile['Statistics']== '4') | (readFile['Statistics']== '5'))]

        print(readFile)

        print(readFile.columns)
        
        #Adding model set and file name columns to the data held within the file for easier data management
        readFile.insert(0, "Model Set", modelSet)
        readFile.insert(10, "File Name", targetFile)
        
        print(readFile['Statistics'].nunique())
        
        #Finding the number of rows present within the file to then assign the model set and model run
        count = readFile['Statistics'].size
        
        count = int(count/readFile['Statistics'].nunique())
        
        print(count)
        
        #Adds the iteration number for the run number into the data
        for i in range(1, count+1):
            for count in range(readFile['Statistics'].nunique()):    
                modelRun.append(i)
            
        print(modelRun)

        #Adding the model run column to the dataset        
        readFile.insert(1, "Model Run", modelRun)
        
        data_list.append(readFile)
        
        #Appending the values for the model sets into the dataset
        modelRun = []
        modelSet += 1
        
        print(modelSet)
           
    #Appends the new data into the array for the compiled array
    compiledData = pd.concat(data_list)
    
    global file1
    file1 = destinationFolder + '/Combined.csv'
    
    path = destinationFolder + '/Images'
    
    try: 
        os.mkdir(path) 
    except OSError as error: 
        print(error)
    
    #Creates and opens the compiled file in the destination folder
    with open(file1, mode='w', newline='') as newTestFile:        
        
        #Renaming the headings for the file
        compiledData.rename(columns = {'Statistics':'Facies Code', 'for':'Name', 'Facies':'Percentage', 'in':'N', 'zone':'Intervals', 'Top':'Min', '-':'Mean', 'Base':'Max', '(Unfiltered)':'Standard Deviation'}, inplace = True)
        
        #Appending the freshly cleaned data into the newly created file
        compiledData.to_csv(file1, index = False)
        
        print(compiledData.columns)
    
####################################################################################################################
##Plotting Graphs Code##
####################################################################################################################

#Function to sort the data from the combined csv file into seperate arrays based on the facies code and model set
def dataSorter():

    globals()['identifierArray'] = []
    globals()['combinedIDArray'] = []
    
    #Reads in the data from the selected CSV file
    df = pd.read_csv(file1, sep = ',')
        
    print(df)
        
    print(df['Model Run'].max())
    print(df['Model Run'].min())
    
    print(df['Model Set'].max())  
    print(df['Model Set'].min())
    
    #Determines the unique facies codes to split the sand and shale values
    faciesCodes = df['Facies Code'].unique()
    globals()['modelSets'] = df['Model Set'].unique()
    globals()['fileNames'] = df['File Name'].unique()
          
    #Creating the data capture array for regular values
    columnHeadingsVals = {'Identifier', 'File Name', 'Maximum GERG CI', 'Maximum PERG CI', 'Mean PERG CI', 'Mean GERG CI', 'STDEV PERG CI', 'STDEV GERG CI', 'TF PERG CI', 'TF GERG CI', 'Maximum GERG SIG', 'Maximum PERG SIG', 'Mean PERG SIG', 'Mean GERG SIG', 'STDEV PERG SIG', 'STDEV GERG SIG', 'TF PERG SIG', 'TF GERG SIG', 'Mean Periodicity1', 'Maximum Periodicity1', 'STDEV Periodicity1', 'TF Periodicity1', 'Mean Periodicity5', 'Maximum Periodicity5', 'STDEV Periodicity5', 'TF Periodicity5', 'AVG Ergodic', 'MAX Ergodic', 'STDEV Ergodic', 'TF Ergodic', 'Mean CI', 'Maximum CI', 'STDEV CI', 'TF CI', 'Mean Avg', 'Maximum Avg', 'STDEV Avg', 'TF Avg', 'Mean Min', 'Maximum Min', 'STDEV Min', 'TF Min', 'Identifier', 'Mean Max', 'Maximum Max', 'STDEV Max', 'TF Max', 'Mean Range', 'Maximum Range', 'STDEV Range', 'TF Range', 'Mean Nugget', 'Maximum Nugget', 'STDEV Nugget', 'TF Nugget', 'Mean Sill', 'Maximum Sill', 'STDEV Sill', 'TF Sill'}
    columnHeadingsPlots = {'Identifier', 'AVG BW Plot', 'MAX BW Plot', 'STDEV BW Plot', 'TF BW Plot', 'AVG LINE Plot', 'MAX LINE Plot', 'STDEV LINE Plot', 'TF LINE Plot', 'AVG HST Plot', 'MAX HST Plot', 'STDEV HST Plot', 'TF HST Plot', 'AVG GERG HST Plot', 'MAX GERG HST Plot', 'STDEV GERG HST Plot', 'TF GERG HST Plot', 'AVG PERG HST Plot', 'MAX PERG HST Plot', 'STDEV PERG HST Plot', 'TF PERG HST Plot', 'AVG SV Plot', 'MAX SV Plot', 'STDEV SV Plot', 'TF SV Plot', 'Identifier', 'AVG PDG Plot', 'MAX PDG Plot', 'STDEV PDG Plot', 'TF PDG Plot', 'AVG Goovaerts Plot', 'MAX Goovaerts Plot', 'STDEV Goovaerts Plot', 'TF Goovaerts Plot'}
    
    globals()['INDIVIDUALVALS'] = pd.DataFrame(columns = columnHeadingsVals)
    globals()['INDIVIDUALPLOTS'] = pd.DataFrame(columns = columnHeadingsPlots)
    
    #For loop for collating a list of identifiers - with a general composition of 'facies number'.'model set'
    for facies in faciesCodes:
        print(facies)
        for files in globals()['fileNames']:
            fileN = files[:-4]
            identifier = str(facies) + '.' + fileN
            globals()['data%s' % identifier] = []
            for index, row in df.iterrows():               
                if ((row['Facies Code'] == facies) and (row['File Name'] == files)):
                    print(row)
                    globals()['data%s' % identifier].append(row)
            identifierArray.append(identifier)
            
            #Creation of the path for the images folder to be generated into
            path = destinationFolder + '/Images/' + identifier
            
            #Attempts to create the folder in the selected path, if this is not possible, then an error is printed
            try: 
                os.mkdir(path) 
            except OSError as error: 
                print(error)
     
    #Creating two blank two-dimensional arrays with only the identifier and file name, ready for the output values and graph locations to be assigned to
    for i in identifierArray:
        globals()['INDIVIDUALVALS'] = globals()['INDIVIDUALVALS'].append(pd.Series(), ignore_index = True)
        globals()['INDIVIDUALVALS'].loc[[len(globals()['INDIVIDUALVALS'])-1],'Identifier'] = i
        globals()['INDIVIDUALVALS'].loc[[len(globals()['INDIVIDUALVALS'])-1],'File Name'] = globals()['data%s' % i][0]['File Name']
        print(globals()['INDIVIDUALVALS']['File Name']) 
        globals()['INDIVIDUALPLOTS'] = globals()['INDIVIDUALPLOTS'].append(pd.Series(), ignore_index = True)
        globals()['INDIVIDUALPLOTS'].loc[[len(globals()['INDIVIDUALPLOTS'])-1],'Identifier'] = i
        print(globals()['INDIVIDUALPLOTS']['Identifier']) 
        
    print('END PLOT')
    
####################################################################################################################
##Histogram and Box Plot Code##
####################################################################################################################
#Function for the creation of the histograms to look at the distribution of the dataset
def histogramBoxPlot(dataArray, dataVal, i, path):
    
    print(dataArray)
    
    histogram, (ax_Histogram) = plt.subplots(1, sharex=True)  
    
    ax_Histogram = sns.distplot(x=dataArray, kde=True, hist_kws=dict(edgecolor='b', linewidth=0))
            
    ax_Histogram.set(ylabel = 'KDE - Kernel Density Estimation')
            
    if dataVal == 0:
        ax_Histogram.set(xlabel = 'Thickness (m)')
        plt.title('Average Thickness Histogram (m)')
    elif dataVal == 1:
        ax_Histogram.set(xlabel = 'Standard Deviation of Thickness')
        plt.title('Standard Deviation of Thickness Histogram')
    elif dataVal == 2:
        ax_Histogram.set(xlabel = 'Target Fraction (%)')
        plt.title('Values of Modelled Target Fraction Histogram (%)')
                          
    print("FC" + str(i) + " 95% CI")
    print(st.t.interval(alpha=0.95, df=len(dataArray)-1, loc=np.mean(dataArray), scale=st.sem(dataArray)))
    
    plt.tight_layout()
    
    plt.show()
            
    histogramPath = path + 'Histogram.png'
    
    plt.savefig(histogramPath)	
                
    plt.clf()
    
    plt.close()

####################################################################################################################
##Periodogram Code##
####################################################################################################################
#Function for the creation of the periodogram graphs to determine the PERG value
def periodogramPlot(dataArray, dataVal, i, path):

    periodogram, (ax_Periodogram) = plt.subplots(1, sharex=True)             
                     
    t = list(range(1, len(dataArray)+1))
    
    #Detrending the signal
    detrended = signal.detrend(dataArray)
        
    print(i)
    
    #One-dimensional gaussian filter to smooth the dataset           
    smoothed = (gaussian_filter1d(detrended, sigma=2))
                   
    df = pd.DataFrame(smoothed)
    Y = df[0]
        
    dy = 0.5 + 0.5 * np.random.random(len(dataArray))
    omega = np.linspace(50, 100, 500)
    
    #Adding the 5%, 1% and 0.5% significance markers to the periodogram    
    sig = np.array([0.05, 0.01, 0.005])
    PS, z = lomb_scargle(df.index, Y, dy, omega, generalized=True, significance=sig)
        
    omega = omega/ (2* len(dataArray))
        
    print(omega)
    
    frequency, power = LombScargle(t, smoothed).autopower(minimum_frequency = 0, maximum_frequency = 0.5)
    ax_Periodogram.plot(frequency, power) 
    
    xlim = (0, 0.5)
    for zi, pi in zip(z, sig):
        plt.plot(xlim, (zi, zi), ':k', lw=1)
        plt.text(xlim[-1] - 0.001, zi - 0.02, "$%.1g$" % pi, ha='right', va='top')
        
    ax_Periodogram.legend(labels=['Detrended and Smoothed Data'], facecolor='white')
    
    print(power)
    
    maxPower = 0
    PERG = 0
    
    #Locating the peak power within the singal and then determining its position in relation to the number of realisations through 1/power
    for p in range(1, len(power)-1):
                 
        if (power[p] > maxPower) & (frequency[p] > (1/(0.9*(len(dataArray))))):
            maxPower = power[p]
            plotPERG = frequency[p]
            PERG = 1/frequency[p]
            print(PERG)
            print(maxPower)
            
        #If PERG is larger than the number of datapoints, set the PERG to the number of datapoints
        if PERG > len(dataArray):
            PERG = len(dataArray)
    
    ax_Periodogram.axvline(x=plotPERG, color='m', ls='--', lw=1)
    ax_Periodogram.text(plotPERG + 0.01, maxPower, 'PERG Value', fontsize='x-small')
    
    ax_Periodogram.set(xlabel = 'Frequency')
    ax_Periodogram.set(ylabel = 'Lomb-Scargle Power')
        
    if dataVal == 0:
        plt.title('Average Thickness Value Periodogram')
    elif dataVal == 1:
        plt.title('Standard Deviation of Thickness Periodogram')
    elif dataVal == 2:
        plt.title('Values of Modelled Target Fraction Periodogram')
        
    ax_Periodogram.set_xlim(xmin=0)
    ax_Periodogram.set_ylim(ymin=0)

    plt.show()

    periodogramPath = path +'Periodogram.png'

    plt.savefig(periodogramPath)
        
    plt.clf()
    
    plt.close()
    
    return PERG

####################################################################################################################
##Goovaerts Plot Code##
####################################################################################################################
#Function to plot the Goovaerts graph in order to determine the GOO value
def goovaertsPlot(dataArray, dataVal, i, path):
    
    goovaerts, (ax_Goovaerts) = plt.subplots(1, sharex=True) 
    
    lag = []
    StDev = []
            
    #Determining the standard deviation at every single lag distance based upon the size of the array increasing by one each loop
    for p in range(1, len(dataArray)):
        DATA = []
        for k in range(0, p):
            DATA.append(dataArray[k])
        StDev.append(statistics.pstdev(DATA))    
        lags = p + 1
        lag.append(lags)
        
    print(lag)
    print(StDev)
    
    #Application of a 1-dimensional Gaussian filter for smoothing of the trace
    ysmoothed = gaussian_filter1d(StDev, sigma=2)
    
    ysmoothed[0] = 0
    
    data = {'Lag': lag, 'StDev': StDev, 'Smoothed': ysmoothed}
    
    df = pd.DataFrame(data)
    
    xi = np.linspace(0, len(StDev), 100)
    
    cof_u, cov = curve_fit(f, lag, StDev)

    yi = list(map(lambda x: spherical(lag, *cof_u), xi))
       
    #Best fit spherical line to the dataset
    cof, cov = curve_fit(f, lag, StDev, p0=[3., 14.])
    
    yi = list(map(lambda x: spherical(x, *cof), xi))
    
    sns.lineplot(x=xi, y=yi, color='b')
    
    print('Max: ' + str(max(yi)))
    
    count = 0 
    
    for a in yi:
        print(a)
        if a == max(yi):
            position = count
            print(position)
            break
        count += 1
        
    print(position)
    
    print(xi[position])
    
    GERG = xi[position]
    
    print(GERG)

    ###############################################################################
    
    
    ax_Goovaerts = sns.lineplot(x = 'Lag', y = 'StDev', data=df, color='r', lw=3)
    ax_Goovaerts = sns.lineplot(x = 'Lag', y = 'Smoothed', data=df, color='g', lw=3)
    
    ax_Goovaerts.legend(labels=['Spherical Variogram Fit', 'Non-Smoothed','Smoothed'], facecolor='white')
    
    ax_Goovaerts.axvline(x=GERG, color='m', ls='--', lw=1)
    ax_Goovaerts.text(GERG, 0.1, 'GOO Value', fontsize='x-small')
    
    ax_Goovaerts.set(xlabel = 'Number of Realizations')
    ax_Goovaerts.set(ylabel = 'Standard Deviation')
    
    if dataVal == 0:
        plt.title('Average Thickness Value')
    elif dataVal == 1:
        plt.title('Standard Deviation of Thickness')
    elif dataVal == 2:
        plt.title('Values of Modelled Target Fraction')
        
    ax_Goovaerts.set_xlim(xmin=0)
    ax_Goovaerts.set_ylim(ymin=0)
    
    plt.show()
            
    goovaertsPath = path +'Goovaerts.png'
            
    plt.savefig(goovaertsPath)
            
    plt.clf()
    
    plt.close()
    
    print(GERG)
    
    
    return GERG

####################################################################################################################
##Graphing Function Code##
####################################################################################################################
#Function to create the similarity plot to compare the intersectional area percentage of the full dataset and the PERG and GOO values
def similarityPlot(dataArray, dataVal, i, GERG, PERG, path): 
    
    dataArrayPERG = []
    dataArrayGERG = []
    
    comparison, (ax_Comparison) = plt.subplots(1, sharex=True) 
    
    print(PERG)
    print(GERG)
    
    for x in range(0, int(PERG)):
       dataArrayPERG.append(dataArray[x])

    
    for x in range(0, int(GERG)):
        dataArrayGERG.append(dataArray[x])
                        
    ax_Comparison = sns.distplot(x = dataArray, kde=True, color='c', hist_kws=dict(edgecolor='c', linewidth=0))
    ax_Comparison = sns.distplot(x = dataArrayGERG, kde=True, color='r', hist_kws=dict(edgecolor='r', linewidth=0))
    ax_Comparison = sns.distplot(x = dataArrayPERG, kde=True, color='g', hist_kws=dict(edgecolor='g', linewidth=0))
    ax_Comparison.legend(labels=['Full Dataset','GOO', 'PERG'], facecolor='white')
    
    ax_Comparison.set(ylabel='KDE - Kernel Density Estimation')
            
    if dataVal == 0:
        ax_Comparison.set(xlabel = 'Thickness (m)')
        plt.title('Average Thickness Histogram (m)')
    elif dataVal == 1:
        ax_Comparison.set(xlabel = 'Standard Deviation of Thickness')
        plt.title('Standard Deviation of Thickness Histogram')
    elif dataVal == 2:
        ax_Comparison.set(xlabel = 'Target Fraction (%)')
        plt.title('Values of Modelled Target Fraction Histogram (%)')
                       
    ax_Comparison.set_ylim(ymin=0)
        
    plt.show()
                           
    print("FC" + str(i) + " 95% CI")
    print(st.t.interval(alpha=0.95, df=len(dataArray)-1, loc=np.mean(dataArray), scale=st.sem(dataArray)))        
                
    comparisonPath = path +'Comparison.png'
    
    plt.savefig(comparisonPath)	
    
    plt.clf()
    
    plt.close()
    
    similarity, (ax_gerg, ax_perg) = plt.subplots(2, sharex=True)
    
    if dataVal == 0:
        ax_perg.set(xlabel='Thickness (m)')
        plt.title('Average Thickness Similarity Plot')
    elif dataVal == 1:
        ax_perg.set(xlabel='Standard Deviation of Thickness')
        plt.title('Standard Deviation of Thickness Similarity Plot')
    elif dataVal == 2:
        ax_perg.set(xlabel='Target Fraction (%)')
        plt.title('Values of Modelled Target Fraction Similarity Plot')
    
    kdeFull = gaussian_kde(dataArray, bw_method=0.3)
    
    kdePERG = gaussian_kde(dataArrayPERG, bw_method=0.3)
    
    xminP = min(min(dataArray), min(dataArrayPERG))
    xmaxP = max(max(dataArray), max(dataArrayPERG))
    
    dxP = 0.2 * (xmaxP - xminP)
    
    xminP -= dxP
    xmaxP += dxP
    
    #Plotting the PERG and Full dataset intersectional overlap percentages
    x = np.linspace(xminP, xmaxP, 500)
    kdeFull_x = kdeFull(x)
    kdePERG_x = kdePERG(x)
    intersectionPERG_x = np.minimum(kdeFull_x, kdePERG_x)
    
    ax_perg.plot(x, kdeFull_x, color='c', label='Full')
    ax_perg.fill_between(x, kdeFull_x, 0, color='c', alpha=0.2)
    ax_perg.plot(x, kdePERG_x, color='g', label='PERG')
    ax_perg.fill_between(x, kdePERG_x, 0, color='g', alpha=0.2)
    ax_perg.plot(x, intersectionPERG_x, color='m')
    ax_perg.fill_between(x, intersectionPERG_x, 0, facecolor='none', edgecolor='m', hatch='xx', label='intersection')
    
    area_intersectionPERG_x = np.trapz(intersectionPERG_x, x)
    area_intersectionPERG_x = round(area_intersectionPERG_x * 100, 1)
     
    ax_perg.legend(labels=['Full Dataset','PERG: ' + str(int(PERG)), 'Intersection: ' + str(area_intersectionPERG_x) + '%'], facecolor='white')
    plt.tight_layout()
    
    plt.ylim(ymin=0)
    
    kdeGERG = gaussian_kde(dataArrayGERG, bw_method=0.3)
    
    xminG = min(min(dataArray), min(dataArrayGERG))
    xmaxG = max(max(dataArray), max(dataArrayGERG))
    
    dxG = 0.2 * (xmaxG - xminG)
    
    xminG -= dxG
    xmaxG += dxG
    
    #Plotting the PERG and Full dataset intersectional overlap percentages
    x = np.linspace(xminG, xmaxG, 500)
    kdeFull_x = kdeFull(x)
    kdeGERG_x = kdeGERG(x)
    intersectionGERG_x = np.minimum(kdeFull_x, kdeGERG_x)
    
    ax_gerg.plot(x, kdeFull_x, color='c', label='Full')
    ax_gerg.fill_between(x, kdeFull_x, 0, color='c', alpha=0.2)
    ax_gerg.plot(x, kdeGERG_x, color='r', label='GERG')
    ax_gerg.fill_between(x, kdeGERG_x, 0, color='r', alpha=0.2)
    ax_gerg.plot(x, intersectionGERG_x, color='m')
    ax_gerg.fill_between(x, intersectionGERG_x, 0, facecolor='none', edgecolor='m', hatch='xx', label='intersection')
    
    area_intersectionGERG_x = np.trapz(intersectionGERG_x, x)
    area_intersectionGERG_x = round(area_intersectionGERG_x * 100, 1)
     
    ax_gerg.legend(labels=['Full Dataset','GOO: ' + str(int(GERG)), 'Intersection: ' + str(area_intersectionGERG_x) + '%'], facecolor='white')
    
    ax_perg.set(ylabel='KDE - Kernel Density Estimation')
    
    print(area_intersectionGERG_x)
    print(area_intersectionPERG_x)
        
    similarityPath = path +'Similarity.png'

    ax_perg.set_ylim([0, None])
    ax_gerg.set_ylim([0, None])
    
    plt.tight_layout()
    
    plt.show()
    
    plt.savefig(similarityPath)	
    
    plt.clf()
       
    #################################################################################################
    
    
    potentialVals = []

    high_area_intersectionPotential_x = 0    
    
    print('GOO: ' + str(GERG))
    print('PERG: ' + str(PERG))
    
    #Incrementally increasing the values between the GERG and the PERG and comparing this to the Full dataset intersectional overlap percentages
    if GERG < PERG:
        for x in range(int(GERG), int(PERG)+1):
            potentialVals.append(x)
        
        print('Potential Vals: ' + str(potentialVals))
        
        for pV in potentialVals:
            potentialArray = []
            for r in range(0, int(pV)):
                potentialArray.append(dataArray[r])
            
            print('Potential Array: ' + str(potentialArray))
            
            kdeFull = gaussian_kde(dataArray, bw_method=0.3)
            kdePotential = gaussian_kde(potentialArray, bw_method=0.3)
            
            xminP = min(min(dataArray), min(potentialArray))
            xmaxP = max(max(dataArray), max(potentialArray))
            
            dxP = 0.2 * (xmaxP - xminP)
            
            xminP -= dxP
            xmaxP += dxP
           
            x = np.linspace(xminP, xmaxP, 500)
            kdeFull_x = kdeFull(x)
            kdePotential_x = kdePotential(x)
            intersectionPotential_x = np.minimum(kdeFull_x, kdePotential_x)
            
            area_intersectionPotential_x = np.trapz(intersectionPotential_x, x)
            area_intersectionPotential_x = round(area_intersectionPotential_x * 100, 1)
            
            print(area_intersectionPotential_x)
            
            if area_intersectionPotential_x > high_area_intersectionPotential_x:
                high_area_intersectionPotential_x = area_intersectionPotential_x
                RR = pV
                print(RR)
                print('HERE')
    else:
        RR = int(PERG)+1
   
    similarity, (ax_RR) = plt.subplots(1, sharex=True) 
   
    dataArrayRR = []

    for r in range(0, int(RR)):
        dataArrayRR.append(dataArray[r])    
    
    kdeRR = gaussian_kde(dataArrayRR, bw_method=0.3)
    
    xminG = min(min(dataArray), min(dataArrayRR))
    xmaxG = max(max(dataArray), max(dataArrayRR))
    
    dxG = 0.2 * (xmaxG - xminG)
    
    xminG -= dxG
    xmaxG += dxG
    
    x = np.linspace(xminG, xmaxG, 500)
    kdeFull_x = kdeFull(x)
    kdeRR_x = kdeRR(x)
    intersectionRR_x = np.minimum(kdeFull_x, kdeRR_x)
    
    ax_RR.plot(x, kdeFull_x, color='c', label='Full')
    ax_RR.fill_between(x, kdeFull_x, 0, color='c', alpha=0.2)
    ax_RR.plot(x, kdeRR_x, color='r', label='MK')
    ax_RR.fill_between(x, kdeRR_x, 0, color='r', alpha=0.2)
    ax_RR.plot(x, intersectionRR_x, color='m')
    ax_RR.fill_between(x, intersectionRR_x, 0, facecolor='none', edgecolor='m', hatch='xx', label='intersection')
    
    area_intersectionRR_x = np.trapz(intersectionRR_x, x)
    area_intersectionRR_x = round(area_intersectionRR_x * 100, 1)
     
    ax_RR.legend(labels=['Full Dataset','RR: ' + str(RR), 'Intersection: ' + str(area_intersectionRR_x) + '%'], facecolor='white')
    
    ax_RR.set(ylabel='KDE - Kernel Density Estimation')
    
    if dataVal == 0:
        ax_RR.set(xlabel='Thickness (m)')
        plt.title('Average Thickness Similarity Recommended Runs')
    elif dataVal == 1:
        ax_RR.set(xlabel='Standard Deviation of Thickness')
        plt.title('Standard Deviation of Thickness Similarity Recommended Runs')
    elif dataVal == 2:
        ax_RR.set(xlabel='Target Fraction (%)')
        plt.title('Values of Modelled Target Fraction Similarity Recommended Runs')
        
    similarityPathRR = path +'SimilarityRR.png'
    
    plt.ylim(ymin=0)
    
    plt.tight_layout()
    
    plt.show()  
    
    plt.savefig(similarityPathRR)	
    
    plt.clf()
    
    ##########################################################################################################

    dataArrayCurrent = []

    similarity, (ax_Current) = plt.subplots(1, sharex=True) 

    for x in range(0, 20):
        dataArrayCurrent.append(dataArray[x])

    if dataVal == 0:
        ax_Current.set(xlabel='Thickness (m)')
        plt.title('Average Thickness Similarity to 20 Runs')
    elif dataVal == 1:
        ax_Current.set(xlabel='Standard Deviation of Thickness')
        plt.title('Standard Deviation of Thickness Similarity to 20 Runs')
    elif dataVal == 2:
        ax_Current.set(xlabel='Target Fraction (%)')
        plt.title('Values of Modelled Target Fraction Similarity to 20 Runs')
    
    ax_Current.set(ylabel='KDE - Kernel Density Estimation')
    
    kdeFull = gaussian_kde(dataArray, bw_method=0.3)
    
    kdeCurrent = gaussian_kde(dataArrayCurrent, bw_method=0.3)
       
    xminP = min(min(dataArray), min(dataArrayCurrent))
    xmaxP = max(max(dataArray), max(dataArrayCurrent))
    
    dxP = 0.2 * (xmaxP - xminP)
    
    xminP -= dxP
    xmaxP += dxP
    
    #Plotting the first 20 and Full dataset intersectional overlap percentages
    x = np.linspace(xminP, xmaxP, 500)
    kdeFull_x = kdeFull(x)
    kdeCurrent_x = kdeCurrent(x)
    intersectionCurrent_x = np.minimum(kdeFull_x, kdeCurrent_x)
    
    ax_Current.plot(x, kdeFull_x, color='c', label='Full')
    ax_Current.fill_between(x, kdeFull_x, 0, color='c', alpha=0.2)
    ax_Current.plot(x, kdeCurrent_x, color='g', label='PERG')
    ax_Current.fill_between(x, kdeCurrent_x, 0, color='g', alpha=0.2)
    ax_Current.plot(x, intersectionCurrent_x, color='m')
    ax_Current.fill_between(x, intersectionCurrent_x, 0, facecolor='none', edgecolor='m', hatch='xx', label='intersection')
    
    area_intersectionCurrent_x = np.trapz(intersectionCurrent_x, x)
    area_intersectionCurrent_x = round(area_intersectionCurrent_x * 100, 1)
     
    ax_Current.legend(labels=['Full Dataset','Current Suggested Runs: 20', 'Intersection: ' + str(area_intersectionCurrent_x) + '%'], facecolor='white')
    
    ax_Current.set_ylim(ymin=0)
    
    plt.tight_layout()
    
    plt.show()  
    
    CurrentPath = path +'CurrentSuggestedRuns.png'
    
    plt.savefig(CurrentPath)	
    
    plt.clf()
    
    return RR, area_intersectionRR_x, area_intersectionCurrent_x

####################################################################################################################
##Compiled Final Values Code##
####################################################################################################################
#Function to output the GOO, PERG, RR and Similarity percentages to a new CSV file
def OutputValues(outputDF):
    
    fileLoc = destinationFolder + '/Output.csv'

    #Creates and opens the compiled file in the destination folder
    with open(fileLoc, mode='w', newline='') as newTestFile:        
        
        #Appending the freshly cleaned data into the newly created file
        outputDF.to_csv(fileLoc, index = False, columns=['File Name', 'Mean PERG', 'Mean GOO', 'Mean RR', 'Mean Similarity', 'Mean GOO Similarity', 'StDEV PERG', 'StDEV GOO', 'StDEV RR', 'StDEV Similarity', 'StDev GOO Similarity', 'TF PERG', 'TF GOO', 'TF RR', 'TF Similarity', 'TF GOO Similarity'])

####################################################################################################################
##Graphing Function Code##
####################################################################################################################
#Function to decide which process functions are run
def graphingFunction():
    
    sns.set(style="whitegrid")
    
    outputColumns = {'File Name', 'Mean PERG', 'Mean GOO', 'Mean RR', 'Mean Similarity', 'Mean GOO Similarity', 'StDEV PERG', 'StDEV GOO', 'StDEV RR', 'StDEV Similarity', 'StDev GOO Similarity', 'TF PERG', 'TF GOO', 'TF RR', 'TF Similarity', 'TF GOO Similarity'}
    
    indexVals = []
    
    for x in range(0, len(identifierArray)):
        indexVals.append(x)

    outputDF = pd.DataFrame(columns = outputColumns, index=indexVals)
    
    print(outputDF)
    
    indexVal = 0
    
    for i in identifierArray:
        PERG = 0
        GERG = 0
        dataArrayMean = []
        dataArraySTDEV = []
        dataArrayTF = []
        
        outputDF.loc[indexVal]['File Name'] = i
               
        for j in range(len(globals()['data%s' % i])):
            dataArrayMean.append(globals()['data%s' % i][j]['Mean'])
            dataArraySTDEV.append(globals()['data%s' % i][j]['Standard Deviation'])
            dataArrayTF.append(globals()['data%s' % i][j]['Percentage'])
            
        dataVal = 0
        while dataVal < 3:
            dataArray = []
            
            ## INSERT FIGURE CREATION HERE
            
            if dataVal == 0:
                dataArray = dataArrayMean
                path = destinationFolder + '/Images/' + i + '/' + i + '_Mean_'
            elif dataVal == 1:
                dataArray = dataArraySTDEV
                path = destinationFolder + '/Images/' + i + '/' + i + '_STDEV_'
            elif dataVal == 2:
                dataArray = dataArrayTF
                path = destinationFolder + '/Images/' + i + '/' + i + '_TF_'

            ##### RENAME THESE// NEEDED?????
            histogramBoxPlot(dataArray, dataVal, i, path)
            PERG = periodogramPlot(dataArray, dataVal, i, path)
            GERG = goovaertsPlot(dataArray, dataVal, i, path)
                
            RR, percentageSim, percentageSimGOO = similarityPlot(dataArray, dataVal, i, GERG, PERG, path)    
            
            if dataVal == 0:
                outputDF.iloc[indexVal]['Mean PERG'] = PERG
                outputDF.iloc[indexVal]['Mean GOO'] = GERG
                outputDF.iloc[indexVal]['Mean RR'] = RR
                outputDF.iloc[indexVal]['Mean Similarity'] = percentageSim
                outputDF.iloc[indexVal]['Mean GOO Similarity'] = percentageSimGOO
            elif dataVal == 1:
                outputDF.iloc[indexVal]['StDEV PERG'] = PERG
                outputDF.iloc[indexVal]['StDEV GOO'] = GERG
                outputDF.iloc[indexVal]['StDEV RR'] = RR 
                outputDF.iloc[indexVal]['StDEV Similarity'] = percentageSim
                outputDF.iloc[indexVal]['StDev GOO Similarity'] = percentageSimGOO
            elif dataVal == 2:
                outputDF.iloc[indexVal]['TF PERG'] = PERG
                outputDF.iloc[indexVal]['TF GOO'] = GERG
                outputDF.iloc[indexVal]['TF RR'] = RR
                outputDF.iloc[indexVal]['TF Similarity'] = percentageSim
                outputDF.iloc[indexVal]['TF GOO Similarity'] = percentageSimGOO
                
            
            dataVal = dataVal + 1
            
            plt.close('all')
            
        indexVal = indexVal + 1
        
        plt.close('all')
    
    OutputValues(outputDF)
    
    plt.close('all')

#Creates a spherical fit to the input data for the GOO value to be determined
def f(h, a, b):
    return spherical(h, a, b)
    
####################################################################################################################
##Get Target Folder Code##
####################################################################################################################
#Globally sets the location of the target folder for the various input files to be found within
def getTargetFolder():   
    global targetFolder 
    targetFolder = filedialog.askdirectory(parent=root, title=ttl)
    messagebox.showinfo("The target folder selected is: ", targetFolder)
    v5.set(1)
    targetFiles.append(os.listdir(targetFolder))
    print(targetFiles)
    
####################################################################################################################
##Get Destination Folder Code##
####################################################################################################################
#Globally sets the location of the destination folder for the various files and outputs to be assigned to
def getDestinationFolder():
    global destinationFolder
    destinationFolder = filedialog.askdirectory(parent=root, title=ttl)
    messagebox.showinfo("The destination folder selected is: ", destinationFolder)
    v6.set(1)

####################################################################################################################
##UI Design##
####################################################################################################################

#Create the buttons, labels and tickboxes for the UI
#Assigns functionality to the various buttons
btnClear=Button(window, text="Clear", fg='black', command=clearForm)
btnClear.place(x=50, y=250)

btnStart=Button(window, text="Start", fg='black', command=startProcess)
btnStart.place(x=150, y=250)

btnTarget=Button(window, text="Select Target Folder", fg='black', command=getTargetFolder)
btnTarget.place(x=50, y=110)

btnDestination=Button(window, text="Select Destination Folder", fg='black', command=getDestinationFolder)
btnDestination.place(x=50, y=160)

label = ttk.Label(text="Algorithm Type:")
label.place(x=250, y=250)

algorithmCombo=ttk.Combobox(window)
algorithmCombo.place(x=350, y=250)
algorithmCombo['values']=('SIS', 'OBM', 'MPS')
algorithmCombo['state'] = 'readonly'

lblTitle=Label(window, text="StReAMS", fg='black', font=("Helvetica", 20))
lblTitle.place(x=50, y=15)


v5 = IntVar()
v6 = IntVar()

C5 = Checkbutton(window, variable = v5)
C6 = Checkbutton(window, variable = v6)

#Placement of the tickboxes

C5.place(x=225, y=110)
C6.place(x=225, y=160)

window.mainloop()