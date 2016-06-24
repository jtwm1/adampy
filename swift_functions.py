from __future__ import print_function
import numpy as np
import sys
import requests
import os
from scipy.optimize import curve_fit
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt
import emcee

###################################################
##### Download All Swift GRB XRT Light Curves #####
###################################################

def all_lc_data(directory):
    """ Downloads and saves all Swift detected GRB XRT afterglow data to
    directory of your choice """
    """ This first part gets the GRB name and trigger list """
    r = requests.get('http://swift.gsfc.nasa.gov/archive/grb_table/table.php?obs=Swift&year=All+Years&restrict=none&grb_trigger=1').text
    soup = BeautifulSoup(r,'lxml')
    table = soup.find('table')
    rows = table.findAll('tr')
    grb_names = []
    trigger_nos = []
    c = 0

    for row in rows:
    
        if c >= 1:
            cols = row.findAll('td')
            cols = [ele.text.strip() for ele in cols]

            if len(cols) == 2:
                grb_name = cols[0]            
                trigger_no = cols[1]

                if len(trigger_no) == 6:
                    grb_names.append(grb_name)
                    trigger_nos.append(trigger_no)
                    
                if len(trigger_no) == 26:
                    grb_names.append(grb_name)
                    trigger_no = str(trigger_no)[21:26]
                    trigger_nos.append(trigger_no)

                if len(trigger_no) == 21:
                    grb_names.append(grb_name)
                    
                    if str(trigger_no)[0] == 'G':
                        trigger_no = str(trigger_no)[15:21]
                        trigger_nos.append(trigger_no)

                    else:
                        trigger_no = str(trigger_no)[0:6]
                        trigger_nos.append(trigger_no)

                if len(trigger_no) == 12:
                    grb_names.append(grb_name)
                    trigger_no = str(trigger_no)[6:12]
                    trigger_nos.append(trigger_no)
                
                if len(trigger_no) == 20:
                    grb_names.append(grb_name)
                    trigger_no = str(trigger_no)[15:20]
                    trigger_nos.append(trigger_no)

                if len(trigger_no) == 32:
                    grb_names.append(grb_name)
                    trigger_no = str(trigger_no)[27:32]
                    trigger_nos.append(trigger_no)

                if trigger_no == 'BATSS':
                    grb_names.append(grb_name)
                    trigger_nos.append(trigger_no)

                if trigger_no == 'Ground Analysis':
                    grb_names.append(grb_name)
                    trigger_nos.append(trigger_no)
        c += 1
    
    """ Downloads all available data if it isn't already downloaded """  
    for i in range(0,len(trigger_nos),1):
        lines = []
        d = 0
    
        if len(trigger_nos[i]) >= 6:
            xray_lc_data_url = 'http://www.swift.ac.uk/xrt_curves/00' + trigger_nos[i] + '/curve.qdp'

        if len(trigger_nos[i]) == 5:
            xray_lc_data_url = 'http://www.swift.ac.uk/xrt_curves/000' + trigger_nos[i] + '/curve.qdp'

        xray_lc_file = directory + 'GRB' + grb_names[i] + '.qdp'

        if os.path.exists(xray_lc_file):
            print('Already have XRT afterglow data for GRB', grb_names[i])
        
        else:  
            r2 = requests.get(xray_lc_data_url).text
            for line in r2.split('\n'):
                lines.append(line)
                d += 1
            
            if lines[0] == '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">':
                print('Cannot find XRT afterglow Data for GRB', grb_names[i])

            else:
                orig_stdout = sys.stdout 
                xray_lc_out = open(xray_lc_file, "w")
                sys.stdout = xray_lc_out
                for j in range(0,len(lines),1):
                    if lines[j] != '':
                        print(lines[j])

                sys.stdout = orig_stdout
                xray_lc_out.close()
                print('Retrieved afterglow data for GRB', grb_names[i])

#########################################
##### X-Ray Afterglow File Plotting #####
#########################################

def xray_afterglow_plot(afterglow_file):
    """ Extracts data from QDP file and plots the x-ray afterglow. Also makes
    the data useful to analyse later """
    f = open(afterglow_file,'r')
   
    """ Create plotting environment for x-ray afterglow """
    fig,ax = plt.subplots()
    plt.title('X-Ray Afterglow')
    plt.xlabel('Time since GRB trigger (s)',fontsize=13)
    plt.ylabel('Count Rate (0.3-10 keV) (s$^{-1}$)',fontsize=13)
    plt.yscale('log')
    plt.xscale('log')

    """ Counter/Array Set Up """
    c = 0
    wtslew_start = 0
    wtslew_end = 0
    wt_start = 0
    wt_end = 0
    pc_start = 0
    pc_end = 0
    uplim_start = 0
    uplim_end = 0
    endoffile = 0
    
    data = []
    wtslew_time = []
    wtslew_timepos = []
    wtslew_timeneg = []
    wtslew_rate = []
    wtslew_ratepos = []
    wtslew_rateneg = []
    wt_time = []
    wt_timepos = []
    wt_timeneg = []
    wt_rate = []
    wt_ratepos = []
    wt_rateneg = []
    pc_time = []
    pc_timepos = []
    pc_timeneg = []
    pc_rate = []
    pc_ratepos = []
    pc_rateneg = []
    uplim_time = []
    uplim_timepos = []
    uplim_timeneg = []
    uplim_rate = []
    uplim_ratepos = []
    uplim_rateneg = []

    """ Read the file and plot the contents """
    for line in f:
        temp = line.strip()
        data.append(temp)
        c += 1        
       
        if line.strip() == '! WTSLEW data':
            wtslew_start = c + 1     
      
        if line.strip() == '! WT data':
            wt_start = c + 1
            wtslew_end = wt_start - 3
        
        if line.strip() == '! PC data':
            pc_start = c + 1
            wt_end = pc_start - 3
            
            if wt_start == 0:
               wtslew_end = pc_start - 3 
        
        if line.strip() == '! PC Upper limit':
            uplim_start = c
            pc_end = uplim_start - 2
            uplim_end = uplim_start + 1
            
            if pc_start == 0:
                wt_end = uplim_start - 3
                
    if f.readline() == '':
        endoffile = c
            
    if wtslew_end == 0:
        wtslew_end == endoffile
        
    """ WT Slew Data """    
    if wtslew_start != 0:   
        for p in range (wtslew_end - wtslew_start):
        
            wtslew_time += [float(data[wtslew_start + p].split()[0])]
            wtslew_timepos +=[float(data[wtslew_start + p].split()[1])]           
            wtslew_timeneg += [-1*float(data[wtslew_start + p].split()[2])]
            wtslew_rate += [float(data[wtslew_start + p].split()[3])]
            wtslew_ratepos += [float(data[wtslew_start + p].split()[4])]
            wtslew_rateneg += [-1*float(data[wtslew_start + p].split()[5])]
                        
    if wtslew_start != 0:
        
       wtslew_timeerr = [wtslew_timeneg,wtslew_timepos]
       wtslew_rateerr = [wtslew_rateneg,wtslew_ratepos]
       ax.scatter(wtslew_time,wtslew_rate,color='skyblue',marker='.',zorder=1)
       ax.errorbar(wtslew_time,wtslew_rate,xerr=wtslew_timeerr,yerr=wtslew_rateerr,capsize=0,fmt='o',color='skyblue',marker='.',zorder=1)
            
    if wt_end == 0:
        wt_end = endoffile        

    """ WT Data """   
    if wt_start != 0:
        for q in range (wt_end - wt_start):
        
            wt_time += [float(data[wt_start + q].split()[0])]
            wt_timepos += [float(data[wt_start + q].split()[1])]           
            wt_timeneg += [-1*float(data[wt_start + q].split()[2])]
            wt_rate += [float(data[wt_start + q].split()[3])]
            wt_ratepos += [float(data[wt_start + q].split()[4])]
            wt_rateneg += [-1*float(data[wt_start + q].split()[5])]
            
    if wt_start != 0:
        
        wt_timeerr = [wt_timeneg,wt_timepos]
        wt_rateerr = [wt_rateneg,wt_ratepos]
        ax.scatter(wt_time,wt_rate,color='blue',marker='.',zorder=1)
        ax.errorbar(wt_time,wt_rate,xerr=wt_timeerr,yerr=wt_rateerr,capsize=0,fmt='o',color='blue',marker='.',zorder=1)
          
    if pc_end == 0:
        pc_end = endoffile        

    """ PC Data """
    if pc_start != 0:    
        for r in range (pc_end - pc_start):
        
            pc_time += [float(data[pc_start + r].split()[0])]
            pc_timepos += [float(data[pc_start + r].split()[1])]           
            pc_timeneg += [-1*float(data[pc_start + r].split()[2])]
            pc_rate += [float(data[pc_start + r].split()[3])]
            pc_ratepos += [float(data[pc_start + r].split()[4])]
            pc_rateneg += [-1*float(data[pc_start + r].split()[5])]
    
    if pc_start != 0:
        
        pc_timeerr = [pc_timeneg,pc_timepos]
        pc_rateerr = [pc_rateneg,pc_ratepos]
        ax.scatter(pc_time,pc_rate,color='red',marker='.',zorder=1)
        ax.errorbar(pc_time,pc_rate,xerr=pc_timeerr,yerr=pc_rateerr,capsize=0,fmt='o',color='red',marker='.',zorder=1)             
    
    if uplim_end == 0:
        uplim_end == endoffile        

    """ Upper Limit Data """    
    if uplim_start != 0:
        for s in range (uplim_end - uplim_start):
        
            uplim_time += [float(data[uplim_start + s].split()[0])]
            uplim_timepos += [float(data[uplim_start + s].split()[1])]           
            uplim_timeneg += [-1*float(data[uplim_start + s].split()[2])]
            uplim_rate += [float(data[uplim_start + s].split()[3])]
            uplim_ratepos += [float(data[uplim_start + s].split()[4])]
            uplim_rateneg += [float(data[uplim_start + s].split()[5])]
     
    if uplim_start != 0:
        
        uplim_timeerr = [uplim_timeneg,uplim_timepos]
        uplim_rateerr = [uplim_rateneg,uplim_ratepos]
        symbols = [u'\u2193']      
        for i, symbol in enumerate(symbols):
                ax.errorbar(uplim_time,uplim_rate,xerr=uplim_timeerr,yerr=uplim_rateerr,capsize=0,fmt='o',color='red',marker='',zorder=1)     

                for uplim_time, uplim_rate in zip(uplim_time,uplim_rate):
                    ax.text(uplim_time,uplim_rate,symbol,color='red',fontname='STIXGeneral',size=20,va='top',ha='center',clip_on=True,zorder=1)
                    
    f.close()
    plt.draw()

    """ Create arrays to make use of data later """
    time = wtslew_time + wt_time + pc_time
    time = np.asarray(time)
    timepos = wtslew_timepos + wt_timepos + pc_timepos
    timepos = np.asarray(timepos)
    timeneg = wtslew_timeneg + wt_timeneg + pc_timeneg
    timeneg = np.asarray(timeneg)
    rate = wtslew_rate + wt_rate + pc_rate
    rate = np.asarray(rate)
    ratepos = wtslew_ratepos + wt_ratepos + pc_ratepos
    ratepos = np.asarray(ratepos)
    rateneg = wtslew_rateneg + wt_rateneg + pc_rateneg
    rateneg = np.asarray(rateneg)

    return(time,timepos,timeneg,rate,ratepos,rateneg)

#####################################################
##### Download All Swift GRB XRT Hardness Files ##### 
#####################################################

def all_hr_data(directory):
    """ Downloads and saves all Swift detected GRB XRT hardness data to
    directory of your choice """
    """ This first part gets the GRB name and trigger list """
    r = requests.get('http://swift.gsfc.nasa.gov/archive/grb_table/table.php?obs=Swift&year=All+Years&restrict=none&grb_trigger=1').text
    soup = BeautifulSoup(r,'lxml')
    table = soup.find('table')
    rows = table.findAll('tr')
    grb_names = []
    trigger_nos = []
    c = 0

    for row in rows:
    
        if c >= 1:
            cols = row.findAll('td')
            cols = [ele.text.strip() for ele in cols]

            if len(cols) == 2:
                grb_name = cols[0]            
                trigger_no = cols[1]
                
                if len(trigger_no) == 6:
                    grb_names.append(grb_name)
                    trigger_nos.append(trigger_no)
                    
                if len(trigger_no) == 26:
                    grb_names.append(grb_name)
                    trigger_no = str(trigger_no)[21:26]
                    trigger_nos.append(trigger_no)

                if len(trigger_no) == 21:
                    grb_names.append(grb_name)
                    
                    if str(trigger_no)[0] == 'G':
                        trigger_no = str(trigger_no)[15:21]
                        trigger_nos.append(trigger_no)

                    else:
                        trigger_no = str(trigger_no)[0:6]
                        trigger_nos.append(trigger_no)

                if len(trigger_no) == 12:
                    grb_names.append(grb_name)
                    trigger_no = str(trigger_no)[6:12]
                    trigger_nos.append(trigger_no)
                
                if len(trigger_no) == 20:
                    grb_names.append(grb_name)
                    trigger_no = str(trigger_no)[15:20]
                    trigger_nos.append(trigger_no)

                if len(trigger_no) == 32:
                    grb_names.append(grb_name)
                    trigger_no = str(trigger_no)[27:32]
                    trigger_nos.append(trigger_no)

                if trigger_no == 'BATSS':
                    grb_names.append(grb_name)
                    trigger_nos.append(trigger_no)

                if trigger_no == 'Ground Analysis':
                    grb_names.append(grb_name)
                    trigger_nos.append(trigger_no)
        c += 1
        
    """ Downloads all available data if it isn't already downloaded """ 
    for i in range(0,len(trigger_nos),1):
        lines = []
        d = 0
    
        if len(trigger_nos[i]) >= 6:
            xray_hr_data_url = 'http://www.swift.ac.uk/xrt_curves/00' + trigger_nos[i] + '/hardrat.qdp'

        if len(trigger_nos[i]) == 5:
            xray_hr_data_url = 'http://www.swift.ac.uk/xrt_curves/000' + trigger_nos[i] + '/hardrat.qdp'

        xray_hr_file = directory + 'GRB' + grb_names[i] + '.qdp'

        if os.path.exists(xray_hr_file):
            print('Already have XRT Hardness data for GRB', grb_names[i])
        
        else:  
            r = requests.get(xray_hr_data_url).text
            for line in r.split('\n'):
                lines.append(line)
                d += 1
            
            if lines[0] == '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">':
                print('Cannot find XRT Hardness Data for GRB', grb_names[i])

            else:
                orig_stdout = sys.stdout 
                xray_hr_out = open(xray_hr_file, "w")
                sys.stdout = xray_hr_out
                for j in range(0,len(lines),1):
                    if lines[j] != '':
                        print(lines[j])

                sys.stdout = orig_stdout
                xray_hr_out.close()
                print('Retrieved hardness data for GRB', grb_names[i])

########################################
##### X-ray Hardness File Plotting #####
########################################

def xray_hardness_plot(hardness_file):
    """ Opens and plots the hardness ratio qdp file and makes arrays for
    further analysis of the values """
    f=open(hardness_file,'r')

    """ Create plotting environment for hard x-ray, soft x-ray and
    ratio data """
    fig = plt.figure()
    
    ax1 = fig.add_subplot(311)
    ax1.set_title('Hardness Ratio')
    ax1.set_ylabel('(1.51-10 keV) c s$^{-1}$')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    ax2 = fig.add_subplot(312,sharex=ax1)
    ax2.set_ylabel('(0.3-1.51) keV c s$^{-1}$')
    ax2.set_yscale('log')

    ax3 = fig.add_subplot(310,sharex=ax1)
    ax3.set_ylabel('Ratio')
    ax3.set_xlabel('Time since GRB trigger (s)')
    
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    """ Counter/Array Set Up """

    c = 0
    endoffile = 0
    wtsthard_start = 0
    wtsthard_stop = 0
    wtstsoft_start = 0
    wtstsoft_stop = 0
    wtstratio_start = 0
    wtstratio_stop = 0
    wthard_start = 0
    wthard_stop = 0
    wtsoft_start = 0
    wtsoft_stop = 0
    wtratio_start = 0
    wtratio_stop = 0
    pchard_start = 0
    pchard_stop = 0
    pcsoft_start = 0
    pcsoft_stop = 0
    pcratio_start = 0
    pcratio_stop = 0

    data = []
    wtsthard_time = []
    wtsthard_timepos = []
    wtsthard_timeneg = []
    wtsthard_rate = []
    wtsthard_ratepos = []
    wtsthard_rateneg = []
    wtstsoft_time = []
    wtstsoft_timepos = []
    wtstsoft_timeneg = []
    wtstsoft_rate = []
    wtstsoft_ratepos = []
    wtstsoft_rateneg = []
    wtstratio_time = []
    wtstratio_timepos = []
    wtstratio_timeneg = []
    wtstratio_rate = []
    wtstratio_ratepos = []
    wtstratio_rateneg = []
    wthard_time = []
    wthard_timepos = []
    wthard_timeneg = []
    wthard_rate = []
    wthard_ratepos = []
    wthard_rateneg = []
    wtsoft_time = []
    wtsoft_timepos = []
    wtsoft_timeneg = []
    wtsoft_rate = []
    wtsoft_ratepos = []
    wtsoft_rateneg = []
    wtratio_time = []
    wtratio_timepos = []
    wtratio_timeneg = []
    wtratio_rate = []
    wtratio_ratepos = []
    wtratio_rateneg = []
    pchard_time = []
    pchard_timepos = []
    pchard_timeneg = []
    pchard_rate = []
    pchard_ratepos = []
    pchard_rateneg = []
    pcsoft_time = []
    pcsoft_timepos = []
    pcsoft_timeneg = []
    pcsoft_rate = []
    pcsoft_ratepos = []
    pcsoft_rateneg = []
    pcratio_time = []
    pcratio_timepos = []
    pcratio_timeneg = []
    pcratio_rate = []
    pcratio_ratepos = []
    pcratio_rateneg = []

    """ Read the file and plot the contents """
    for line in f:
        temp = line.strip()
        data.append(temp)
        c += 1

        if line.strip() == '! WTST -- hard data':
            wtsthard_start = c + 1

        if line.strip() == '! WTST -- soft data':
            wtstsoft_start = c + 1
            wtsthard_stop = wtstsoft_start - 3

        if line.strip() == '! WTST -- hardness ratio':
            wtstratio_start = c + 1
            wtstsoft_stop = wtstratio_start - 3

        if line.strip() == '! WT -- hard data':
            wthard_start = c + 1
            wtstratio_stop = wthard_start - 3

        if line.strip() == '! WT -- soft data':
            wtsoft_start = c + 1
            wthard_stop = wtsoft_start - 3

        if line.strip() == '! WT -- hardness ratio':
            wtratio_start = c + 1
            wtsoft_stop = wtratio_start - 3

        if line.strip() == '! PC -- hard data':
            pchard_start = c + 1
            wtratio_stop = pchard_start - 3

            if wthard_start == 0:
                wtstratio_stop = pchard_start - 3

        if line.strip() == '! PC -- soft data':
            pcsoft_start = c + 1
            pchard_stop = pcsoft_start - 3

        if line.strip() == '! PC -- hardness ratio':
            pcratio_start = c + 1
            pcsoft_stop = pcratio_start - 3

    if f.readline() == '':
            endoffile = c

    if wtstratio_stop == 0:
        wtstratio_stop = endoffile

    """ WT Slew Data """
    if wtsthard_start != 0:   
        for p in range (wtsthard_stop - wtsthard_start):
        
            wtsthard_time += [float(data[wtsthard_start + p].split()[0])]
            wtsthard_timepos +=[float(data[wtsthard_start + p].split()[1])]           
            wtsthard_timeneg += [-1*float(data[wtsthard_start + p].split()[2])]
            wtsthard_rate += [float(data[wtsthard_start + p].split()[3])]
            wtsthard_ratepos += [float(data[wtsthard_start + p].split()[4])]
            wtsthard_rateneg += [float(data[wtsthard_start + p].split()[4])]

    if wtsthard_start != 0:

       wtsthard_timeerr = [wtsthard_timeneg,wtsthard_timepos]
       wtsthard_rateerr = [wtsthard_rateneg,wtsthard_ratepos]
       ax1.scatter(wtsthard_time,wtsthard_rate,color='skyblue',marker='.',zorder=1)
       ax1.errorbar(wtsthard_time,wtsthard_rate,xerr=wtsthard_timeerr,yerr=wtsthard_rateerr,capsize=0,fmt='o',color='skyblue',marker='.',zorder=1)

    if wtstsoft_start != 0:   
        for q in range (wtstsoft_stop - wtstsoft_start):
        
            wtstsoft_time += [float(data[wtstsoft_start + q].split()[0])]
            wtstsoft_timepos +=[float(data[wtstsoft_start + q].split()[1])]           
            wtstsoft_timeneg += [-1*float(data[wtstsoft_start + q].split()[2])]
            wtstsoft_rate += [float(data[wtstsoft_start + q].split()[3])]
            wtstsoft_ratepos += [float(data[wtstsoft_start + q].split()[4])]
            wtstsoft_rateneg += [float(data[wtstsoft_start + q].split()[4])]

    if wtstsoft_start != 0:
        
       wtstsoft_timeerr = [wtstsoft_timeneg,wtstsoft_timepos]
       wtstsoft_rateerr = [wtstsoft_rateneg,wtstsoft_ratepos]
       ax2.scatter(wtstsoft_time,wtstsoft_rate,color='skyblue',marker='.',zorder=1)
       ax2.errorbar(wtstsoft_time,wtstsoft_rate,xerr=wtstsoft_timeerr,yerr=wtstsoft_rateerr,capsize=0,fmt='o',color='skyblue',marker='.',zorder=1)

    if wtstratio_start != 0:   
        for r in range (wtstratio_stop - wtstratio_start):
        
            wtstratio_time += [float(data[wtstratio_start + r].split()[0])]
            wtstratio_timepos +=[float(data[wtstratio_start + r].split()[1])]           
            wtstratio_timeneg += [-1*float(data[wtstratio_start + r].split()[2])]
            wtstratio_rate += [float(data[wtstratio_start + r].split()[3])]
            wtstratio_ratepos += [float(data[wtstratio_start + r].split()[4])]
            wtstratio_rateneg += [float(data[wtstratio_start + r].split()[4])]

    if wtstratio_start != 0:
        
       wtstratio_timeerr = [wtstratio_timeneg,wtstratio_timepos]
       wtstratio_rateerr = [wtstratio_rateneg,wtstratio_ratepos]
       ax3.scatter(wtstratio_time,wtstratio_rate,color='skyblue',marker='.',zorder=1)
       ax3.errorbar(wtstratio_time,wtstratio_rate,xerr=wtstratio_timeerr,yerr=wtstratio_rateerr,capsize=0,fmt='o',color='skyblue',marker='.',zorder=1)

    if wtratio_stop == 0:
        wtratio_stop = endoffile

    """ WT Data """
    if wthard_start != 0:   
        for p in range (wthard_stop - wthard_start):
        
            wthard_time += [float(data[wthard_start + p].split()[0])]
            wthard_timepos +=[float(data[wthard_start + p].split()[1])]           
            wthard_timeneg += [-1*float(data[wthard_start + p].split()[2])]
            wthard_rate += [float(data[wthard_start + p].split()[3])]
            wthard_ratepos += [float(data[wthard_start + p].split()[4])]
            wthard_rateneg += [float(data[wthard_start + p].split()[4])]

    if wthard_start != 0:

       wthard_timeerr = [wthard_timeneg,wthard_timepos]
       wthard_rateerr = [wthard_rateneg,wthard_ratepos]
       ax1.scatter(wthard_time,wthard_rate,color='blue',marker='.',zorder=1)
       ax1.errorbar(wthard_time,wthard_rate,xerr=wthard_timeerr,yerr=wthard_rateerr,capsize=0,fmt='o',color='blue',marker='.',zorder=1)

    if wtsoft_start != 0:   
        for q in range (wtsoft_stop - wtsoft_start):
        
            wtsoft_time += [float(data[wtsoft_start + q].split()[0])]
            wtsoft_timepos +=[float(data[wtsoft_start + q].split()[1])]           
            wtsoft_timeneg += [-1*float(data[wtsoft_start + q].split()[2])]
            wtsoft_rate += [float(data[wtsoft_start + q].split()[3])]
            wtsoft_ratepos += [float(data[wtsoft_start + q].split()[4])]
            wtsoft_rateneg += [float(data[wtsoft_start + q].split()[4])]

    if wtsoft_start != 0:
        
       wtsoft_timeerr = [wtsoft_timeneg,wtsoft_timepos]
       wtsoft_rateerr = [wtsoft_rateneg,wtsoft_ratepos]
       ax2.scatter(wtsoft_time,wtsoft_rate,color='blue',marker='.',zorder=1)
       ax2.errorbar(wtsoft_time,wtsoft_rate,xerr=wtsoft_timeerr,yerr=wtsoft_rateerr,capsize=0,fmt='o',color='blue',marker='.',zorder=1)

    if wtratio_start != 0:   
        for r in range (wtratio_stop - wtratio_start):
        
            wtratio_time += [float(data[wtratio_start + r].split()[0])]
            wtratio_timepos +=[float(data[wtratio_start + r].split()[1])]           
            wtratio_timeneg += [-1*float(data[wtratio_start + r].split()[2])]
            wtratio_rate += [float(data[wtratio_start + r].split()[3])]
            wtratio_ratepos += [float(data[wtratio_start + r].split()[4])]
            wtratio_rateneg += [float(data[wtratio_start + r].split()[4])]

    if wtratio_start != 0:
        
       wtratio_timeerr = [wtratio_timeneg,wtratio_timepos]
       wtratio_rateerr = [wtratio_rateneg,wtratio_ratepos]
       ax3.scatter(wtratio_time,wtratio_rate,color='blue',marker='.',zorder=1)
       ax3.errorbar(wtratio_time,wtratio_rate,xerr=wtratio_timeerr,yerr=wtratio_rateerr,capsize=0,fmt='o',color='blue',marker='.',zorder=1)

    if pcratio_stop == 0:
        pcratio_stop = endoffile

    """ PC Data """
    if pchard_start != 0:   
        for p in range (pchard_stop - pchard_start):
        
            pchard_time += [float(data[pchard_start + p].split()[0])]
            pchard_timepos +=[float(data[pchard_start + p].split()[1])]           
            pchard_timeneg += [-1*float(data[pchard_start + p].split()[2])]
            pchard_rate += [float(data[pchard_start + p].split()[3])]
            pchard_ratepos += [float(data[pchard_start + p].split()[4])]
            pchard_rateneg += [float(data[pchard_start + p].split()[4])]

    if pchard_start != 0:

       pchard_timeerr = [pchard_timeneg,pchard_timepos]
       pchard_rateerr = [pchard_rateneg,pchard_ratepos]
       ax1.scatter(pchard_time,pchard_rate,color='red',marker='.',zorder=1)
       ax1.errorbar(pchard_time,pchard_rate,xerr=pchard_timeerr,yerr=pchard_rateerr,capsize=0,fmt='o',color='red',marker='.',zorder=1)

    if pcsoft_start != 0:   
        for q in range (pcsoft_stop - pcsoft_start):
        
            pcsoft_time += [float(data[pcsoft_start + q].split()[0])]
            pcsoft_timepos +=[float(data[pcsoft_start + q].split()[1])]           
            pcsoft_timeneg += [-1*float(data[pcsoft_start + q].split()[2])]
            pcsoft_rate += [float(data[pcsoft_start + q].split()[3])]
            pcsoft_ratepos += [float(data[pcsoft_start + q].split()[4])]
            pcsoft_rateneg += [float(data[pcsoft_start + q].split()[4])]

    if pcsoft_start != 0:
        
       pcsoft_timeerr = [pcsoft_timeneg,pcsoft_timepos]
       pcsoft_rateerr = [pcsoft_rateneg,pcsoft_ratepos]
       ax2.scatter(pcsoft_time,pcsoft_rate,color='red',marker='.',zorder=1)
       ax2.errorbar(pcsoft_time,pcsoft_rate,xerr=pcsoft_timeerr,yerr=pcsoft_rateerr,capsize=0,fmt='o',color='red',marker='.',zorder=1)

    if pcratio_start != 0:   
        for r in range (pcratio_stop - pcratio_start):
        
            pcratio_time += [float(data[pcratio_start + r].split()[0])]
            pcratio_timepos +=[float(data[pcratio_start + r].split()[1])]           
            pcratio_timeneg += [-1*float(data[pcratio_start + r].split()[2])]
            pcratio_rate += [float(data[pcratio_start + r].split()[3])]
            pcratio_ratepos += [float(data[pcratio_start + r].split()[4])]
            pcratio_rateneg += [float(data[pcratio_start + r].split()[4])]

    if pcratio_start != 0:
        
       pcratio_timeerr = [pcratio_timeneg,pcratio_timepos]
       pcratio_rateerr = [pcratio_rateneg,pcratio_ratepos]
       ax3.scatter(pcratio_time,pcratio_rate,color='red',marker='.',zorder=1)
       ax3.errorbar(pcratio_time,pcratio_rate,xerr=pcratio_timeerr,yerr=pcratio_rateerr,capsize=0,fmt='o',color='red',marker='.',zorder=1)
       
    ax3.set_ylim(ymin=0)
    plt.draw()

    """ Creates arrays of the data """
    time_hard = wtsthard_time + wthard_time + pchard_time
    timepos_hard = wtsthard_timepos + wthard_timepos + pchard_timepos
    timeneg_hard = wtsthard_timeneg + wthard_timeneg + pchard_timeneg
    rate_hard = wtsthard_rate + wthard_rate + pchard_rate
    ratepos_hard = wtsthard_ratepos + wthard_ratepos + pchard_ratepos
    rateneg_hard = wtsthard_rateneg + wthard_rateneg + pchard_rateneg

    time_soft = wtstsoft_time + wtsoft_time + pcsoft_time
    timepos_soft = wtstsoft_timepos + wtsoft_timepos + pcsoft_timepos
    timeneg_soft = wtstsoft_timeneg + wtsoft_timeneg + pcsoft_timeneg
    rate_soft = wtstsoft_rate + wtsoft_rate + pchard_rate
    ratepos_soft = wtstsoft_ratepos + wtsoft_ratepos + pcsoft_ratepos
    rateneg_soft = wtstsoft_rateneg + wtsoft_rateneg + pcsoft_rateneg

    time_ratio = wtstratio_time + wtratio_time + pcratio_time
    timepos_ratio = wtstratio_timepos + wtratio_timepos + pcratio_timepos
    timeneg_ratio = wtstratio_timeneg + wtratio_timeneg + pcratio_timeneg
    rate_ratio = wtstratio_rate + wtratio_rate + pcratio_rate
    ratepos_ratio = wtstratio_ratepos + wtratio_ratepos + pcratio_ratepos
    rateneg_ratio = wtstratio_rateneg + wtratio_rateneg + pcratio_rateneg

    return(time_hard,timepos_hard,timeneg_hard,rate_hard,ratepos_hard,rateneg_hard,
           time_soft,timepos_soft,timeneg_soft,rate_soft,ratepos_soft,rateneg_soft,
           time_ratio,timepos_ratio,timeneg_ratio,rate_ratio,ratepos_ratio,rateneg_ratio)

############################################
##### Refine Afterglow Data for Models #####
############################################

def refine_afterglow_data(afterglow_file):
    """ Takes afterglow file data and removes any flaring (manually) so the
    data can be used to fit decay models """
    time,timepos,timeneg,rate,ratepos,rateneg = xray_afterglow_plot(afterglow_file)
    time_min = 0
    time_max = len(time) - 1
    rateerr = (ratepos+rateneg)/2
    
    """ Ask if there is any flaring present """
    flaring = float(input('How many flares are present? '))

    if flaring >= 1:
        time_fstart1 = float(input('\nTime first flaring started: '))
        time_fstop1 = float(input('Time first flaring stopped: '))

        flare_mask1 = np.ma.masked_inside(time,time_fstart1,time_fstop1)
        time = time[flare_mask1 != False]
        rate = rate[flare_mask1 != False]
        time_min = 0
        time_max = len(time) - 1
        rateerr = rateerr[flare_mask1 != False]
        plt.axvline(time_fstart1,color='grey')
        plt.axvline(time_fstop1,color='grey')
        plt.axvspan(time_fstart1,time_fstop1,alpha=0.5,color='grey')

        if flaring >= 2:
            time_fstart2 = float(input('\nTime second flaring started: '))
            time_fstop2 = float(input('Time second flaring stopped: '))

            flare_mask2 = np.ma.masked_inside(time,time_fstart2,time_fstop2)
            time = time[flare_mask2 != False]
            rate = rate[flare_mask2 != False]
            time_min = 0
            time_max = len(time) - 1
            rateerr = rateerr[flare_mask2 != False]
            plt.axvline(time_fstart2,color='grey')
            plt.axvline(time_fstop2,color='grey')
            plt.axvspan(time_fstart2,time_fstop2,alpha=0.5,color='grey')

            if flaring >= 3:
                time_fstart3 = float(input('\nTime third flaring started: '))
                time_fstop3 = float(input('Time third flaring stopped: '))

                flare_mask3 = np.ma.masked_inside(time,time_fstart3,time_fstop3)
                time = time[flare_mask3 != False]
                rate = rate[flare_mask3 != False]
                time_min = 0
                time_max = len(time) - 1
                rateerr = rateerr[flare_mask3 != False]
                plt.axvline(time_fstart3,color='grey')
                plt.axvline(time_fstop3,color='grey')
                plt.axvspan(time_fstart3,time_fstop3,alpha=0.5,color='grey')
                
    return(time,rate,rateerr,time_min,time_max)

##################################
##### Choose Power Law Model #####
##################################

def choose_model():
    breaks = float(input('\nHow many breaks in the afterglow decay? '))
    power_laws = breaks + 1
    print('\nPlease input initial estimated parameters:')
    
    if power_laws == 1:
        def func(x,n,a):
            return n * (x**a)

        n = float(input('Norm: '))
        a = float(input('Slope: '))
        guess = [n,a]
        param_names = ['norm','a']
        
    if power_laws == 2:
        def func(x,n,a1,xb,a2):
            return n * ((x**a1)*(x<=xb) + ((xb**(a1-a2))*(x**a2))*(x>xb))

        n = float(input('Norm: '))
        a1 = float(input('First Slope: '))
        xb = float(input('Time Break: '))
        a2 = float(input('Second Slope: '))
        guess = [n,a1,xb,a2]
        param_names = ['norm','a1','tb1','a2']
        
    if power_laws == 3:    
        def func(x,n,a1,xb1,a2,xb2,a3):
            return n * (((x**a1)*(x<=xb1)) + ((xb1**(a1-a2))*(x**a2)*(x>xb1)*(x<xb2)) + ((xb1**(a1-a2))*(xb2**(a2-a3))*(x**a3)*(x>=xb2)))

        n = float(input('Norm: '))
        a1 = float(input('First Slope: '))
        xb1 = float(input('First Time Break: '))
        a2 = float(input('Second Slope: '))
        xb2 = float(input('Second Time Break: '))
        a3 = float(input('Third Slope: '))
        guess = [n,a1,xb1,a2,xb2,a3]
        param_names = ['norm','a1','tb1','a2','tb2','a3']
        
    return(func,guess,param_names)
      
############################################
##### Scipy Optimize Power Law Fitting #####
############################################

def decay_fitting(afterglow_file):
    """ Takes the data from refined afterglow data, chooses a power law model,
    optimises the parameters using scipy curve fit and plots the data. It also
    removes any flaring if necessary """
    time,rate,rateerr,time_min,time_max = refine_afterglow_data(afterglow_file)
            
    """ Choose a power law with initial parameters and fit and optimize
    the curve (least squares) """
    func,guess,param_names = choose_model()   
    params, pcov = curve_fit(func,xdata=time,ydata=rate,p0=guess,
                            sigma=rateerr)
    
    perrs = np.sqrt(np.diag(pcov))

    """ Calculate Chi Squared Value for Model """
    dof = (len(time) - len(params))
    chi_squared = np.sum((func(time,*params) - rate)**2/rate)
    reduced_chi_squared = chi_squared/dof
    
    """ Plot the data and model using matplotlib """
    x = np.linspace(time[time_min],time[time_max],10000)
    modely = func(x,*params)
    plt.plot(x,modely,color='black',zorder=2)
    print('')

    for i in range(0,len(params)):
        print("""{0} = {1} (± {2})""".format(param_names[i],params[i],perrs[i]))
    
    print('\nChi Squared Value = ',chi_squared)
    print('Degrees of Freedom = ',dof)
    print('Reduced Chi Squared Value = ',reduced_chi_squared)
    
    return(params,perrs,chi_squared,dof)

#########################################
##### Single Power Law MCMC fitting #####
#########################################

def pl1_mc(afterglow_file):
    """ Takes the data from refined afterglow data and runs a Markov-Chain Monte
    Carlo to find the parameters of the fit and plots the data. It also
    removes any flaring if necessary """
    time,rate,rateerr,time_min,time_max = refine_afterglow_data(afterglow_file)
    x = time
    y = rate
    yerr = rateerr

    """ Give initial parameters """
    norm = float(input('Norm: '))
    slope = float(input('Slope: '))

    def prior(p):
        """ Define boundaries that the walkers should not leave """
        n, a = p
        if -10 < a < 10 and n > 0:
            return 0.0
        return -np.inf

    def err_func(p, x, y, yerr):
        """ Log-likelihood normal distribution function """
        n, a = p
        model = n * (x**a)
        sigma = yerr
        return -0.5*((np.sum((((y-model)/sigma)**2) + np.log(sigma**2))) + len(y)*np.log(2*np.pi))
    
    def prob(p, x, y, yerr):
        """ Probability function for MCMC """
        pr = prior(p)
        if not np.isfinite(pr):
            return -np.inf
        return pr + err_func(p, x, y, yerr)

    """ Set up sampler and run the MCMC. Remove burn in period and
    flatten samples array """
    guess = [norm,slope]
    nwalkers, ndim, nsteps = 100, 2, 10000
    p0 = [guess+(0.01*np.random.randn(ndim)) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, prob, args=(x, y, yerr))

    print('Running MCMC...')
    sampler.run_mcmc(p0, nsteps)
    print('MCMC Done!\n')
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

    """ Print results and uncretainties (16-84 quartile) """
    n, a = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
    print('norm = ',n[0],'(+',n[1],'-',n[2],')')
    print('a = ',a[0],'(+',a[1],'-',a[2],')')

    """ Plot Optimised Model over Data """
    xx = np.linspace(x[time_min],x[time_max],10000)
    ymod = n[0]*(xx**a[0])
    plt.plot(xx,ymod,color='black')

    """ Calculate and print Chi-squared values """
    dof = (len(x) - ndim)
    chi_squared = np.sum(((n[0]*(x**a[0])) - y)**2/y)
    reduced_chi_squared = chi_squared/dof
    
    print('\nChi Squared Value = ',chi_squared)
    print('Degrees of Freedom = ',dof)
    print('Reduced Chi Squared Value = ',reduced_chi_squared)

    return(n,a,chi_squared,dof)

#########################################
##### Broken Power Law MCMC fitting #####
#########################################

def pl2_mc(afterglow_file):
    """ Takes the data from refined afterglow data and runs a Markov-Chain Monte
    Carlo to find the parameters of the fit and plots the data. It also
    removes any flaring if necessary """
    time,rate,rateerr,time_min,time_max = refine_afterglow_data(afterglow_file)
    x = time
    y = rate
    yerr = rateerr

    """ Give initial parameters """
    norm = float(input('Norm: '))
    slope1 = float(input('First Slope: '))
    breakpoint = float(input('Breakpoint: '))
    slope2 = float(input('Second Slope: '))

    def prior(p):
        """ Define boundaries that the walkers should not leave """
        n, a1, xb, a2 = p
        if n > 0 and -7 < a1 < 1 and -7 < a2 < 1 and xb > 0:
            return 0.0
        return -np.inf

    def err_func(p, x, y, yerr):
        """ Log-likelihood normal distribution function """
        n, a1, xb, a2 = p
        model = n * ((x**a1)*(x<=xb) + ((xb**(a1-a2))*(x**a2))*(x>xb))
        sigma = yerr
        return -0.5*((np.sum((((y-model)/sigma)**2) + np.log(sigma**2))) + len(y)*np.log(2*np.pi))
    
    def prob(p, x, y, yerr):
        """ Probability function for MCMC """
        pr = prior(p)
        if not np.isfinite(pr):
            return -np.inf
        return pr + err_func(p, x, y, yerr)

    """ Set up sampler and run the MCMC. Remove burn in period and
    flatten samples array """
    guess = [norm,slope1,breakpoint,slope2]
    nwalkers, ndim, nsteps = 100, 4, 10000
    p0 = [guess+(0.01*np.random.randn(ndim)) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, prob, args=(x, y, yerr))

    print('Running MCMC...')
    sampler.run_mcmc(p0, nsteps)
    print('MCMC Done!\n')
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

    """ Print results and uncretainties (16-84 quartile) """
    n, a1, xb, a2 = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
    print('norm = ',n[0],'(+',n[1],'-',n[2],')')
    print('a1 = ',a1[0],'(+',a1[1],'-',a1[2],')')
    print('xb = ',xb[0],'(+',xb[1],'-',xb[2],')')
    print('a2 = ',a2[0],'(+',a2[1],'-',a2[2],')')

    """ Plot Optimised Model over Data """
    xx = np.linspace(time[time_min],time[time_max],10000)
    ymod = n[0] * ((xx**a1[0])*(xx<=xb[0]) + (xb[0]**(a1[0]-a2[0])*(xx**a2[0])*(xx>xb[0])))
    plt.plot(xx,ymod,color='black')

    """ Calculate and print Chi-squared values """
    dof = (len(time) - ndim)
    chi_squared = np.sum(((n[0] * ((x**a1[0])*(x<=xb[0]) + (xb[0]**(a1[0]-a2[0])*(x**a2[0])*(x>xb[0])))) - y)**2/y)
    reduced_chi_squared = chi_squared/dof
  
    print('\nChi Squared Value = ',chi_squared)
    print('Degrees of Freedom = ',dof)
    print('Reduced Chi Squared Value = ',reduced_chi_squared)

    return(n,a1,xb,a2,chi_squared,dof)

################################################
##### Doubly Broken Power Law MCMC fitting #####
################################################

def pl3_mc(afterglow_file):
    """ Takes the data from refined afterglow data and runs a Markov-Chain Monte
    Carlo to find the parameters of the fit and plots the data. It also
    removes any flaring if necessary """
    time,rate,rateerr,time_min,time_max = refine_afterglow_data(afterglow_file)
    x = time
    y = rate
    yerr = rateerr

    """ Give initial parameters """
    norm = float(input('Norm: '))
    slope1 = float(input('First Slope: '))
    breakpoint1 = float(input('Breakpoint: '))
    slope2 = float(input('Second Slope: '))
    breakpoint2 = float(input('Second Breakpoint: '))
    slope3 = float(input('Third Slope: '))

    def prior(p):
        """ Define boundaries that the walkers should not leave """
        n, a1, xb1, a2, xb2, a3 = p
        if n > 0 and -7 < a1 < 1 and -7 < a2 < 1 and -7 < a3 < 1 and xb1 > 0 and xb2 > 0:
            return 0.0
        return -np.inf

    def err_func(p, x, y, yerr):
        """ Log-likelihood normal distribution function """
        n, a1, xb1, a2, xb2, a3 = p
        model = n * (((x**a1)*(x<=xb1)) + ((xb1**(a1-a2))*(x**a2)*(x>xb1)*(x<xb2)) + ((xb1**(a1-a2))*(xb2**(a2-a3))*(x**a3)*(x>=xb2)))
        sigma = yerr
        return -0.5*((np.sum((((y-model)/sigma)**2) + np.log(sigma**2))) + len(y)*np.log(2*np.pi))
    
    def prob(p, x, y, yerr):
        """ Probability function for MCMC """
        pr = prior(p)
        if not np.isfinite(pr):
            return -np.inf
        return pr + err_func(p, x, y, yerr)

    """ Set up sampler and run the MCMC. Remove burn in period and
    flatten samples array """
    guess = [norm,slope1,breakpoint1,slope2,breakpoint2,slope3]
    nwalkers, ndim, nsteps = 100, 6, 10000
    p0 = [guess+(0.01*np.random.randn(ndim)) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, prob, args=(x, y, yerr))

    print('Running MCMC...')
    sampler.run_mcmc(p0, nsteps)
    print('MCMC Done!\n')
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

    """ Print results and uncretainties (16-84 quartile) """
    n, a1, xb1, a2, xb2, a3 = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
    print('norm = ',n[0],'(+',n[1],'-',n[2],')')
    print('a1 = ',a1[0],'(+',a1[1],'-',a1[2],')')
    print('xb1 = ',xb1[0],'(+',xb1[1],'-',xb1[2],')')
    print('a2 = ',a2[0],'(+',a2[1],'-',a2[2],')')
    print('xb2 = ',xb2[0],'(+',xb2[1],'-',xb2[2],')')
    print('a3 = ',a3[0],'(+',a3[1],'-',a3[2],')')

    """ Plot Optimised Model over Data """
    xx = np.linspace(time[time_min],time[time_max],10000)
    ymod = n[0] * ((xx**a1[0])*(xx<=xb1[0]) + (xb1[0]**(a1[0]-a2[0])*(xx**a2[0])*(xb1[0]<xx)*(xx<xb2[0])) + ((xb1[0]**(a1[0]-a2[0]))*(xb2[0]**(a2[0]-a3[0]))*(xx**a3[0])*(xx>=xb2[0])))
    plt.plot(xx,ymod,color='black')

    """ Calculate and print Chi-squared values """
    dof = (len(time) - ndim)
    chi_squared = np.sum(((n[0] * ((x**a1[0])*(x<=xb1[0]) + (xb1[0]**(a1[0]-a2[0])*(x**a2[0])*(xb1[0]<x)*(x<xb2[0])) + ((xb1[0]**(a1[0]-a2[0]))*(xb2[0]**(a2[0]-a3[0]))*(x**a3[0])*(x>=xb2[0])))) - y)**2/y)
    reduced_chi_squared = chi_squared/dof
  
    print('\nChi Squared Value = ',chi_squared)
    print('Degrees of Freedom = ',dof)
    print('Reduced Chi Squared Value = ',reduced_chi_squared)

    return(n,a1,xb1,a2,xb2,a3,chi_squared,dof)
