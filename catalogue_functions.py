import os
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import requests
import json
import sys
from lxml import html
from bs4 import BeautifulSoup

###################################
##### Swift GRB Position Info ##### 
###################################

def swift_info(swift_file_name):
    """ Gets the XRT info and saves it to a file or if webpage can't be found
    it uses the last updated file (Thanks to UK Swift Science Data Centre) """

    try:
        r = requests.get("http://www.swift.ac.uk/xrt_positions/?txt=1").text
    
        data_list = []

        for line in r.split('\n'):
            fields = line.split('|')

            if len(fields) < 5:

                continue
        
            grbname = fields[0]
            grbcoord = SkyCoord(fields[1],fields[2],unit=(u.hourangle, u.deg))
            grbra = float(grbcoord.ra.degree)
            grbdec = float(grbcoord.dec.degree)
            data = {'Name': grbname, 'RA': grbra, 'Dec': grbdec}
            data_list.append(data)

        with open(swift_file_name, 'w') as outfile:
            json.dump(data_list, outfile, indent=0)

    except requests.exceptions.RequestException as e:
        print('Cannot access Swift UK XRT positions webpage so using back up file!\n')
        with open(swift_file_name) as infile:
            data_list = json.load(infile)
            
    return data_list

###################################
##### Fermi GRB Position Info #####
###################################

def fermi_info(fermi_file_name):
    """ Gets the Fermi LAT info and saves it to a file or if webpage can't be
    found it uses the last updated file (Thanks to Fermi/GSFC Website) """
    
    try:
        r = requests.get("http://slac.stanford.edu/~glast/LATBA/PublicTableGRBs.xml").content

        data_list = []
        grballinfo = html.fromstring(r)

        for i in range(0,len(grballinfo),1):
            grbname = grballinfo[i][1].text
            grbra = grballinfo[i][5].text
            grbdec = grballinfo[i][6].text
            grbra = float(grbra)
            grbdec = float(grbdec)
            data = {'Name': grbname, 'RA': grbra, 'Dec': grbdec}
            data_list.append(data)

        with open(fermi_file_name, 'w') as outfile:
            json.dump(data_list, outfile, indent=0)

    except requests.exceptions.RequestException as e:
        print('Cannot access Stanford Fermi LAT GRB table webpage so using back up file!\n')
        with open(fermi_file_name) as infile:
            data_list = json.load(infile)
            
    return data_list

###################################
##### BATSE GRB Position Info ##### 
###################################

def batse_info(batse_file_name):
    """ Gets the BATSE info and saves it to a file or if webpage can't be
    found it uses the last updated file (Thanks to Gamma-Ray Astrophysics
    NSSTC webpage) """
    
    try:
        r = requests.get("http://gammaray.nsstc.nasa.gov/batse/grb/catalog/4b/tables/4br_grossc.basic").text
        data_list = []

        for line in r.split('\n'):
    
            if '#' not in line:
                fields = line.split()

                if len(fields) == 13:
                    grbname = 'GRB' + fields[2]
                    grbcoord = SkyCoord(fields[5],fields[6],unit=(u.hourangle, u.deg))
                    grbra = float(grbcoord.ra.degree)
                    grbdec = float(grbcoord.dec.degree)
                    data = {'Name': grbname, 'RA': grbra, 'Dec': grbdec}
                    data_list.append(data)

        with open(batse_file_name, 'w') as outfile:
            json.dump(data_list, outfile, indent=0)

    except requests.exceptions.RequestException as e:
        print('Cannot access NASA BATSE catalogue webpage so using back up file!\n')
        with open(batse_file_name) as infile:
            data_list = json.load(infile)
            
    return data_list

######################################
##### INTEGRAL GRB Position Info ##### 
######################################

def integral_info(integral_file_name):
    """ Gets the INTEGRAL info and saves it to a file or if webpage can't be found
    it uses the last updated file (Thanks to IBAS Results Site) """
    
    try:
        r = requests.get("http://ibas.iasf-milano.inaf.it/IBAS_Results.html").text
        root = html.fromstring(r)
        data_list = []
        c = 0
        rows = root.xpath('//table')[1].findall('.//tr')

        for row in rows:

            if c >= 1:
                tabledata = ([col.text for col in row.getchildren()])
                grbname = 'GRB' + tabledata[0]
                grbcoord = SkyCoord(tabledata[1],tabledata[2],unit=(u.hourangle, u.deg))
                grbra = float(grbcoord.ra.degree)
                grbdec = float(grbcoord.dec.degree)
                data_list.append({'Name': grbname, 'RA': grbra, 'Dec': grbdec})

            c+= 1

        with open(integral_file_name, 'w') as outfile:
            json.dump(data_list, outfile, indent=0)

    except requests.exceptions.RequestException as e:
        print('Cannot access IBAS results webpage so using back up file!\n')
        with open(integral_file_name) as infile:
            data_list = json.load(infile)
            
    return data_list

#####################################
##### Gaia Alerts Position Info ##### 
#####################################

def gaia_info(gaia_file_name):
    """ Gets the Gaia Alert info and saves it to a file or if webpage can't be
    found it uses the last updated file (ESA Gaia, DPAC and the Photometric Science
    Alerts Team) """
    
    try:
        r = requests.get("http://gsaweb.ast.cam.ac.uk/alerts/alerts.csv").text
        data_list = []
        c = 0
    
        for line in r.split('\n'):
            fields = line.split(',')

            if len(fields) < 3:

                continue
        
            c += 1
            
            if c >= 2:
                gaianame = fields[0]
                gaiara = fields[2]
                gaiadec = fields[3]
                data = {'Name': gaianame, 'RA': gaiara, 'Dec': gaiadec}
                data_list.append(data)

        with open(gaia_file_name, 'w') as outfile:
            json.dump(data_list, outfile, indent=0)
    
    except requests.exceptions.RequestException as e:
        print('Cannot access Gaia Transient Alerts webpage so using back up file!\n')
        with open(gaia_file_name) as infile:
            data_list = json.load(infile)
    
    return data_list

########################################
##### Galactic Novae Position Info ##### 
########################################

def galnovae_info(gaia_file_name):
    """ Gets the Galactic Novae info and saves it to a file or if webpage can't be
    found it uses the last updated file (Thanks to Project Pluto Website) """
    
    try:
        r = requests.get("http://projectpluto.com/galnovae/galnovae.txt").text
        data_list = []

        for line in r.split('\n'):

            if 'NOVA' in line:
                continue
    
            if line == '':
                continue

            lines = str(line)
            novaename = lines[0:13]
            if novaename == ('*            ' or '             '):
                novaename = 'Unknown Novae'
                
            novaeras = lines[44:56]
            novaedecs = lines[59:71]
            novaecoord = SkyCoord(novaeras,novaedecs,unit=(u.hourangle, u.deg))
            novaera = float(novaecoord.ra.degree)
            novaedec = float(novaecoord.dec.degree)
            data_list.append({'Name': novaename, 'RA': novaera, 'Dec': novaedec})
            
        with open(galnovae_file_name, 'w') as outfile:
            json.dump(data_list, outfile, indent=0)

    except requests.exceptions.RequestException as e:
        print('Cannot access Project Pluto Galactice Novae table webpage so using back up file!\n')
        with open(galnovae_file_name) as infile:
            data_list = json.load(infile)

    return data_list

####################################
##### Supernovae Position Info ##### 
####################################

def sn_info(sn_file_name):
    """ Gets the SN info and saves it to a file or if webpage can't be
    found it uses the last updated file (Thanks to IAC) """
    
    try:
        r = requests.get("http://www.cbat.eps.harvard.edu/lists/Supernovae.html").text
        soup = BeautifulSoup(r,'lxml')
        table = soup.find('pre')

        for script in table(["script", "style"]):
            script.extract()

        tabletext = table.get_text()
        c=0
        data_list = []

        for line in tabletext.split('\n'):
    
            if (line != '' and c >= 2):
                lines = str(line)
                ra = str(lines[87:98])
                dec = str(lines[99:110])
        
                if (ra and dec) != '           ':
                    snname = lines[0:6]
                    if snname == '      ':
                        snname = 'Unknown SN'
                        
                    sncoords = SkyCoord(lines[87:98],lines[99:110],unit=(u.hourangle,u.deg))
                    snra = float(sncoords.ra.degree)
                    sndec = float(sncoords.dec.degree)
                    data_list.append({'Name': snname, 'RA': snra, 'Dec': sndec})
            c += 1
            
        with open(sn_file_name, 'w') as outfile:
            json.dump(data_list, outfile, indent=0)

    except requests.exceptions.RequestException as e:
        print('Cannot access IAU SN table webpage so using back up file!\n')
        with open(sn_file_name) as infile:
            data_list = json.load(infile)

    return data_list

##################################
##### Magnetar Position Info ##### 
##################################

def magnetar_info(magnetar_file_name):
    """ Gets the magnetar info and saves it to a file or if webpage can't be found
    it uses the last updated file (Thanks to McGill Online Magnetar Catalogue) """

    try:
        r = requests.get("http://www.physics.mcgill.ca/~pulsar/magnetar/TabO1.txt").text
        data_list = []
        c = 0

        for line in r.split('\n'):
            fields = line.split('|')

            if (c >= 4 and len(fields) == 26):
                magname = fields[1]
                ra = str(fields[22])[0:12]
                dec = str(fields[23])[0:12]

                if (ra != '            ' and dec != '            '):
                    magcoords = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
                    magra = float(magcoords.ra.degree)
                    magdec = float(magcoords.dec.degree)
                    data_list.append({'Name': magname, 'RA': magra, 'Dec': magdec})
        
            c += 1
        
        with open(magnetar_file_name, 'w') as outfile:
            json.dump(data_list, outfile, indent=0)

    except requests.exceptions.RequestException as e:
        print('Cannot access McGill Magnetar table webpage so using back up file!\n')
        with open(magnetar_file_name) as infile:
            data_list = json.load(infile)
            
    return data_list

########################################
##### PESSTO Objects Position Info ##### 
########################################

def pessto_info(pessto_file_name):
    """ Gets the PESSTO Transient info and saves it to a file or if webpage can't be
    found it uses the last updated file (Thanks to PESSTO website) """
    
    try:
        r = requests.get("http://www.pessto.org/index.py?pageName=classifications").text
        soup = BeautifulSoup(r,'lxml')
        data_list = []
        c = 0
        table = soup.find('table')
        rows = table.findAll('tr')

        for row in rows:
    
            if c >= 1:
                cols = row.findAll('td')
                cols = [ele.text.strip() for ele in cols]
                data_list.append({'Name': cols[0], 'RA': cols[1], 'Dec': cols[2]})
                
            c += 1

        with open(pessto_file_name, 'w') as outfile:
            json.dump(data_list, outfile, indent=0)

    except requests.exceptions.RequestException as e:
        print('Cannot access PESSTO Transient table webpage so using back up file!\n')
        with open(pessto_file_name) as infile:
            data_list = json.load(infile)
            
    return data_list

##################################
##### Cross Catalogue Search ##### 
##################################

def catalogue_cross_ref(trigger_name,trigger_ra,trigger_dec,
                        obj_name,obj_ra,obj_dec,error_radius):
    ''' Convert error radius from arcseconds to degrees '''
    error_radius_deg = error_radius * 0.000277778

    """ Iterate over all triggers in the candidate file """
    for i in range (len(trigger_name)):
        a = trigger_ra[i] - error_radius_deg
        b = trigger_ra[i] + error_radius_deg
        c = trigger_dec[i] - error_radius_deg
        d = trigger_dec[i] + error_radius_deg

        obj_ra = obj_ra.astype(np.float)
        obj_dec = obj_dec.astype(np.float)
        
        mask1 = np.ma.masked_outside(obj_ra,a,b)
        obj_dec_m = obj_dec[mask1 != False]
        obj_ra_m = obj_ra[mask1 != False]
        obj_name_m = obj_name[mask1 != False]

        if len(obj_name_m) != 0:
            mask2 = np.ma.masked_outside(obj_dec_m,c,d)
            obj_dec_mm = obj_dec_m[mask2 != False]
            obj_ra_mm = obj_ra_m[mask2 != False]
            obj_name_mm = obj_name_m[mask2 != False]
            
            """ If there are objects within the error region print the number, names and coordinates """
            
            if len(obj_name_mm) != 0:
                print(trigger_name[i],'matched with',len(obj_name_mm),
                      'catalogue objects within', error_radius_as,'arc second(s)')
                print(trigger_name[i],'RA and Dec =',trigger_ra[i],trigger_dec[i],'\n')                

                for j in range(len(obj_name_mm)):
                    print(obj_name_mm[j],'RA and Dec =',obj_ra_mm[j],obj_dec_mm[j],'\n')
