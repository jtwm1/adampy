import os
import numpy as np
import math
import sys

#######################################
##### Get Polarisation Parameters #####
#######################################

def get_pol_params(folder_path):
    """ Load relevant files """
    file_name_ord0 = 'angle0_ord.1'
    file_name_ord22 = 'angle225_ord.1'
    file_name_ord45 = 'angle45_ord.1'
    file_name_ord67 = 'angle675_ord.1'
    file_name_ext0 = 'angle0_exord.1'
    file_name_ext22 = 'angle225_exord.1'
    file_name_ext45 = 'angle45_exord.1'
    file_name_ext67 = 'angle675_exord.1'

    ordin_0 = os.path.join(folder_path,file_name_ord0)
    ordin_22 = os.path.join(folder_path,file_name_ord22)
    ordin_45 = os.path.join(folder_path,file_name_ord45)
    ordin_67 = os.path.join(folder_path,file_name_ord67)
    extra_0 = os.path.join(folder_path,file_name_ext0)
    extra_22 = os.path.join(folder_path,file_name_ext22)
    extra_45 = os.path.join(folder_path,file_name_ext45)
    extra_67 = os.path.join(folder_path,file_name_ext67)

    """ Defines two lists of ordinary and extra ordinary files to extract data """
    ordinary_beam = [ordin_0,ordin_22,ordin_45,ordin_67]
    extra_beam = [extra_0,extra_22,extra_45,extra_67]

    def beam_data(beam_angle):
        """ Extracts data for all targets per angle of selected beam """
        total_data = {}
    
        f = open(beam_angle,'r')
        data = np.genfromtxt(beam_angle,delimiter='',dtype=float,skip_header=2,
                             usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13])
    
        if np.shape(data) == (13,):

            x_centre = data[0]
            if np.isnan(x_centre) == True:
                print('X coord of source in',beam_angle,'contains INDEF!!!')
            y_centre = data[1]
            if np.isnan(y_centre) == True:
                print('Y coord of source in',beam_angle,'contains INDEF!!!')
            mag = data[2]
            if np.isnan(mag) == True:
                print('Magnitude of source in',beam_angle,'contains INDEF!!!')
            mag_err = data[3]
            if np.isnan(mag_err) == True:
                print('Mag Error of source in',beam_angle,'contains INDEF!!!')
            perr = data[4]
            msky = data[5]
            if np.isnan(msky) == True:
                print('M_sky of source in',beam_angle,'contains INDEF!!!')
            stdev = data[6]
            if np.isnan(stdev) == True:
                print('Standard Deviation of source in',beam_angle,'contains INDEF!!!')
            nsky = data[7]
            if np.isnan(nsky) == True:
                print('N_sky of source in',beam_angle,'contains INDEF!!!')
            nsrej = data[8]
            if np.isnan(nsrej) == True:
                print('Nsreg of source in',beam_angle,'contains INDEF!!!')
            serr = data[9]
            area = data[10]
            if np.isnan(area) == True:
                print('Area of source in',beam_angle,'contains INDEF!!!')
            flux = data[11]
            if np.isnan(flux) == True:
                print('Flux of source in',beam_angle,'contains INDEF!!!')
            totsum = data[12]
            if np.isnan(totsum) == True:
                print('Totsum of source in',beam_angle,'contains INDEF!!!')
            name = 'Source 0'
            total_data[name] = {'x': x_centre, 'y': y_centre, 'mag': mag, 'mag err': mag_err, 'p err': perr, 'm sky': msky, 'st dev': stdev,
                                'n sky': nsky, 'nsrej': nsrej, 'serr': serr, 'area': area, 'flux': flux, 'sum': totsum}
        
        else:
            for i in range(0,len(data),1):
            
                x_centre = data[i,0]
                if np.isnan(x_centre) == True:
                    print('X coord of source',[i],'in',beam_angle,'contains INDEF!!!')
                y_centre = data[i,1]
                if np.isnan(y_centre) == True:
                    print('Y coord of source',[i],'in',beam_angle,'contains INDEF!!!')
                mag = data[i,2]
                if np.isnan(mag) == True:
                    print('Magnitude of source',[i],'in',beam_angle,'contains INDEF!!!')
                mag_err = data[i,3]
                if np.isnan(mag_err) == True:
                    print('Mag Error of source',[i],'in',beam_angle,'contains INDEF!!!')
                perr = data[i,4]
                msky = data[i,5]
                if np.isnan(msky) == True:
                    print('M_sky of source',[i],'in',beam_angle,'contains INDEF!!!')
                stdev = data[i,6]
                if np.isnan(stdev) == True:
                    print('Standard Deviation of source',[i],'in',beam_angle,'contains INDEF!!!')
                nsky = data[i,7]
                if np.isnan(nsky) == True:
                    print('N_sky of source',[i],'in',beam_angle,'contains INDEF!!!')
                nsrej = data[i,8]
                if np.isnan(nsrej) == True:
                    print('Nsreg of source',[i],'in',beam_angle,'contains INDEF!!!')
                serr = data[i,9]
                area = data[i,10]
                if np.isnan(area) == True:
                    print('Area of source',[i],'in',beam_angle,'contains INDEF!!!')
                flux = data[i,11]
                if np.isnan(flux) == True:
                    print('Flux of source',[i],'in',beam_angle,'contains INDEF!!!')
                totsum = data[i,12]
                if np.isnan(totsum) == True:
                    print('Totsum of source',[i],'in',beam_angle,'contains INDEF!!!')
                name = 'Source '+ str(i)
                total_data[name] = {'x': x_centre, 'y': y_centre, 'mag': mag, 'mag err': mag_err, 'p err': perr, 'm sky': msky, 'st dev': stdev,
                                    'n sky': nsky, 'nsrej': nsrej, 'serr': serr, 'area': area, 'flux': flux, 'sum': totsum}
        
        f.close()
        return total_data

    def target_list(ordin_data_0):
        """ Creates a target list with main target and relevant
        number of field stars """
        target_list = []

        if np.shape(ordin_data_0) == (13,):
            target_list.append('Source 0')

        else:
            for i in range(0,len(ordin_data_0),1):
                name = 'Source '+ str(i)
                target_list.append(name)

        return target_list

    """ Data from file for each angle of each beam stored in these dictionaries
    and target list stored in an array """
    ordin_data_0 = beam_data(ordinary_beam[0])
    ordin_data_22 = beam_data(ordinary_beam[1])
    ordin_data_45 = beam_data(ordinary_beam[2])
    ordin_data_67 = beam_data(ordinary_beam[3])

    extra_data_0 = beam_data(extra_beam[0])
    extra_data_22 = beam_data(extra_beam[1])
    extra_data_45 = beam_data(extra_beam[2])
    extra_data_67 = beam_data(extra_beam[3])

    target_list = target_list(ordin_data_0)

    """ Ensure all angles in both ordinary and extraordinary beams have the
    same number of sources """
    if (len(ordin_data_0) or len(ordin_data_22) or len(ordin_data_45) or len(ordin_data_67)) != (len(extra_data_0) != len(extra_data_22) or len(extra_data_45) or len(extra_data_67)):
        print('One or more data files have unequal numbers of sources!')
        sys.exit()

    def flux_error(beam_info,target_list):
        """ Calculates the flux uncertainty for each source per angle of each beam
        using data gathered above """
        flux_error = []
        gain = 1.1
        k = 1
        nd = 1
        eta = 1
    
        for i in range(0,len(target_list),1):
            flux_err1 = beam_info[target_list[i]]['flux']/(gain*eta*nd)
            flux_err2 = beam_info[target_list[i]]['area']* beam_info[target_list[i]]['st dev']*beam_info[target_list[i]]['st dev']
            flux_err3 = (k/beam_info[target_list[i]]['n sky'])*(beam_info[target_list[i]]['area']*beam_info[target_list[i]]['st dev'])**2
        
            flux_error_calc = math.sqrt(flux_err1 + flux_err2 + flux_err3)
            flux_error.append(flux_error_calc)

        return flux_error

    """ Flux errors stored in following arrays """
    ordin_fluxerr_0 = flux_error(ordin_data_0,target_list)
    ordin_fluxerr_22 = flux_error(ordin_data_22,target_list)
    ordin_fluxerr_45 = flux_error(ordin_data_45,target_list)
    ordin_fluxerr_67 = flux_error(ordin_data_67,target_list)

    extra_fluxerr_0 = flux_error(extra_data_0,target_list)
    extra_fluxerr_22 = flux_error(extra_data_22,target_list)
    extra_fluxerr_45 = flux_error(extra_data_45,target_list)
    extra_fluxerr_67 = flux_error(extra_data_67,target_list)

    def norm_flux(ordin_beam,extra_beam,ordin_fluxerr,extra_fluxerr,target_list):
        """ Calculates the normalised flux per angle for each beam and the error on the flux """
        norm_flux_value = []
        norm_flux_err = []
    
        for i in range(0,len(target_list),1):
            norm_flux = (ordin_beam[target_list[i]]['flux']-extra_beam[target_list[i]]['flux'])/(ordin_beam[target_list[i]]['flux']+extra_beam[target_list[i]]['flux'])
            norm_flux_value.append(norm_flux)
            a = math.sqrt((ordin_fluxerr[i]**2)+(extra_fluxerr[i]**2))
            norm_flux_err.append(norm_flux*math.sqrt(((a/(ordin_beam[target_list[i]]['flux']-extra_beam[target_list[i]]['flux']))**2)
                                                     +((a/(ordin_beam[target_list[i]]['flux']+extra_beam[target_list[i]]['flux']))**2)))

        return(norm_flux_value,norm_flux_err)

    """ Normalised flux values and errors stored in fllowing arrays """
    norm_flux_0,norm_flux_err_0 = norm_flux(ordin_data_0,extra_data_0,ordin_fluxerr_0,extra_fluxerr_0,target_list)
    norm_flux_22,norm_flux_err_22 = norm_flux(ordin_data_22,extra_data_22,ordin_fluxerr_22,extra_fluxerr_22,target_list)
    norm_flux_45,norm_flux_err_45 = norm_flux(ordin_data_45,extra_data_45,ordin_fluxerr_45,extra_fluxerr_45,target_list)
    norm_flux_67,norm_flux_err_67 = norm_flux(ordin_data_67,extra_data_67,ordin_fluxerr_67,extra_fluxerr_67,target_list)

    """ Coversion factor from degrees to radians """
    dtr = math.pi/180

    def pol_parameters(norm_flux_0,norm_flux_22,norm_flux_45,norm_flux_67,target_list):
        """ Calculates the measured Q, U, and P values of the objects """
        q_values = []
        u_values = []
        p_values = []

        for i in range(0,len(target_list),1):
            q = (0.5*norm_flux_0[i]*math.cos(4*0*dtr))+(0.5*norm_flux_22[i]*math.cos(4*22.5*dtr))+(0.5*norm_flux_45[i]*math.cos(4*45*dtr))+(0.5*norm_flux_67[i]*math.cos(4*67.5*dtr))
            u = (0.5*norm_flux_0[i]*math.sin(4*0*dtr))+(0.5*norm_flux_22[i]*math.sin(4*22.5*dtr))+(0.5*norm_flux_45[i]*math.sin(4*45*dtr))+(0.5*norm_flux_67[i]*math.sin(4*67.5*dtr))
            p = math.sqrt((q**2)+(u**2))
            theta = 0.5*math.atan(u/q)*(1/dtr)
            q_values.append(q)
            u_values.append(u)
            p_values.append(p)
        
        return(q_values,u_values,p_values)

    """ Q, U, and P values in following arrays """
    q_values,u_values,p_values = pol_parameters(norm_flux_0,norm_flux_22,norm_flux_45,norm_flux_67,target_list)

    """ Need to make more accurate from future measurements but it will do for now
    (Note that the curve is for R Band measurements so could be different for other
    filters in the optical) """
    par_ang = float(input('Parallactic Angle (degrees): '))
    curveq = 4.29118958*np.cos(np.deg2rad((2*par_ang + 81.91428479)))/100
    curveu = -4.3095486*np.cos(np.deg2rad((2*par_ang + 170.75451418)))/100

    real_q = [x - curveq for x in q_values]
    real_u = [x - curveu for x in u_values]
    real_p = []
    for i in range(0,len(target_list),1):
        real_p.append(math.sqrt((real_q[i]**2)+(real_u[i]**2)))
    
    def calc_theta(real_u,real_q):
        """ Calculates theta for all objects """
        theta_values = []

        for i in range(0,len(target_list),1):
            theta_values.append(0.5*math.atan(real_u[i]/real_q[i])*(1/dtr))

        return theta_values

    def position_angle(theta_values,real_q,real_u,target_list):
        """ Calculate proper position angles """
        corr_theta_values = []
    
        for i in range(0,len(target_list),1):

            if real_q[i] < 0:
                corr_theta_values.append(theta_values[i]+90)

            if real_q[i] > 0 and real_u[i] > 0:
                corr_theta_values.append(theta_values[i]+0)

            if real_q[i] > 0 and real_u[i] < 0:
                corr_theta_values.append(theta_values[i]+180)

        return corr_theta_values

    """ Theta value arrays """
    theta_values = calc_theta(real_u,real_q)
    corr_theta_values = position_angle(theta_values,real_q,real_u,target_list)

    def parameter_errors(norm_flux_err_0,norm_flux_err_22,norm_flux_err_45,norm_flux_err_67,real_p,ordin_data_0,extra_data_0,
                         ordin_data_22,extra_data_22,ordin_data_45,extra_data_45,ordin_data_67,extra_data_67,target_list):
        """ Calculate errors on Q, U and uses them to calculate error on P, Theta and
        also calculates simple P error from flux intensity SNR assuming all
        errors are gaussian (Patat & Romaniello, 2006 (section 3) """
        q_errors = []
        u_errors = []
        sig_p = []
        p_abs = []
        theta_errors = []

        for i in range(0,len(target_list),1):
            q_errors.append(math.sqrt(((0.5*norm_flux_err_0[i]*math.cos(0*dtr))**2)+((0.5*norm_flux_err_22[i]*math.cos(22.5*dtr))**2)+
                                      ((0.5*norm_flux_err_45[i]*math.cos(45*dtr))**2)+((0.5*norm_flux_err_67[i]*math.cos(67.5*dtr))**2)))
            u_errors.append(math.sqrt(((0.5*norm_flux_err_0[i]*math.sin(0*dtr))**2)+((0.5*norm_flux_err_22[i]*math.sin(22.5*dtr))**2)+
                                      ((0.5*norm_flux_err_45[i]*math.sin(45*dtr))**2)+((0.5*norm_flux_err_67[i]*math.sin(67.5*dtr))**2)))
            p_abs.append(1/(math.sqrt((ordin_data_0[target_list[i]]['flux']+extra_data_0[target_list[i]]['flux']+ordin_data_22[target_list[i]]['flux']
                                       +extra_data_22[target_list[i]]['flux']+ordin_data_45[target_list[i]]['flux']+extra_data_45[target_list[i]]['flux']
                                       +ordin_data_67[target_list[i]]['flux']+extra_data_67[target_list[i]]['flux'])/4)))


        for j in range(0,len(target_list),1):
            sig_p.append((0.5/real_p[j])*math.sqrt((((2*(real_q[j])**2)*(q_errors[j]/real_q[j]))**2)+
                                                   (((2*(real_u[j])**2)*(u_errors[j]/real_u[j]))**2)))

        for k in range(0,len(target_list),1):
            theta_errors.append(p_abs[k]/(2*real_p[k]*dtr))

        return(q_errors,u_errors,sig_p,p_abs,theta_errors)   
    
    """ Store Q, U and P, Theta and simplified P error values in following arrays """
    q_errors,u_errors,sig_p,p_abs,theta_errors = parameter_errors(norm_flux_err_0,norm_flux_err_22,norm_flux_err_45,norm_flux_err_67,real_p,ordin_data_0,extra_data_0,
                                        ordin_data_22,extra_data_22,ordin_data_45,extra_data_45,ordin_data_67,extra_data_67,target_list)

    def estimated_polarisation(ordin_data_0,extra_data_0,ordin_data_22,extra_data_22,ordin_data_45,extra_data_45,ordin_data_67,extra_data_67,ordin_fluxerr_0,ordin_fluxerr_22,
                               ordin_fluxerr_45,ordin_fluxerr_67,extra_fluxerr_0,extra_fluxerr_22,extra_fluxerr_45,extra_fluxerr_67,real_p,p_abs,target_list):
        """ Calculate etap from snr of flux and then use W-K to estimate
        polarisation where applicable """
        snr_f0 = []
        snr_f22 = []
        snr_f45 = []
        snr_f67 = []
        snr_fav = []
        eta = []
        wk_est = []
        p_corr = []

        for i in range(0,len(target_list),1):
            snr_f0.append((ordin_data_0[target_list[i]]['flux']+extra_data_0[target_list[i]]['flux'])/math.sqrt((ordin_fluxerr_0[i]**2)+(extra_fluxerr_0[i]**2)))
            snr_f22.append((ordin_data_22[target_list[i]]['flux']+extra_data_22[target_list[i]]['flux'])/math.sqrt((ordin_fluxerr_22[i]**2)+(extra_fluxerr_22[i]**2)))
            snr_f45.append((ordin_data_45[target_list[i]]['flux']+extra_data_45[target_list[i]]['flux'])/math.sqrt((ordin_fluxerr_45[i]**2)+(extra_fluxerr_45[i]**2)))
            snr_f67.append((ordin_data_67[target_list[i]]['flux']+extra_data_67[target_list[i]]['flux'])/math.sqrt((ordin_fluxerr_67[i]**2)+(extra_fluxerr_67[i]**2)))

        for j in range(0,len(target_list),1):
            snr_fav.append((snr_f0[j]+snr_f22[j]+snr_f45[j]+snr_f67[j])/4)

        for k in range(0,len(target_list),1):
            eta.append(real_p[k]*snr_fav[k])
    
        for l in range(0,len(target_list),1):
    
            if eta[l] >= 2:
                wk = math.sqrt(1-((p_abs[l]/real_p[l])**2))
                wk_est.append(wk)
                p_corr.append(real_p[l]*wk)
            
            elif eta[l] < 2:
                wk_est.append(str('  N/A  '))
                print(target_list[l],'has eta < 2 !!!')
                p_corr.append(0)
            
        return(eta,wk_est,p_corr)

    """ Eta, W-K estimator and corrected P values stored in following arrays """
    eta_values,wk_est_values,p_corr_values = estimated_polarisation(ordin_data_0,extra_data_0,ordin_data_22,extra_data_22,ordin_data_45,extra_data_45,ordin_data_67,extra_data_67,ordin_fluxerr_0,ordin_fluxerr_22,
                                           ordin_fluxerr_45,ordin_fluxerr_67,extra_fluxerr_0,extra_fluxerr_22,extra_fluxerr_45,extra_fluxerr_67,real_p,p_abs,target_list)

    """ Convert polarisations into percentages (*100) and round to proportionate
    significant decimal places """
    q_values = [x * 100 for x in q_values]
    real_q = [x * 100 for x in real_q]
    q_errors = [x * 100 for x in q_errors]
    u_values = [x * 100 for x in u_values]
    real_u = [x * 100 for x in real_u]
    u_errors = [x * 100 for x in u_errors]
    p_values = [x * 100 for x in p_values]
    real_p = [x * 100 for x in real_p]
    p_corr_values = [x * 100 for x in p_corr_values]
    for i in range(0,len(target_list),1):
        if p_corr_values[i] == 0:
            p_corr_values[i] = str('  ???  ')
        
    p_abs = [x * 100 for x in p_abs]

    """ Round all values to 5 significant Figures for ease when looking at printed
    results (Final result arrays!) """
    q = []
    q_r = []
    q_err = []
    u = []
    u_r = []
    u_err = []
    p = []
    p_r = []
    p_corr = []
    p_err = []
    theta = []
    theta_err = []
    eta = []
    wk_est = []

    for i in range(0,len(target_list),1):
        q.append(round(q_values[i],5))
        q_r.append(round(real_q[i],5))
        q_err.append(round(q_errors[i],5))
        u.append(round(u_values[i],5))
        u_r.append(round(real_u[i],5))
        u_err.append(round(u_errors[i],5))
        p.append(round(p_values[i],5))
        p_r.append(round(real_p[i],5))
        if p_corr_values[i] != '  ???  ':
            p_corr.append(round(p_corr_values[i],5))
        if p_corr_values[i] == '  ???  ':
            p_corr.append(p_corr_values[i])
        p_err.append(round(p_abs[i],5))
        theta.append(round(corr_theta_values[i],5))
        theta_err.append(round(theta_errors[i],5))
        eta.append(round(eta_values[i],5))
        if wk_est_values[i] != '  N/A  ':
            wk_est.append(round(wk_est_values[i],5))
        if wk_est_values[i] == '  N/A  ':
            wk_est.append(wk_est_values[i])

    return(q,q_r,q_err,u,u_r,u_err,p,p_r,p_corr,p_err,theta,theta_err,eta,wk_est,target_list)

##############################################
##### Write Polarisation Results to File #####
##############################################

def params_to_file(folder_path):
    """ Create and write to file to save the results """
    q,q_r,q_err,u,u_r,u_err,p,p_r,p_corr,p_err,theta,theta_err,eta,wk_est,target_list = get_pol_params(folder_path)    
    orig_stdout = sys.stdout
    result_file=folder_path+'source_results.txt'
    resultf = open(result_file, "w")
    sys.stdout = resultf

    print('                           ### RESULTS ###                            ')
    print('----------------------------------------------------------------------')
    print('')
    print(' Qm(%)  Qr(%)  Q Err(%)  Um(%)  Ur(%)  U Err(%)')
    for i in range(0,len(target_list),1):
        print(q[i], q_r[i], q_err[i], u[i], u_r[i], u_err[i])
    print('')
    print(' Eta    W-K Est   Pm(%)  Pr(%)  Pcorr(%) P Err(%)')
    for j in range(0,len(target_list),1):
        print(eta[j], wk_est[j], p[j], p_r[j], p_corr[j], p_err[j])
    print('')
    print(' Ang    Ang Err')
    for k in range(0,len(target_list),1):
        print(theta[k], theta_err[k])
    print('')
    print('----------------------------------------------------------------------')

    sys.stdout = orig_stdout
    resultf.close()
