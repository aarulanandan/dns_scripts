import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
from os import listdir
from os.path import isfile, join, basename, isdir
from numpy.fft import fft, fftfreq

#creates non-normalized .csv files for average and fluctuation values of each individual variable
    #variable directories (var_dirs) must be list of full file paths to location of extracted data files
    #variable names (var_names) must be list of identifiable characters/names for individual variables (eg. T or Temperature) in order in which directories are given in var_dirs
    #output directory (out_path) must be subdirectory in which individual directories for each term will be created and populated
def read_data(Xgrid,Ygrid,var_dirs,var_names,out_path):
    #example: read_data(4096,4096,['/localdata/scratch/Haripriya/Detonation_DNS/2_H2O2Ar_L16/RawVars/Temperature','/localdata/scratch/Haripriya/Detonation_DNS/2_H2O2Ar_L16/RawVars/CHrelease'],['T','Q'],'/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16')

    #resizing array if necessary
    Xsz = np.arange(Xgrid)
    Ysz = np.arange(Ygrid)

    nvars = len(var_names)

    if not os.path.isdir(os.path.join(out_path,'Averages')):
        os.mkdir(os.path.join(out_path,'Averages'))
    if not os.path.isdir(os.path.join(out_path,'Fluctuations')):
        os.mkdir(os.path.join(out_path,'Fluctuations'))

    #calculating averages and fluctuations of each term
    for n in range(0,nvars):

        dinfo = [f for f in os.listdir(var_dirs[n]) if os.path.isfile(os.path.join(var_dirs[n],f))]
        nfiles = len(dinfo)

        q_avg = 0
        for fi in range(0,nfiles):

            q_raw = np.asarray(pd.read_csv(os.path.join(var_dir[n],dinfo[fi]),header=None))
            q_raw2 = q_raw[Ysz,...]
            q = q_raw2[...,Xsz]

            q_avg = q_avg+(q/nfiles)

        filename = 'Averages/'+var_names[n][0]+'avg.csv'
        np.savetxt(os.path.join(out_path,filename),q_avg,delimiter=',')

        for fi in range(0,nfiles):

            q_raw = np.asarray(pd.read_csv(os.path.join(var_dirs[n],dinfo[fi]),header=None))
            q_raw2 = q_raw[Ysz,...]
            q = q_raw2[...,Xsz]

            q_fluc = q-q_avg

            if not os.path.isdir(os.path.join(out_path,'Fluctuations',var_names[n])):
                os.mkdir(os.path.join(out_path,'Fluctuations',var_names[n]))

            filename = 'Fluctuations/'+var_names[n]+'/'+var_names[n][0]+(''.join(filter(str.isdigit,str(dinfo[fi]))))+'.csv'
            np.savetxt(os.path.join(out_path,filename),q_fluc,delimiter=',')

#creates normalized .csv files for source terms
    #shock.csv must exist (find_shock function must be run already)
    #terms (terms) must be list in format ['Xavg*Yfluc',etc.]
    #variable names (var_names) must be list of identifiable characters/names that match those used in terms (eg. T or Q) and are given in the order in which average and fluctuation files are given in avg_files and fluc_files
    #Average files (avg_files) must be list of filenames with full path
    #Fluctuation files (fluc_files) must be list of paths to directories containing files
    #Average and fluctuation file paths must be listed in order in which they appear in terms
    #output directory (out_path) must be directory in which individual subdirectories for each term will be created
def source(Xgrid,Ygrid,cell_width,terms,var_names,avg_files,fluc_files,out_path):
    #example: source(4096,4096,1.45,['Tfluc*Qavg','Tfluc*Qfluc'],['T','Q'],['/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/Averages/Tavg.csv','/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/Averages/Qavg.csv'],['/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/Fluctuations/Temperature/','/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/Fluctuations/Heat_Release'],'/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16')

    nterms = len(terms)
    nvars = len(var_names)
    Xsz = np.arange(Xgrid)
    Ysz = np.arange(Ygrid)
    dinfon = []

    shock = np.asarray(pd.read_csv(os.path.join(out_path,'shock.csv')))
    (vars_norm,L_norm,t_norm) = normalize(Xgrid,Ygrid,avg_files,cell_width,'',shock,out_path)

    for m in range(0,nvars):

        dinfo = [f for f in os.listdir(fluc_files[m]) if os.path.isfile(os.path.join(fluc_files[m],f))]
        dinfon = np.append(dinfon,dinfo)

    nfiles = len(dinfon[0])

    for n in range(0,nterms):

        term_ind = terms[n].split('*')
        nterms_ind = len(term_ind)

        indx = []

        for m in range(0,nterms_ind):
            indx_m = var_names.index(term_ind[m][0])
            indx = np.append(indx,indx_m)

        os.mkdir(os.path.join(out_path+(''.join(term_ind))))
        
        #term n calculation
        for fi in range(0,nfiles):

            term_n = np.ones([Ygrid,Xgrid])

            for m in range(0,nterms_ind):
                if 'avg' in term_ind[m]:
                    q = np.asarray(pd.read_csv(avg_files[(indx[m])],header=None))
                    q_1 = q[Ysz,...]
                    term_q = q_1[...,Xsz]
                else:
                    q = np.asarray(pd.read_csv(os.path.join(fluc_files[(indx[m])],dinfon[(indx[m])][fi]),header=None))
                    q_1 = q[Ysz,...]
                    term_q = q_1[...,Xsz]

                term_n = term_n*(term_q/vars_norm[(indx[m])])

            filename = ''.join(term_ind)+'/'+(''.join(filter(str.isdigit,str(dinfo[fi]))))+'.csv'

            np.savetxt(os.path.join(out_path,filename),term_n,delimiter=',')

#Plots non-normalized spatial distribution of extracted quantities at each time step
    #Lengths (Xlen, Ylen) must be real dimensions of channel (HxW) in cm
    #Variable directories (var_dirs) must be list of full paths to directories containing extracted data files for each variable of interest
    #Variable names (var_names) must be list of identifiable characters or names (eg. T or Temperature) in order of given variable directories
    #Variable units (var_units) must be list of real units for each quantity
    #Output path (out_path) must be full file path to directory in which subdirectories 'plots' and individual subdirectories for each variable will be created and populated
def raw_plots(Xlen,Ylen,Xgrid,Ygrid,var_dirs,var_names,var_units,out_path):
    #example: raw_plots(16,4,4096,4096,['/localdata/scratch/Haripriya/Detonation_DNS/3_H2O2Ar_L13.65/RawVars/Temperature/','/localdata/scratch/Haripriya/Detonation_DNS/3_H2O2Ar_L13.65/RawVars/Pressure/'],['Temperature','Pressure'],['T [K]','P [bar]'],'/localdata/scratch/adhav/Source_Comp_3_H2O2Ar_L13.65/')

    Xsz = np.arange(Xgrid)
    Ysz = np.arange(Ygrid)
    nvars = len(var_names)

    X = (Xsz/Xgrid)*Xlen
    Y = (Ysz/Ygrid)*Ylen

    if not os.path.isdir(os.path.join(out_path,'plots')):
            os.mkdir(os.path.join(out_path,'plots'))

    for n in range(0,nvars):

        dinfo = [f for f in os.listdir(var_dirs[n]) if os.path.isfile(os.path.join(var_dirs[n],f))]
        nfiles = len(dinfo)

        os.mkdir(os.path.join(out_path,'plots',var_names[n]))

        for fi in range(0,nfiles):

            q_raw = np.asarray(pd.read_csv(os.path.join(var_dirs[n],dinfo[fi])))
            q1_raw = q_raw[Ysz,...]
            q = q1_raw[...,Xsz]
           
            q_fig,q_ax = plt.subplots()
            q_plot = q_ax.pcolormesh(X,Y,q, cmap='jet', shading='gouraud',vmin=(np.nanmin(q)), vmax=(np.nanmax(q)))
            q_bar = q_fig.colorbar(q_plot)
            q_bar.set_label(var_units[n],fontsize=12)
            q_ax.set_title(var_names[n],fontsize=14)
            q_ax.set_xlabel('X [cm]',fontsize=12)
            q_ax.set_ylabel('Y [cm]',fontsize=12)
            filename = 'plots/'+var_names[n]+'/'+var_names[n][0]+(''.join(filter(str.isdigit,str(dinfo[fi]))))+'.jpg'
            q_fig.savefig(os.path.join(out_path,filename))

#Plots normalized spatial distribution of average values for each variable and fluctuation values at each time step for each variable
    #shock.csv must exist (find_shock must already be run)
    #Lengths (Xlen, Ylen) must be real lengths in cm
    #Cell width (cell_width) must be in cm
    #Variable names (var_names) must be list of strings containing identifiable names for each variable (eg. Temperature)
    #Variable units (var_units) must be list of strings that will be placed on colorbar
    #Average file paths (avg_files) must be list of full file paths to files containing average values, in order in which variable names are given
    #Fluctuation file paths (fluc_files) must be list of full paths to directories containing fluctuation files, in order in which variable names are given
    #Output path (out_path) must be full path to directory in which plots are stored
def proc_plots(Xgrid,Ygrid,Xlen,Ylen,cell_width,var_names,var_units,avg_files,fluc_files,out_path):
#example: proc_plots(4096,4096,16,4,1.45,['Temperature','Heat Release'],[r'$\frac{T}{T_{ref}}$',r'$\frac{Q}{Q_{ref}}$'],['/localdata/scratch/adhav/Source_Comp_H2O2Ar_L16/Averages/Tavg.csv','/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/Averages/Qavg.csv'],['/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/Fluctuations/Temperature/','/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/Fluctuations/Heat_Release/'],'/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/')

    Xsz = np.arange(Xgrid)
    Ysz = np.arange(Ygrid)
    nvars = len(var_names)

    shock = np.asarray(pd.read_csv(os.path.join(out_path,'shock.csv')))
    (vars_norm,L_norm,t_norm) = normalize(Xgrid,Ygrid,avg_files,cell_width,'',shock,out_path)

    X = ((Xsz/Xgrid)*Xlen)/L_norm
    Y = ((Ysz/Ygrid)*Ylen)/L_norm

    if not os.path.isdir(os.path.join(out_path,'plots')):
        os.mkdir(os.path.join(out_path,'plots'))

    for n in range(0,nvars):

        q = np.asarray(pd.read_csv(avg_files[n],header=None))
        q_1 = q[Ysz,...]
        q_avg = q_1[...,Xsz]

        q_avg = q_avg/vars_norm[n]

        qav_fig,qav_ax = plt.subplots()
        qav_plot = qav_ax.pcolormesh(X,Y,q_avg, cmap='jet', shading='gouraud',vmin=(np.nanmin(q_avg)), vmax=(np.nanmax(q_avg)))
        qav_bar = qav_fig.colorbar(qav_plot)
        qav_bar.set_label(var_units[n],fontsize=12)
        qav_ax.set_title('Average '+var_names[n],fontsize=14)
        qav_ax.set_xlabel(r'$\frac{X}{\lambda}$',fontsize=12)
        qav_ax.set_ylabel(r'$\frac{Y}{\lambda}$',fontsize=12)
        filename = 'plots/Average_'+var_names[n]+'.jpg'
        qav_fig.savefig(os.path.join(out_path,filename))


        dinfo = [f for f in os.listdir(fluc_files[n]) if os.path.isfile(os.path.join(fluc_files[n],f))]
        nfiles = len(dinfo)
        for m in range(0,nfiles):
            q = np.asarray(pd.read_csv(os.path.join(fluc_files[n],dinfo[m])))
            q_1 = q[Ysz,...]
            q_fluc = q_1[...,Xsz]

            q_fluc = q_fluc/vars_norm[n]

            if not os.path.isdir(os.path.join(out_path,'plots',var_names[n])):
                os.mkdir(os.path.join(out_path,'plots',var_names[n]))

            qfl_fig,qfl_ax = plt.subplots()
            qfl_plot = qfl_ax.pcolormesh(X,Y,q_fluc,norm=colors.TwoSlopeNorm(vmin=(np.nanmin(q_fluc)),vcenter=0,vmax=(np.nanmax(q_fluc))), cmap='bwr', shading='gouraud')
            qfl_bar = qfl_fig.colorbar(qfl_plot)
            qfl_bar.set_label(var_units[n],fontsize=12)
            qfl_ax.set_title(var_names[n]+' Fluctuation',fontsize=14)
            qfl_ax.set_xlabel(r'$\frac{X}{\lambda}$',fontsize=12)
            qfl_ax.set_ylabel(r'$\frac{Y}{\lambda}$',fontsize=12)
            filename = 'plots/'+var_names[n]+'/Fluctuation_'+var_names[n][0]+(''.join(filter(str.isdigit,str(dinfo[m]))))+'.jpg'
            qfl_fig.savefig(os.path.join(out_path,filename))    

#creates .csv file of average spatial temperature distribution and shock location (x-coordinate at each y-coordinate)
    #temperature directory (t_dir) and output path (out_path) inputs must be full file paths to directory containing extracted temperature data files and directory in which average temperature and shock files are saved, respectively
    #only works for temperature
    #returns array representing shock x position at each y coordinate as output
    #RUN FIRST before other functions involving normalization
def find_shock(Xgrid,Ygrid,t_dir,out_path):
    #example: shock = find_shock(4096,4096,'/localdata/scratch/Haripriya/Detonation_DNS/2_H2O2Ar_L16/Temperature/','/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/')

    Xsz = np.arange(Xgrid)
    Ysz = np.arange(Ygrid)

    dinfo =  [f for f in os.listdir(t_dir) if os.path.isfile(os.path.join(t_dir,f))]

    t_avg = 0
    for fi in range(0,nfiles):

        t_raw = np.asarray(pd.read_csv(os.path.join(t_dir,dinfo[fi]),header=None))
        t_raw2 = t_raw[Ysz,...]
        t_r = t_raw2[...,Xsz]

        t_avg = t_avg+(t_r/nfiles)

    filename = 'Tavg.csv'
    np.savetxt(os.path.join(out_path,filename),t_avg,delimiter=',')    

    shock = np.zeros(Ygrid)
    for i in range(0,(Ygrid-1)):
        for j in range(0,(Xgrid-1)):
            if t_avg[i,j] > 400:
                shock[i] = j
                break

    filename = 'shock.csv'
    np.savetxt(os.path.join(out_path,filename),shock,delimiter=',')

    return shock

#calculates values used for non-dimensionalization of each variable
    #average files (avg_files) must be list of full file paths to average distribution files
    #length (L) and time (t) normalization values can either be set or left empty (''). If left empty, will default to 1
    #shock (shock) must be array of shock x-locations at each y coordinate, run find_shock first
    #output path (out_path) must be full file path to directory used for output files
    #returns list of normalization values for each variable for which an average file is provided (vars_norm), length normalization (L_norm), and time normalization (t_norm) values
def normalize(Xgrid,Ygrid,avg_files,L,t,shock,out_path):
    #example: (vars_norm,L_norm,t_norm) = normalize(4096,4096,['/localdata/scratch/Haripriya/Detonation_DNS/2_H2O2Ar_L16/Averages/Temperature/Tavg.csv','/localdata/scratch/Haripriya/Detonation_DNS/2_H2O2Ar_L16/Averages/Pressure/Pavg.csv','/localdata/scratch/Haripriya/Detonation_DNS/2_H2O2Ar_L16/Averages/CHrelease/Qavg_old.csv'],1.45,'','/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/shock.csv','/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/')

    nvars = len(avg_files)
    vars_norm = []

    for n in range(0,nvars):
        q_avg = np.asarray(pd.read_csv(os.path.join(out_path,avg_files[n])))
        q_mean = np.zeros(Ygrid)
        for i in range(0,(Ygrid-1)):
            xinds = np.arange(shock[i]-1,shock[i]+1).astype('int')
            qi = q_avg[i,xinds]
            q_mean[i] = np.mean(qi)
        q_norm = np.mean(q_mean)
        vars_norm.append(q_norm)

    if not L == '':
        L_norm = L
    else:
        L_norm = 1

    if not t == '':
        t_norm = t
    else:
        t_norm = 1

    return (vars_norm,L_norm,t_norm)


#creates .csv of location of edge of heat release range (x-coordinate for each y-coordinate) for a given threshold
    #average heat release file input (avgq_file)  must be full file path
    #only works for heat release
    #Threshold (thresh) is value of second-derivative of heat release at which the edge is determined
    #kernel size (ksz) determines kernel used for smoothing
    #output path (out_path) must be full file path to directory in which .csv is created
    #RUN BEFORE ANY VOLUME INTEGRATION
    #returns array representing end of heat release generation region
def find_qrange(Xgrid,Ygrid,avgq_file,thresh,ksz,out_path):
    #example: qrange = find_qrange(4096,4096,'/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/Averages/Qavg.csv',10**(-3),100,'/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/')

    Xsz = np.arange(Xgrid)
    Ysz = np.arange(Ygrid)

    Q_raw = np.asarray(pd.read_csv(avgq_file,header=None))
    Q1_raw = Q_raw[Ysz,...]
    Q = Q1_raw[...,Xsz]

    shock = np.asarray(pd.read_csv(os.path.join(out_path,'shock.csv')))
    (Q_norm,L_norm,t_norm) = normalize(Xgrid,Ygrid,avgq_file,'','',shock,out_path)

    Q = Q/Q_norm
    Q_mean = np.mean(Q,axis=0)
    gradq = np.gradient(Q_mean)
    k = np.ones(ksz)/ksz
    gradq_conv = np.convolve(gradq,k,mode='same')
    grad2q = np.gradient(gradq_conv)
    qmax = np.nanargmax(Q_mean)

    for j in range(qmax,(Xgrid-1)):
        if (abs(grad2q[j])/(np.nanmax(abs(grad2q[qmax:Xgrid])))) < thresh:
            qrange = j
            break
    
    qrange = np.repeat(qrange,Ygrid)

    filename = 'qrange.csv'
    np.savetxt(os.path.join(out_path,filename),qrange,delimiter=',')

    return qrange
            
#creates plot containing volume-integrated source terms at each time step and plots of the discrete Fourier transforms of each term over the time range
    #shock and qrange .csv files must exist (shock and find_qrange must have been run)
    #term directories (term_dirs) must be a list of full file paths to directories containing source term .csv files (source must have been run)
    #term names (term_names) must be a list of strings containing names of each source term (eg. ['Tfluc*Qfluc','Tfluc*Qavg'])
    #output directory (out_path) must be full path to directory containing shock and qrange files and in which source term fft plots and final plot will be saved
def integrated_plot(Xgrid,Ygrid,term_dirs,term_names,out_path):
    #example: integrated_plot(4096,4096,['/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/TflucQfluc/','/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/TflucQavg/'],['Tfluc*Qfluc','Tfluc*Qavg'],'/localdata/scratch/adhav/Source_Comp_2_H2O2Ar_L16/')

    shock = np.asarray(pd.read_csv(os.path.join(out_path,'shock.csv')))
    qrange = np.asarray(pd.read_csv(os.path.join(out_path,'qrange.csv')))

    r_fig,r_ax = plt.subplots()

    nterms = len(term_dirs)
    for n in range(0,nterms):
        dinfo = [f for f in os.listdir(term_dirs[n]) if os.path.isfile(os.path.join(term_dirs[n],f))]
        nfiles = len(dinfo)

        term_ind = term_names[n].split('*')
        nterms_ind = len(term_ind)

        r = []
        times = []
        r = np.array(r)
        times = np.array(times)

        #volume integration of Rayleigh source term using trapezoidal rule, appended each iteration to array
        for fi in range(0,nfiles):
        
            term_raw = np.asarray(pd.read_csv(os.path.join(term_dirs[n],dinfo[fi]),header=None))
            term_int = np.zeros(Ygrid)
        
            for i in range(0,(Ygrid-1)):
        
                sz = qrange[i]-shock[i]
                xinds = np.arange(shock[i]-1,qrange[i]-1).astype('int')
                termi = term_raw[i,xinds]
                term_int[i] = np.trapz(termi)

            rt = np.trapz(term_int)
            r = np.append(r,rt)

            t = int(''.join(filter(str.isdigit,str(dinfo[fi]))))
            times = np.append(times,t)

        #sorting
        ind_arr = np.argsort(times,axis=0)
        times = np.take_along_axis(times,ind_arr,axis=0)
        r = np.take_along_axis(r,ind_arr,axis=0)
   
        R = fft(r)
        N = len(R)
        n1 = np.arange(N)
        freq = n1/N
        n1_one = N//2
        freq_one = freq[:n1_one]

        label_name = []
        for m in range(0,nterms_ind):
            if 'avg' in term_ind[m]:
                label_name[m] = r'$\overline{%s}$' % term_ind[m][0]
            else:
                label_name[m] = r"$%s'$" % term_ind[m][0]
        label_full = ''.join(label_name)

        f_fig,f_ax = plt.subplots()
        f_ax.stem(freq_one,np.abs(R[:n1_one]),'b',markerfmt='',basefmt='-b')
        f_ax.set_title('Discrete Fourier Transform',fontsize=14)
        f_ax.set_xlabel('Frequency',fontsize=12)
        f_ax.set_ylabel(label_full+' Amplitude',fontsize=12)
        filename = term_name+'_integrated_fft.jpg'
        f_fig.savefig(os.path.join(out_path,filename))

        #plotting both source terms at each time step on same plot
        r_ax.plot(times,r,label=label_full)
        
    r_ax.set_title('Volume-Integrated Heat Release Rayleigh Source',fontsize=14)
    r_ax.set_xlabel('Time Step',fontsize=12)
    r_ax.set_ylabel('Integrated Terms',fontsize=12)
    r_ax.legend(loc='lower right')
    filename = 'integrated_sources.jpg'
    r_fig.savefig(os.path.join(out_path,filename)) 

