import os,numpy,sys

om_l=0.120
om_h=0.155
omb_l=0.0215
omb_h=0.0235
ns_l=0.85
ns_h=1.05
s8_l=0.616
s8_h=0.9
w_l=-1.30
w_h=-0.70
num=2
n=1

om_mid=(om_l+om_h)/2
omb=(omb_l+omb_h)/2
ns_mid=(ns_l+ns_h)/2
s8_mid=(s8_l+s8_h)/2
w_mid=(w_l+w_h)/2


f=open('table.txt',mode='r')

L=f.readlines()
om=numpy.array([])
omb=numpy.array([])
ns=numpy.array([])
w=numpy.array([])
s8=numpy.array([])
h=numpy.array([])

for l in L:
    line=l.split()
    om=numpy.append(om,eval(line[2]))
    omb=numpy.append(omb,eval(line[4]))
    ns=numpy.append(ns,eval(line[6]))
    w=numpy.append(w,-eval(line[8]))
    s8=numpy.append(s8,eval(line[10]))
    h=numpy.append(h,eval(line[12]))

n=s8.shape[0]
id=numpy.where(s8<s8_l)
s8[id]=s8_l
n=37

directory="Config"
if not os.path.exists(directory):
        os.makedirs(directory)

directory_new='Params'    
if not os.path.exists(directory_new):
        os.makedirs(directory_new)

for i in xrange(2,37):
    h0=h[i]*100.+4.
    if (h0>84):
        h0=84
    
    
    filename="Config/params_%s.config"%i
    f=open(filename,mode='w')
    param_out="out_params_%s.dat"%i
    f.write("empow_%s.dat \n"%i)
    f.write("%f \n"%omb[i])
    f.write("%f \n"%om[i])
    f.write("%f \n"%ns[i])
    f.write("%f \n"%h0)
    f.write("%f \n"%w[i])
    f.write("%f \n"%s8[i])
    f.write("%i \n"%num)
    f.write(param_out)
    f.close()
    os.system('emu_new.exe %s'%filename)
    os.system('mv kfile.dat kfile_new.dat')'''
    '''f=open(filename,mode='w')
    param_out="out_params_%s_old.dat"%i
    f.write("empow_%s_old.dat \n"%i)
    f.write("%f \n"%om[i])
    f.write("%f \n"%omb[i])
    f.write("%f \n"%ns[i])
    f.write("%f \n"%s8[i])
    f.write("%f \n"%w[i])
    f.write("%i \n"%num)
    f.write(param_out)
    f.close()
    os.system('emu.exe %s'%filename)
    #sys.exit()
    
    #H0=h0*100.
    
    #now change parameters in camb
    diff_s8=100.
    eps=0.001
    const=3.65

    directory='Results_redshifts_%s'%i
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    
    amp=(s8[i]**2)*const*1.e-9
    print amp,s8[i]
    file0=open("params.x",mode='r')
    LINES=file0.readlines()
    filecamb="Params/params_%s.x"%i
    fileout=open(filecamb,mode='w')
    for line in LINES:
        #print line
        if "output_root = " in line:
            fileout.write("output_root = Results_redshifts_%s/Matter_%s \n"%(i,i))
        elif "ombh2          = " in line:
            fileout.write("ombh2          = %s \n"%(omb[i]))
        elif "omch2          = " in line:
            fileout.write("omch2          = %s \n"%(om[i]-omb[i]))                   
        elif "hubble" in line:
            fileout.write("hubble         = %s \n"%h0)                   
        elif "scalar_amp(1)             = " in line:
            fileout.write("scalar_amp(1)             = %s \n"%amp)                   
        elif "scalar_spectral_index(1)  = " in line:
            fileout.write("scalar_spectral_index(1)  =%s \n"%ns[i])                   
        elif "w              = " in line:
            fileout.write("w              = %s \n"%w[i])                   
    
        else:
            fileout.write(line)
    file0.close()
    fileout.close()
    camout="camb_%s.dat"%i
    os.system("camb_new1.0 %s>%s"%(filecamb,camout))
    file=open(camout,mode='r')
    L=file.readlines()
    for line in L:
        if "sigma8" in line:
            l0=numpy.array(line.split())
            n0=l0.shape[0]
            s8cam=eval(l0[n0-1])
           

    diff_s8=(s8[i]-s8cam)
    print "Difference in sigma_8=",diff_s8
    print "Sigma_8_camb  ",s8cam
    if (numpy.abs(diff_s8)>eps):
        amp_new=0.999*(amp/s8cam**2)*s8[i]**2    
    print "Amp_new=  ",amp_new
    file0=open("params_redshifts.x",mode='r')
    LINES=file0.readlines()
    filecamb="Params/params_%s.x"%i
    fileout=open(filecamb,mode='w')
    
    for line in LINES:
        #print line                                                                                                                           
        if "output_root = " in line:
            fileout.write("output_root = Results_redshifts_%s/Matter_%s \n"%(i,i))
        elif "ombh2          = " in line:
            fileout.write("ombh2          = %s \n"%(omb[i]))
        elif "omch2          = " in line:
            fileout.write("omch2          = %s \n"%(om[i]-omb[i]))
        elif "hubble" in line:
            fileout.write("hubble         = %s \n"%h0)
        elif "scalar_amp(1)             = " in line:
            fileout.write("scalar_amp(1)             = %s \n"%amp_new)
        elif "scalar_spectral_index(1)  = " in line:
            fileout.write("scalar_spectral_index(1)  =%s \n"%ns[i])
        elif "w              = -1" in line:
            fileout.write("w              = %s \n"%w[i])
            
        else:
            fileout.write(line)
    file0.close()
    fileout.close()
    camout="camb_%s.dat"%i
    os.system("camb_new1.0 %s>%s"%(filecamb,camout))
    file=open(camout,mode='r')
    L=file.readlines()
    for line in L:
        if "sigma8" in line:
            l0=numpy.array(line.split())
            n0=l0.shape[0]
            s8cam=eval(l0[n0-1])



    print "Sigma_8 Emu= ",s8[i], "Sigma_8 Camb= ",s8cam
    

