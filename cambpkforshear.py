import numpy,scipy


z=numpy.loadtxt('zfile.dat')


nz=z.shape[0]



for i in xrange(1,2):
    fileconcat='/home/sdeb/EmuCambComp1.0/CambConcat/cambzcom_%s.dat'%i

    file=open(fileconcat,mode='w')
    for j in xrange(nz,0,-1):
    
        file_matter='Results_redshifts_%s/Matter_%s_matterpower_%s.dat'%(i,i,j)
    
        f=open(file_matter,mode='r')
        L=f.readlines()
        file.writelines(L)
        f.close()

    file.close()
