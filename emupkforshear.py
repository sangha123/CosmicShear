import numpy,scipy,sys

k=numpy.loadtxt('kfile_new.dat')
z=numpy.loadtxt('zfile.dat')

for m in xrange(1,13):

    filename='empow_%s'%m+'.dat'
    l=filename[6:]

    Pkemu=numpy.loadtxt(filename)

    fileout='/home/sdeb/EmuCambComp1.0/EmuOut/FEmu'+l

    f=open(fileout,mode='w')
    
    nk=k.shape[0]
    nz=z.shape[0]

    for i in xrange(nz):
        for j in xrange(nk):
            f.write('%f  %f\n'%(k[j],Pkemu[j+i*nk]))

    f.close()
