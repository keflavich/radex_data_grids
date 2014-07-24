#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
# #! /usr/bin/python
#
import math
import os
import time
import sys
import re

# Sometimes, fortran outputs things like "1.404+106" instead of "1.404E+106"
bad_exp = re.compile("([0-9])\+")

# Run a series of Radex models to estimate temperature & density
# from observed ratios of H2CO 1-1/2-2, 1-1/3-3, and 2-2/3-3 lines
#
# Original comes from the RADEX team:
# https://www.sron.rug.nl/~vdtak/radex/index.shtml#script
#
# This version has been heavily modified by Adam Ginsburg
# (adam.ginsburg@colorado.edu, or keflavich@gmail.com)
# in May 2010
#
# Directions:
# To call this code, run:
# mpirun -np 8 ./radex_grid.py > mpi_radex_grid.log &
# In the above statement, "-np 8" means "use 8 processors".  Output
# is redirected to a log file (recommended - the logging is verbose).
#
# The outputs will include .dat files with names specified by the "acts" list
# and "suffix" below, and columns Temperature, log10(dens), log10(col),
# Tex_low, Tex_hi, TauLow, TauUpp, TrotLow, TrotUpp, FluxLow, FluxUpp
# Additionally, it will create a radex.out file.
#
# The code creates (# processors) subdirectories, reformats that data, and
# removes the temporary subdirectories.

# Grid boundaries
#
tmin = 5.0   # minimum kinetic temperature (K)
tmax = 300.0  # maximum kinetic temperature (K)
nmin = 1e3   # minimum H2 density (cm^-3)
nmax = 1e7   # maximum H2 density (cm^-3)
cmin = 1e11  # minimum molecular column density
cmax = 1e17  # maximum molecular column density
#otopmin = 1e-3 # minimum ortho-to-para ratio
#otopmax = 3.0  # maximum ortho-to-para ratio

ntemp = 60     # number of temperature points
ndens = 21    # number of density points
ncol  = 51    # number of column points
#notop = 7      # number of ortho-to-para ratio points

# user does not need to modify these formulae
# they are equivalent to temperatures = numpy.linspace(tmin,tmax,ntemp).tolist()
temperatures = [ tmin + (ii) / float(ntemp-1) * (tmax-tmin)  for ii in range(ntemp) ]
densities    = [ 10**( math.log10(nmin) + (ii) / float(ndens-1) * (math.log10(nmax)-math.log10(nmin)) )  for ii in range(ndens) ]  # LOG
columns      = [ 10**( math.log10(cmin) + (ii) / float(ncol-1) * (math.log10(cmax)-math.log10(cmin)) )  for ii in range(ncol) ]  # LOG
#orthopara    = [ 1e-3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 ]
# LINEAR densities    = [ nmin + (ii) / float(ndens-1) * (nmax-nmin)  for ii in range(ndens) ] # LINEAR
# LINEAR columns      = [ cmin + (ii) / float(ncol-1) * (cmax-cmin)  for ii in range(ncol) ] # LINEAR

#
# Parameters to keep constant
#
tbg   = 2.73         # background radiation temperature
cdmol_default = 1e12 # low enough to be optically thin
dv    = 5.0          # line width (km/s)

#mole = 'o-h2co_troscompt'  # molecular data name
mole = 'ph2co-h2' # Faure version

#import sys
#if len(sys.argv) > 1:
#    orthopararatio = sys.argv[1] 
#else:
#orthopararatio = 1e-3 # H2 ortho-to-para ratio

# Naming suffix to append to line name defined in "acts" below
#suffix = "_T=5to55_lvg_troscompt_100square_orthopara0"
suffix = "_temperature_para_300Kmax_5kms"
# acts = list of lines to ratio.  Filenames must include ".dat" in order for suffix to be applied
acts = ([218.22,218.47,'303-202_322-221.dat'],
        [218.47,218.76,'322-221_321-220.dat'],
        [218.22,218.76,'303-202_321-220.dat'])
flow = 4.0     # lowest frequency transition to include in output file
fupp = 250.0   # highest frequency transition to include in output file
bw    = 0.01   # "bandwidth": free spectral range around line (used to say which lines get printer)

# can run sphere or lvg
# (note that you must generate these executables and name them yourself,
# and they must be in your path or you can specify the full path)
executable = "radex_lvg"
# executable = "radex_sphere"

# verbosity
# 2 = output 1 line for every RADEX run (redirect to log file!)
# 1 = just output major statements (OK to print to screen)
# 0 = silent (only prints exceptions)
verbose = 2

#
# No user changes needed below this point.
#
def write_input(infile,tkin,nh2,cdmol=cdmol_default):
    """
    Write radex.inp file parameters
    """
    infile.write(mole+'.dat\n')
    infile.write('radex.out\n')
    infile.write(str(flow*(1-bw))+' '+str(fupp/(1-bw))+'\n')
    infile.write(str(tkin)+'\n')
    #infile.write('2\n')
    #infile.write('o-H2\n')
    #infile.write(str(nh2/(float(orthopararatio)+1.0)*float(orthopararatio))+'\n')
    #infile.write('p-H2\n')
    #infile.write(str(nh2/(float(orthopararatio)+1.0))+'\n')
    infile.write('1\n')
    infile.write('H2\n')
    infile.write(str(nh2)+'\n')
    infile.write(str(tbg)+'\n')
    infile.write(str(cdmol)+'\n')
    infile.write(str(dv)+'\n')

def read_radex(file,flow,fupp,bw=bw):
    """ 
    less hack-ey way to read radex.out files
    """
    line = file.readline()
    if line == '':
        return 0
    words = line.split()
    freq = 0
    while len(words) > 1 and words[1] == '--': # this case should never have to happen....
        raw_line = file.readline()
        line = bad_exp.sub("\\1E+",raw_line)
        words = line.split()
    while not(freq*(1-bw) < flow < freq*(1+bw)):
        if words[1] == 'T(kin)':
            tkin = float(words[3])
        elif line.find("Density of H2") != -1:
            dens = float(words[5])
        elif line.find("Density of pH2") != -1:
            pdens = float(words[5])
        elif line.find("Density of oH2") != -1:
            odens = float(words[5])
        elif line.find("Column density") != -1:
            col = float(words[4])
        raw_line = file.readline()
        line = bad_exp.sub("\\1E+",raw_line)
        words = line.split()
        if words[1] == '--':
            freq = float(words[4])
    TexLow   = float(words[6])
    TauLow   = float(words[7])
    TrotLow  = float(words[8])
    FluxLow  = float(words[11])
    line = file.readline()
    words = line.split()
    while  not(freq*(1-bw) < flow < freq*(1+bw)):
        line = file.readline()
        words = line.split()
        if words[1] == '--':
            freq = float(words[4])
    TexUpp   = float(words[6])
    TauUpp   = float(words[7])
    TrotUpp  = float(words[8])
    FluxUpp  = float(words[11])
    while len(words) > 1 and words[1] == '--':
        line = file.readline()
        words = line.split()
    return tkin,dens,col,TexLow,TexUpp,TauLow,TauUpp,TrotLow,TrotUpp,FluxLow,FluxUpp
 
# Begin main program

start = time.time()

# Allow for parallel running.  If mpirun is not used, will operate in
# single-processor mode
try:
    from mpi4py import MPI
    mpirank = MPI.COMM_WORLD.rank  # processor ID
    mpisize = MPI.COMM_WORLD.size  # number of processors
except ImportError:
    print "mpi4py not found.  Using a single processor."
    mpirank = 0
    mpisize = 1
if mpisize > 1:
    # each processor gets 1/n_processors of the temperatures, in order
    # If you want to run in parallel with just 1 temperature, 
    # these lines need to be changed
    splits = [ int( math.floor(ii / float(mpisize) * ntemp) ) for ii in range(mpisize+1) ] 
    temperatures = temperatures[splits[mpirank]:splits[mpirank+1]]

    # Make a separate subdirectory for each temperature
    # ("temp" means temporary, though)
    newdir = "radex_temp_%02i" % mpirank
    pwd = os.getcwd() # will return to PWD later
    try:
        os.mkdir(newdir)
    except OSError:
        print "%s exists, continuing" % newdir
    os.chdir(newdir)

if verbose > 0: print "Running code ",executable," with temperatures ",temperatures," densities ",densities," and columns ",columns

for iact,act in enumerate(acts):
    lowfreq = act[0]
    uppfreq = act[1]
    gfil = act[2].replace(".dat",suffix+".dat")
    
    infile = open('radex.inp','w')
    if verbose > 0: print "Processor %i: Starting " % mpirank,gfil

    for temp in temperatures:
        for dens in densities:
            for col in columns:

                write_input(infile,temp,dens,col)
                if (temp == temperatures[-1] and dens == densities[-1] and col == columns[-1]):
                    infile.write('0\n')
                    infile.close()
                else:
                    infile.write('1\n')

                # DEBUG logging
                if verbose > 1: print "Processor %i " % mpirank,
                if verbose > 1: print "temp : %g" % (temp),
                if verbose > 1: print "Column : %g" % (col),
                if verbose > 1: print "dens : %g" % (dens)

    if verbose > 0: print "Processor %i: Finished writing infiles." % mpirank
    if iact == 0:
        if verbose > 0: print "Processor %i: Starting radex code." % mpirank
        #status = os.system('%s < radex.inp > /dev/null' % executable)
        status = 0
        if status != 0:
            print "Command %s failed with exit status %i" % ('%s < radex.inp > /dev/null' % executable,status)
            import pdb; pdb.set_trace()
        
        if verbose > 0: print "Processor %i: Finished Radex." % mpirank

    if verbose > 0: print "Processor %i: Beginning output parsing." % mpirank
    if verbose > 1: print "Processor %i: Printing to file %s." % (mpirank,gfil)
    grid = open(gfil,'w')
    fmt  = '%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n'
    grid.write(fmt.replace('.3e','s') % ("Temperature","log10(dens)",
        "log10(col)","Tex_low","Tex_hi","TauLow","TauUpp","TrotLow","TrotUpp","FluxLow","FluxUpp"))

    radexfile  = open('radex.out')

    rmin = 100
    rmax = 0.1

    while(1):
        radex_out = read_radex(radexfile,lowfreq,uppfreq)
        if radex_out == 0:
            break
        else:
            temp,dens,col,tlow,tupp,taulow,tauupp,trotlow,trotupp,fluxlow,fluxupp = radex_out

        grid.write(fmt %(temp, math.log10(dens), math.log10(col),
            tlow, tupp, taulow, tauupp, trotlow,trotupp,fluxlow,fluxupp))

    grid.close()
    radexfile.close()
    if verbose > 1: print "Processor %i: Completed output parsing.  Wrote %i temperatures, %i densities, %i columns." % \
            (mpirank,len(temperatures),len(densities),len(columns))

stop = time.time()
dure = stop - start
if verbose > 0: print "Processor %i Run time = %f seconds" % (mpirank,dure)
if mpisize > 1:
    os.chdir(pwd)

print """
python ../plot_grids.py 303-202_321-220_temperature_para_300Kmax_5kms.dat --script
python ../plot_grids.py 303-202_322-221_temperature_para_300Kmax_5kms.dat --script
python ../plot_grids.py 322-221_321-220_temperature_para_300Kmax_5kms.dat --script
"""

if mpisize > 1 and mpirank == 0:
    MPI.COMM_WORLD.Barrier()
    if verbose > 0: print "Processor %i: Starting cleanup" % mpirank
    import glob
    filelist = glob.glob("radex_temp_00/*.dat")
    for file in filelist:
        status = os.system("cp %s %s" % (file,file.replace("radex_temp_00/","") ) )
        if status != 0:
            print "Processor %i: " % mpirank,"Command ",("cp %s %s" % (file,file.replace("radex_temp_00/","") ) )," failed with status ",status
        for ii in xrange(1,mpisize):
            if verbose > 1: 
                wc = os.popen('wc %s' % file.replace("_00","_%02i" % ii)).readline().split()[0]
                print "Processor %i merging file %s with length %s onto %s" % (mpirank, 
                        file.replace("_00","_%02i" % ii),wc,file.replace("radex_temp_00/",""))
            status = os.system("tail -n +2 %s >> %s" % (file.replace("_00","_%02i" % ii),file.replace("radex_temp_00/","") ) )
            if status != 0:
                print "Processor %i: " % mpirank,"Command ",("tail -n +2 %s >> %s" % (file.replace("_00","_%02i" % ii),file.replace("radex_temp_00/","") ) )," failed with status ",status
    radexoutlist = glob.glob("radex_temp_*/radex.out")
    try:
        # if a radex.out file exists, move it to radex.out.old
        os.system("mv radex.out radex.out.old")
        os.system("touch radex.out")
    except OSError:
        pass
    for file in radexoutlist:
        os.system("cat %s >> radex.out" % file)
    os.system("rm -r radex_temp_*")
    if verbose > 0: print "Processor %i: " % mpirank,"Cleanup completed"
