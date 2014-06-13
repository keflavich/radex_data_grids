#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
""" WARNING - THIS IS NOT A CUBE GRID BUT A LINEAR GRID!!"""
# #! /usr/bin/python
#
import math
import os
import time
#
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
#tmin = 5.0   # minimum kinetic temperature (K)
#tmax = 55.0  # maximum kinetic temperature (K)
nmin = 1e-2   # minimum H2 density (cm^-3)
nmax = 1e10   # maximum H2 density (cm^-3)
#cmin = 1e11  # minimum molecular column density
#cmax = 1e17  # maximum molecular column density

orthopararatio = [0.0,0.1,1.0,3.0]
ntemp = 1      # number of temperature points
ndens = 1501   # number of density points
ncol  = 1      # number of column points
abundance = 10**-8.0,10**-8.25,10**-8.5,10**-8.75,10**-9,10**-9.5,10**-10.0
# at 1 km/s linewidth, 1 km/s / pc, we have column = 3.08e18 * density * abundance 
# abund = column - log10(dvdr * pc) - density
# column = abund + log10(dvdr * pc) + density 

# user does not need to modify these formulae
# they are equivalent to temperatures = numpy.linspace(tmin,tmax,ntemp).tolist()
temperatures = [20] #[ tmin + (ii) / float(ntemp-1) * (tmax-tmin)  for ii in range(ntemp) ]
densities    = [ 10**( math.log10(nmin) + (ii) / float(ndens-1) * (math.log10(nmax)-math.log10(nmin)) )  for ii in range(ndens) ]  # LOG
columns      = [ 3.08e18 * a * d for a in abundance for d in densities ]
densities = densities * len(abundance)
# LINEAR densities    = [ nmin + (ii) / float(ndens-1) * (nmax-nmin)  for ii in range(ndens) ] # LINEAR
# LINEAR columns      = [ cmin + (ii) / float(ncol-1) * (cmax-cmin)  for ii in range(ncol) ] # LINEAR

#
# Parameters to keep constant
#
tbg   = 2.73         # background radiation temperature
cdmol_default = 1e12 # low enough to be optically thin (this is just the default passed to write_input; it is overridden)
dv    = 1.0          # line width (km/s)

mole = 'o-h2co_troscompt'  # molecular data name

# Naming suffix to append to line name defined in "acts" below
suffix = "_XH2CO=1e-9_troscompt"
# acts = list of lines to ratio.  Filenames must include ".dat" in order for suffix to be applied
acts = ([4.8,14.5,'1-1_2-2.dat'],[4.8,29.0,'1-1_3-3.dat'],[14.5,29.0,'2-2_3-3.dat'])
flow = 4.0     # lowest frequency transition to include in output file
fupp = 200.0   # highest frequency transition to include in output file
bw    = 0.01   # "bandwidth": free spectral range around line (used to say which lines get printer)

# can run sphere or lvg
# (note that you must generate these executables and name them yourself,
# and they must be in your path or you can specify the full path)
executable = "/Users/adam/repos/radex/bin/radex_lvg"
# executable = "radex_sphere"

# verbosity
# 2 = output 1 line for every RADEX run (redirect to log file!)
# 1 = just output major statements (OK to print to screen)
# 0 = silent (only prints exceptions)
verbose = 2

#
# No user changes needed below this point.
#
def write_input(infile,tkin,nh2,orthopararatio=orthopararatio,cdmol=cdmol_default):
    """
    Write radex.inp file parameters
    """
    infile.write(mole+'.dat\n')
    infile.write('radex.out\n')
    infile.write(str(flow*(1-bw))+' '+str(fupp/(1-bw))+'\n')
    infile.write(str(tkin)+'\n')
    northo = (nh2/(orthopararatio+1.0)*orthopararatio)
    if orthopararatio <= 1e-3 or northo < 1e-3:
        infile.write('1\n')
        infile.write('p-H2\n')
        infile.write(str(nh2/(orthopararatio+1.0))+'\n')
    else:
        infile.write('2\n')
        infile.write('o-H2\n')
        infile.write(str(northo)+'\n')
        infile.write('p-H2\n')
        infile.write(str(nh2/(orthopararatio+1.0))+'\n')
    infile.write(str(tbg)+'\n')
    infile.write(str(cdmol)+'\n')
    infile.write(str(dv)+'\n')

def read_radex(outfile,flow,fupp):
    """
    A hack-ey means of reading a radex.out file.  
    Cycles through it based on fixed format
    """
    line  = outfile.readline()
    words = line.split()
    while (words[-1] != "FLUX"):
        line  = outfile.readline()
        words = line.split()
        if words[1] == "T(kin)":
            temp  = float(words[-1])
        if words[1] == "Density" and words[3] == "H2":
            dens  = float(words[-1])
        if words[1] == "Density" and words[3] == "pH2":
            pdens  = float(words[-1])
        if words[1] == "Density" and words[3] == "oH2":
            odens  = float(words[-1])
        if words[1] == "Column":
            col  = float(words[-1])
    line  = outfile.readline()
    line  = outfile.readline()
    words = line.split()
    ftmp  = float(words[4])
    while ((ftmp < flow*(1-bw)) or (ftmp > flow/(1-bw))):
        line  = outfile.readline()
        words = line.split()
        ftmp  = float(words[4])
    low   = float(words[-2])
    TexLow   = float(words[6])
    TauLow   = float(words[7])
    TrotLow  = float(words[8])
    FluxLow  = float(words[11])
    line  = outfile.readline()
    words = line.split()
    ftmp  = float(words[4])
    while ((ftmp < fupp*(1-bw)) or (ftmp > fupp/(1-bw))):
        line  = outfile.readline()
        words = line.split()
        ftmp  = float(words[4])
    upp   = float(words[-2])
    TexUpp   = float(words[6])
    TauUpp   = float(words[7])
    TrotUpp  = float(words[8])
    FluxUpp  = float(words[11])
    if 'odens' in locals():
        orthopararatio = odens/pdens
    else:
        orthopararatio = 0.0
    return temp,dens,col,TexLow,TexUpp,TauLow,TauUpp,TrotLow,TrotUpp,FluxLow,FluxUpp,orthopararatio
 
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
    # each processor gets 1/n_processors of the densities AND columns, in order
    # If you want to run in parallel with just 1 temperature, 
    # these lines need to be changed
    splits = [ int( math.floor(ii / float(mpisize) * ndens) ) for ii in range(mpisize+1) ] 
    densities = densities[splits[mpirank]:splits[mpirank+1]]
    columns = columns[splits[mpirank]:splits[mpirank+1]]

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

    temp = temperatures[0]
    for opr in orthopararatio:
        for dens,col in zip(densities,columns):

            write_input(infile,temp,dens,opr,col)
            if (temp == temperatures[-1] and dens == densities[-1] and col == columns[-1] and opr == orthopararatio[-1]):
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
        status = os.system('%s < radex.inp > /dev/null' % executable)
        if status != 0:
            print "Command %s failed with exit status %i" % ('%s < radex.inp > /dev/null' % executable,status)
            import pdb; pdb.set_trace()
        if verbose > 0: print "Processor %i: Finished Radex." % mpirank

    if verbose > 0: print "Processor %i: Beginning output parsing." % mpirank
    grid = open(gfil,'w')
    fmt  = '%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n'
    grid.write(fmt.replace('.3e','s') % ("Temperature","log10(dens)",
        "log10(col)","opr","Tex_low","Tex_hi","TauLow","TauUpp","TrotLow","TrotUpp","FluxLow","FluxUpp"))

    outfile  = open('radex.out')

    rmin = 100
    rmax = 0.1

    for opr in orthopararatio:
        for dens,col in zip(densities,columns):

            temp,dens,col,tlow,tupp,taulow,tauupp,trotlow,trotupp,fluxlow,fluxupp,opr = read_radex(outfile,lowfreq,uppfreq)

            grid.write(fmt %(temp, math.log10(dens), math.log10(col), opr,
                tlow, tupp, taulow, tauupp, trotlow,trotupp,fluxlow,fluxupp))

    grid.close()
    outfile.close()

stop = time.time()
dure = stop - start
if verbose > 0: print "Processor %i Run time = %f seconds" % (mpirank,dure)
if mpisize > 1:
    os.chdir(pwd)

    MPI.COMM_WORLD.Barrier()

    if mpisize > 1 and mpirank == 0:
        if verbose > 0: print "Starting cleanup"
        import glob
        filelist = glob.glob("radex_temp_00/*.dat")
        for file in filelist:
            os.system("cp %s %s" % (file,file.replace("radex_temp_00/","") ) )
            for ii in xrange(1,mpisize):
                os.system("tail -n +2 %s >> %s" % (file.replace("_00","_%02i" % ii),file.replace("radex_temp_00/","") ) )
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
        if verbose > 0: print "Cleanup completed"
