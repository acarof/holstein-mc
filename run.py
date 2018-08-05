from utils import *
import os, sys
from itertools import product as iterprod

from input import *



try:
    os.mkdir(path)
except:
    print '%s exists already' % path




mega_list = []
second_list = [sublist[1:] for sublist in list_param]
total_list = list(iterprod(*second_list))
for sublist in total_list:
    subdict = {}
    for index in range(len(sublist)):
        subdict.update({
            list_param[index][0]: sublist[index]
        })
    subdict.update({
        'PATH' : path
    })
    mega_list.append(subdict)

print 'Number of traj:', len(mega_list)

def prerun_mc(dict_):
    mass = dict_.get('MASS', 1)
    coupling = dict_['COUPLING']
    dim = len(coupling)
    coupling_ha = 0.036749305 * np.array(coupling)  # Convert in Ha
    frequency = dict_.get('FREQUENCY', 1600)
    frequency_ha = frequency * 0.00000455633528 #  Convert from cm in Ha
    reorga = dict_.get('REORGA', 0.200)
    reorga_ha = reorga * 0.036749305  # Convert from eV to Ha
    kbT = dict_.get('TEMP', 300) *0.0000031667909 # convert from K to Ha
    nsteps = dict_.get('NSTEPS', 1000000)
    grid = dict_.get('GRID', 0.1)
    length = dict_.get('LENGTH')
    #path = dict_.get('PATH', '.')
    sampling = dict_.get('SAMPLE', ['local'])
    toprint = dict_.get('PROPERTIES', ['energies'])

    massy = dict_.get('MASSY', 0)
    freqy = dict_.get('FREQY', 0)
    beta = dict_.get('BETA', 0)



    properties = run_monte_carlo(mass, coupling_ha, frequency_ha, reorga_ha,
                                 massy, freqy, beta,
                                 kbT, length, nsteps, grid, dim,
                                 sampling=sampling,
                                 toprint=toprint)

    title_set = dict_.get('TITLE')
    title = ''
    for info in sorted(title_set):
        title += '%s_%s-' % (info, dict_[info])
    for prop in toprint:
        write_file(properties[prop], title = '%s/%s-%s.dat' % (path, prop, title))



# RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
try:
    nworker = int(sys.argv[1])
except:
    nworker = 1

print "nworker is %s" % nworker

if nworker == 1:
    print 'Use serial'
    for dict_ in mega_list:
        prerun_mc(dict_)
elif nworker == 0:
    from multiprocessing import Pool, cpu_count
    pool = Pool(cpu_count())
    print 'Use parallel with %s processors' % cpu_count()
    pool.map(prerun_mc, mega_list)
elif nworker > 1:
    from multiprocessing import Pool
    pool = Pool(nworker)
    print 'Use parallel with %s processors' % nworker
    pool.map(prerun_mc, mega_list)
