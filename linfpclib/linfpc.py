# linfpcwrite.py>

import numpy as np
import matplotlib.pyplot as plt
import os
import math

#workflow
#map form params and maps namelist
#pick root
#make fpc
#sweep `out` from root
#verify if this is the root we want

def find_nearest(array, value): #random but very useful function
    """
    Finds index of element in array with value closest to given value

    Paramters
    ---------
    array : 1d array
        array
    value : float
        value you want to approximately find in array

    Returns
    -------
    idx : int
        index of nearest element
    """

    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

class plume_input:
    #class that have dictss with all inputs for each namelist
    #params
    #fpc
    #species_n
    #maps
    #scan_input_n
    #guess_n
    def __init__(self,dataname):
        self.dataname = dataname

    def read_input(): #typically load from sample
        pass

    def set_params(self,betap,kperp,kpar,vtp,nspec,nscan,option,\
                    nroot_max,use_map,writeOut):
        use_map_out = '.true.'
        if(not(use_map)):
            use_map_out = '.false.'

        usewriteOut_out = '.true.'
        if(not(use_map)):
            usewriteOut_out = '.false.'

        self.params = {'betap':betap,
                       'kperp':kperp,
                       'kpar':kpar,
                       'vtp':vtp,
                       'nspec':int(nspec),
                       'nscan':int(nscan),
                       'option':int(option),
                       'nroot_max':int(nroot_max),
                       'use_map':use_map_out,
                       'writeOut':usewriteOut_out
        }

    def set_fpc(self,vperpmin=0,vperpmax=0,vparmin=0,vparmax=0,delv=0,vxmin=None,vxmax=None,vymin=None,vymax=None,vzmin=None,vzmax=None,elecdircontribution=0.):
        self.fpc = {'vperpmin':vperpmin,
                    'vperpmax':vperpmax,
                    'vparmin':vparmin,
                    'vparmax':vparmax,
                    'delv':delv,
                    'elecdircontribution':elecdircontribution
        }
        if(vxmin!=None and vxmax!=None and vymin!=None and vymax!=None and vzmin!=None and vzmax!=None):
           self.fpc['vxmin'] = vxmin
           self.fpc['vxmax'] = vxmax
           self.fpc['vymin'] = vymin 
           self.fpc['vymax'] = vymax 
           self.fpc['vzmin'] = vzmin 
           self.fpc['vzmax'] = vzmax  

    def set_maps(self,loggridw,omi,omf,gami,gamf,positive_roots):
        loggridw_out = '.true.'
        if(not(loggridw)):
            loggridw_out = '.false.'

        positive_roots_out = '.true.'
        if(not(loggridw)):
            positive_roots_out = '.false.'

        self.maps = {'loggridw':loggridw_out,
                     'omi':omi,
                     'omf':omf,
                     'gami':gami,
                     'gamf':gamf,
                     'positive_roots':positive_roots_out
        }

    def make_species(self,tauS, muS, alphS, Qs, Ds, vvS,spec_n=-1):
        tempspecies = {'tauS':tauS,
                       'muS':muS,
                       'alphS':alphS,
                       'Qs':Qs,
                       'Ds':Ds,
                       'vvS':vvS
        }
        if(len(self.species[0]) == 0):
            if(spec_n == -1 or spec_n == 1):
                print("No species found, creating first species...")
                self.species = [tempspecies]
            else:
                print("Warning: spec_n ==",spec_n,"but there are no species here...")
        else:
            num_spec = len(self.species)
            if(spec_n == -1 or num_spec+1 == spec_n):
                print("Appending species to list. Total species is now ",num_spec+1)
                self.species.append(tempspecies)
            elif(num_spec+1 < spec_n):
                print("Warning: there are only ",num_spec,"species but user requested we create spec number",spec_n)
                print("Please input a valid spec_n...")
            else:
                print("Replacing species number",spec_n)
                self.species[spec_n-1] = tempspecies

    def make_scan(self,scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen):
        swlog_out = '.true.'
        if(not(swlog)):
            swlog_out = '.false.'

        heating_out = '.true.'
        if(not(heating)):
            heating_out = '.false.'

        eigen_out = '.true.'
        if(not(eigen)):
            eigen_out = '.false.'

        tempscan = {'scan_type':int(scan_type),
                    'scan_style':int(scan_style),
                    'swi':swi,
                    'swf':swf,
                    'swlog':swlog_out,
                    'ns':int(ns),
                    'nres':int(nres),
                    'heating':heating_out,
                    'eigen':eigen_out
        }

        if(len(self.scan_inputs[0]) == 0):
            self.scan_inputs = [tempscan]
        else:
            self.scan_inputs.append(tempscan)

    def make_guess(self,g_om,g_gam):
        tempguess = {'g_om':g_om,
                     'g_gam':g_gam
        }

        if(len(self.guesses[0]) == 0):
            self.guesses = [tempguess]
        else:
            self.guesses.append(tempguess)

    def check_input_for_errors(self):
        #will not halt things,  but will print warnings

        spec1keysthatshouldbe1 = ['tauS','muS','Qs','Ds']
        _tolerance = 0.0001
        for _k in spec1keysthatshouldbe1:
            if(abs(self.species[0][_k]-1.) < _tolerance):
                print("Warning, reference species should have",_k," = 1")
                print("Results are most certianly incorrect!")

    def write_input(self,flnm,outputname,desc='',verbose=False):

        self.check_input_for_errors()

        _replace_input_aux(flnm,verbose=verbose)

        f = open(str(flnm)+'.in', "w")

        f.write('&params\n')
        for key in self.params.keys():
            line = str(key)+'='+str(self.params[key])+'\n'
            f.write(line)
        f.write('dataname=\''+str(self.dataname)+'\'\n')
        f.write('outputname=\''+str(outputname)+'\'\n')
        f.write('/\n\n')

        if(len(self.fpc) != 0):
            f.write('&fpc\n')
            for key in self.fpc.keys():
                line = str(key)+'='+str(self.fpc[key])+'\n'
                f.write(line)
            f.write('/\n\n')

        specidx = 1
        for spec in self.species:
            f.write('&species_'+str(specidx)+'\n')
            for key in spec.keys():
                line = str(key)+'='+str(spec[key])+'\n'
                f.write(line)
            f.write('/\n\n')
            specidx += 1

        if(len(self.maps) != 0):
            f.write('&maps\n')
            for key in self.maps.keys():
                line = str(key)+'='+str(self.maps[key])+'\n'
                f.write(line)
            f.write('/\n\n')

        if(len(self.scan_inputs[0]) != 0):
            scanidx = 1
            for scan in self.scan_inputs:
                f.write('&scan_input_'+str(scanidx)+'\n')
                for key in scan.keys():
                    line = str(key)+'='+str(scan[key])+'\n'
                    f.write(line)
                f.write('/\n\n')
                scanidx += 1

        if(len(self.guesses[0]) != 0):
            gidx = 1
            for guess in self.guesses:
                f.write('&guess_'+str(gidx)+'\n')
                for key in guess.keys():
                    line = str(key)+'='+str(guess[key])+'\n'
                    f.write(line)
                f.write('/\n\n')
                gidx += 1

        f.close()

    def load_from_file(self,flnm):
        paramkeys = ['betap','kperp','kpar','vtp','nspec','nscan','option','nroot_max','use_map','writeOut','dataName','outputName']
        fpckeys = ['vperpmin','vperpmax','vparmin','vparmax','delv','elecdircontribution']
        specieskeys = ['tauS','muS','alphS','Qs','Ds','vvS']
        mapskeys = ['loggridw','omi','omf','gami','gamf','positive_roots']
        scaninputkeys = ['scan_type','scan_style','swi','swf','swlog','ns','nres','heating','eigen']
        guesskeys = ['g_om','g_gam']

        tempfile = open(flnm, "r")
        _tempfile = []
        for line in tempfile:
            _tempfile.append(line)
        tempfile.close()
        tempfile = _tempfile

        _i = 0
        while _i < len(tempfile):
            parse = tempfile[_i].split()

            if(len(parse) == 0):
                parse = ['']

            if(parse[0] == '&params'):
                tempparamdict = {}
                while _i < len(tempfile):
                    line = tempfile[_i].replace(' ', '')
                    parse = line.split('=')
                    if(parse[0]=='/\n'):
                        break
                    if(parse[0] in paramkeys):
                        if(parse[0] in ['nspec','nscan','option','nroot_max']):
                            tempparamdict[parse[0]] = int(parse[1].split('!')[0].replace('\n',''))
                        else:
                            try:
                                tempparamdict[parse[0]] = float(parse[1].split('!')[0].replace('\n',''))
                            except:
                                tempparamdict[parse[0]] = parse[1].split('!')[0].replace('\n','')
                    _i += 1
                self.params = tempparamdict

            if(parse[0] == '&fpc'):
                tempfpcdict = {}
                while _i < len(tempfile):
                    line = tempfile[_i].replace(' ', '')
                    parse = line.split('=')
                    if(parse[0]=='/\n'):
                        break
                    if(parse[0] in fpckeys):
                        try:
                            tempfpcdict[parse[0]] = float(parse[1].split('!')[0].replace('\n',''))
                        except:
                            tempfpcdict[parse[0]] = parse[1].split('!')[0].replace('\n','')
                    _i += 1
                self.fpc = tempfpcdict

            if(parse[0] == '&maps'):
                tempmapdict = {}
                while _i < len(tempfile):
                    line = tempfile[_i].replace(' ', '')
                    parse = line.split('=')
                    if(parse[0]=='/\n'):
                        break
                    if(parse[0] in mapskeys):
                        try:
                            tempmapdict[parse[0]] = float(parse[1].split('!')[0].replace('\n',''))
                        except:
                            tempmapdict[parse[0]] = parse[1].split('!')[0].replace('\n','')
                    _i += 1
                self.maps = tempmapdict

            parse_2 = parse[0].split('_')
            if(parse_2[0] == '&species'):
                tempspecdict = {}
                while _i < len(tempfile):
                    line = tempfile[_i].replace(' ', '')
                    parse = line.split('=')
                    if(parse[0]=='/\n'):
                        break
                    if(parse[0] in specieskeys):
                        try:
                            tempspecdict[parse[0]] = float(parse[1].split('!')[0].replace('\n',''))
                        except:
                            tempspecdict[parse[0]] = parse[1].split('!')[0].replace('\n','')
                    _i += 1
                if(len(self.species[0]) == 0):
                    self.species = [tempspecdict]
                else:
                    self.species.append(tempspecdict)

            if(parse_2[0] == '&scan'):
                tempscandict = {}
                while _i < len(tempfile):
                    line = tempfile[_i].replace(' ', '')
                    parse = line.split('=')
                    if(parse[0]=='/\n'):
                        break
                    if(parse[0] in scaninputkeys):
                        if(parse[0] in ['scan_type','scan_style','ns','nres']):
                            tempscandict[parse[0]] = int(parse[1].split('!')[0].replace('\n',''))
                        else:
                            try:
                                tempscandict[parse[0]] = float(parse[1].split('!')[0].replace('\n',''))
                            except:
                                tempscandict[parse[0]] = parse[1].split('!')[0].replace('\n','')
                    _i += 1
                if(len(self.scan_inputs[0]) == 0):
                    self.scan_inputs = [tempscandict]
                else:
                    self.scan_inputs.append(tempscandict)

            if(parse_2[0] == '&guess'):
                tempguessdict = {}
                while _i < len(tempfile):
                    line = tempfile[_i].replace(' ', '')
                    parse = line.split('=')
                    if(parse[0]=='/\n'):
                        break
                    if(parse[0] in guesskeys):
                        try:
                            tempguessdict[parse[0]] = float(parse[1].split('!')[0].replace('\n',''))
                        except:
                            tempguessdict[parse[0]] = parse[1].split('!')[0].replace('\n','')
                    _i += 1
                if(len(self.guesses[0]) == 0):
                    self.guesses = [tempguessdict]
                else:
                    self.guesses.append(tempguessdict)


            _i += 1


    #main dict
    namelists = {}
    params = {}
    fpc = {}
    species = [{}]
    maps = {}
    scan_inputs = [{}]
    guesses = [{}]

    dataname = 'default'

def _replace_input_aux(inputflnm,verbose=False):
    from os.path import exists

    file_exists = exists(inputflnm+'.in')
    
    if(file_exists):
        if(verbose):
            print("Rewriting file...")
        cmd = 'rm '+str(inputflnm)+'.in'
        os.system(cmd)
    
    else:
        print("File, "+inputflnm+",does not already exist...")
        print("Carrying on...")

def compute_roots(plume_input,inputflnm,outputname,outlog='outlog',verbose=False):
    """
    Computes roots for given input. Uses newton secant method starting with points generated 
    on grid defined by map in plume input

    Paramters
    ---------
    plume_input : class
        class containing input parameters
    inputflnm : string
        name to call input file name
    outputname : string
        name to call output data name
    outlog : string
        filename to write plume output to (print statements- not data)
    verbose : bool
        if true, write print statements

    Returns
    -------
    roots : 1d array (complex)
        found roots
    """

    print("OVERWRITING OPTION; TODO CHECK THAT plume_input IS CORRECT INSTEAD...")
    plume_input.params['option'] = 0
    plume_input.params['use_map'] = '.true.'

    plume_input.write_input(inputflnm,outputname)

    cmd = 'mkdir data/' + plume_input.dataname
    print(cmd)
    os.system(cmd)

    cmd = './plume.e ' + inputflnm + '.in'
    cmd += ' >> ' + outlog
    print(cmd)
    os.system(cmd)

    #kperp,kpar,betap,vtp,wroots(1:2,j),params(1:6,1:nspec)
    rootflnm = 'data/'+str(plume_input.dataname)+'/dispersion_'+outputname+'.roots'
    if(verbose): print("Reading roots from ",rootflnm)

    roots = []
    tempfile = open(rootflnm, "r")
    for line in tempfile:
        parse = line.split()
        temp_om = float(parse[4])
        temp_gam = float(parse[5])
        roots.append(temp_om+temp_gam*1j)
    tempfile.close()

    return np.asarray(roots)
    
def refine_root_and_calc_eigen_guess(plume_input,root,inputflnm,outputname,outlog='outlog',verbose=False):
    """
    Takes specified solution, refines it using newton secant method, computes eigenvalues, and returns data as
    dictionary

    Paramters
    ---------
    plume_input : class
        class containing input parameters
    root : complex
        root to refine
    inputflnm : string
        name to call input file name
    outputname : string
        name to call output data name
    outlog : string
        filename to write plume output to (print statements- not data)
    verbose : bool
        if true, write print statements

    Returns
    -------
    onepointsweep : dict
        dict of found root and associated eigenfunctions
    """
    outputnametemp1 = outputname+'eigenatpoint'
    inputflnmtemp1 = inputflnm+'eigenatpoint'
    
    if(verbose):
        print("OVERWRITING OPTION...")
    plume_input.params['option'] = 1
    plume_input.params['use_map'] = '.false.'
    
    sweepvar = '0'
    midsweepval = plume_input.params['kperp']
    
    scan_type=sweepvar

    scan_style=0
    swi=plume_input.params['kperp']
    swf=plume_input.params['kperp']
    swlog=False
    ns=1
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)
    
    plume_input.guesses = [{}]
    g_om = root.real
    g_gam = root.imag
    plume_input.make_guess(g_om,g_gam)
    
    plume_input.write_input(inputflnmtemp1,outputnametemp1,verbose=verbose)
    cmd = './plume.e ' + inputflnmtemp1 + '.in'
    cmd += ' >> ' + outlog
    if(verbose):
        print(cmd)
    os.system(cmd)
    
    flnmsweep = 'data/'+plume_input.dataname+'/'+outputnametemp1+'_kperp_'+str(int(swi*1000))+'_'+str(int(swf*1000))+'.mode1'
    if(verbose):
        print("Loading ",flnmsweep,"...")
    onepointsweep = load_plume_sweep(flnmsweep,verbose=verbose)
    
    #remove redundant point info
    for _key in onepointsweep:
        onepointsweep[_key] = np.asarray(onepointsweep[_key][0])
    
    return onepointsweep

def compute_fpc_from_root(plume_input,root,inputflnm,outputname,outlog='outlog',cart=False,verbose=False):
    """
    Takes specified solution, refines it using newton secant method, and computes fpc

    Paramters
    ---------
    plume_input : class
        class containing input parameters
    root : complex
        root to refine
    inputflnm : string
        name to call input file name
    outputname : string
        name to call output data name
    outlog : string
        filename to write plume output to (print statements- not data)
    verbose : bool
        if true, write print statements

    Returns
    -------
    cdatafilenames : 1darray (strings)
        names of all fpc/fs1 data
    """

    #use guess to compute fpc
    if(verbose): print("OVERWRITING OPTION AND NUM GUESS AND USE_MAP...")
    if(not(cart)):
        plume_input.params['option'] = 6
    else:
        plume_input.params['option'] = 7
    plume_input.params['use_map'] = '.false.'
    plume_input.params['nroot_max'] = 1

    if(verbose): print("Forcing input file to use root as initial guess")
    plume_input.guesses = [{}]
    g_om = root.real
    g_gam = root.imag
    plume_input.make_guess(g_om,g_gam)

    plume_input.write_input(inputflnm,outputname)

    cmd = 'mkdir data/' + plume_input.dataname
    print(cmd)
    os.system(cmd)

    cmd = './plume.e ' + inputflnm + '.in'
    cmd += ' >> ' + outlog
    print(cmd)
    os.system(cmd)

    cdatafilenames = []
    for _nspec in range(1,len(plume_input.species)+1):
        cdatafilename = 'data/'+str(plume_input.dataname)+'/'+outputname+'.cpar.specie0'+str(_nspec)+'.mode01' #warning: assumes num of species < 10
        cdatafilenames.append(cdatafilename)
        cdatafilename = 'data/'+str(plume_input.dataname)+'/'+outputname+'.cperp.specie0'+str(_nspec)+'.mode01' #warning: assumes num of species < 10
        cdatafilenames.append(cdatafilename)
        cdatafilenames.append('data/'+str(plume_input.dataname)+'/'+outputname+'.df1gyro.real.specie0'+str(_nspec)+'.mode01') #warning: assumes num of species < 10
        cdatafilenames.append('data/'+str(plume_input.dataname)+'/'+outputname+'.df1gyro.imag.specie0'+str(_nspec)+'.mode01') #warning: assumes num of species < 10

    return cdatafilenames

def make_sweeps_that_branch_from_params(plume_input,stylenum,sweepvarkey,sweepmin,sweepmax,root,inputflnm,outputname,outlog='outlog',nsamps=200,verbose=False,use_ps_split_new=True,swlog=True):
    """
    Makes two sweeps that both start at point specified by params namelist and root and branches out in opposite directions.

    Please make sure that your point falls within the sweepmin and sweepmax range!

    Paramters
    ---------
    plume_input : class
        class containing input parameters
    sytlenum : int
        integer of sweep style (see main PLUME readme)
    sweepvarkey : string
        key of variable to be swept over
    sweepmin/sweepmax : float
        min and max value of variable to be swept over
    root : complex
        root to refine
    inputflnm : string
        name to call input file name
    outputname : string
        name to call output data name
    outlog : string
        filename to write plume output to (print statements- not data)
    verbose : bool
        if true, write print statements
    nsamps : int
        number of samples per sweep
    use_ps_split_new : bool
        use new or old power split (warning- must match .new_low_n. in vars.f90 (recompile if changed))
    swlog : bool
        use log spacing between sweep points

    Returns
    -------
    sweep : dict
        dict containing parallel arrays of sweep values
    """
    if(verbose):
        print("OVERWRITING OPTION AND NUM GUESS AND USE_MAP...")
    plume_input.params['option'] = 1
    plume_input.params['use_map'] = '.false.'
    plume_input.params['nroot_max'] = 1

    sweepvar = -1
    midsweepval = 0.

    if(stylenum == 0):
        if(sweepvarkey == 'kperp'):
                sweepvar = '0'
                midsweepval = plume_input.params['kperp']
        elif(sweepvarkey == 'kpar'):
                sweepvar = '1'
                midsweepval = plume_input.params['kpar']
        elif(sweepvarkey == 'betap'):
                sweepvar = '2'
                midsweepval = plume_input.params['betap']
        elif(sweepvarkey == 'vtp'):
                sweepvar = '3'
                midsweepval = plume_input.params['vtp']
        else:
                print("********************************************************************************")
                print("Please input a valid sweepvarkey for this stylenum!*****************************")
                print("********************************************************************************")

    elif(stylenum >= 1):
        if(sweepvarkey == 'tauS'):
                sweepvar = '0'
                midsweepval = plume_input.species[stylenum-1]['tauS']
        elif(sweepvarkey == 'muS'):
                sweepvar = '1'
                midsweepval = plume_input.species[stylenum-1]['muS']
        elif(sweepvarkey == 'alphS'):
                sweepvar = '2'
                midsweepval = plume_input.species[stylenum-1]['alphS']
        elif(sweepvarkey == 'Qs'):
                sweepvar = '3'
                midsweepval = plume_input.species[stylenum-1]['Qs']
        elif(sweepvarkey == 'Ds'):
                sweepvar = '4'
                midsweepval = plume_input.species[stylenum-1]['Ds']
        elif(sweepvarkey == 'vvS'):
                sweepvar = '5'
                midsweepval = plume_input.species[stylenum-1]['vvS']
        else:
                print("********************************************************************************")
                print("Please input a valid sweepvarkey for this style!********************************")
                print("********************************************************************************")

    else:
        print("Please enter a valid stylenum...")


    plume_input.guesses = [{}]
    g_om = root.real
    g_gam = root.imag
    plume_input.make_guess(g_om,g_gam)

    #do first sweep
    outputnametemp1 = outputname+'sweep1'
    inputflnmtemp1 = inputflnm+'sweep1'
    scan_type=sweepvar
    scan_style=stylenum
    swi=midsweepval
    swf=sweepmax
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)
    plume_input.write_input(inputflnmtemp1,outputnametemp1,verbose=verbose)
    cmd = './plume.e ' + inputflnmtemp1 + '.in'
    cmd += ' >> ' + outlog
    if(verbose):
        print(cmd)
    os.system(cmd)

    #do second sweep
    outputnametemp2 = outputname+'sweep2'
    inputflnmtemp2 = inputflnm+'sweep2'
    scan_type=sweepvar
    scan_style=stylenum
    swi=midsweepval
    swf=sweepmin
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)
    plume_input.write_input(inputflnmtemp2,outputnametemp2,verbose=verbose)
    cmd = './plume.e ' + inputflnmtemp2 + '.in'
    cmd += ' >> ' + outlog
    if(verbose):
        print(cmd)
    os.system(cmd)

    #try:
    #load sweeps
    if(stylenum>=1):
        flnmsweep1 = 'data/'+plume_input.dataname+'/'+outputnametemp1+'_'+sweepvarkey+'_'+'s'+str(stylenum)+'_'+str(int(midsweepval*1000))+'_'+str(int(sweepmax*1000))+'.mode1'
    else:
        flnmsweep1 = 'data/'+plume_input.dataname+'/'+outputnametemp1+'_'+sweepvarkey+'_'+str(int(midsweepval*1000))+'_'+str(int(sweepmax*1000))+'.mode1'
    if(verbose):
        print("Loading ",flnmsweep1,"...")
    sweephigh = load_plume_sweep(flnmsweep1,verbose=verbose,use_ps_split_new=use_ps_split_new)

    if(stylenum>=1):
        flnmsweep2 = 'data/'+plume_input.dataname+'/'+outputnametemp2+'_'+sweepvarkey+'_'+'s'+str(stylenum)+'_'+str(int(midsweepval*1000))+'_'+str(int(sweepmin*1000))+'.mode1'
    else:
        flnmsweep2 = 'data/'+plume_input.dataname+'/'+outputnametemp2+'_'+sweepvarkey+'_'+str(int(midsweepval*1000))+'_'+str(int(sweepmin*1000))+'.mode1'
    if(verbose):
        print("Loading ",flnmsweep2,"...")
    sweeplow = load_plume_sweep(flnmsweep2,verbose=verbose,use_ps_split_new=use_ps_split_new)

    if(verbose):
        print("Combining data and returning as 1 sweep...")
    for _key in sweeplow.keys():
        sweeplow[_key] = np.flip(sweeplow[_key])

    sweep = {}
    for _key in sweeplow.keys():
        sweep[_key] = np.concatenate((sweeplow[_key],sweephigh[_key]))
    return sweep

    # except:
    #     print("Could not automatically load and combine sweeps.")
    #     print("Probable cause of error is described below:")
    #     print("If sweepmax and sweepmin and `middle of sweep' val (variable `midsweepval' in function) are not multiples of .001, PLUME will round these numbers")
    #     print("and thus the output file name will be rounded as well. ")
    #     print("It is recommended that one rounds sweepmax and sweepmin to multiples of .001 manually.")

def loadlinfpccepar(filename,verbose=False):
    """
    Loads fpc data for Epar in gyrotropic coordinates

    Paramters
    ---------
    filename : string
        filename to be loaded
    verbose : bool
        if true, write print statements

    Returns
    -------
    linfpcdata : dict
        dict containing fpc data
    """

    if(verbose):print("Opening " + filename)
    try:
        f = open(filename)
    except:
        print("Couldnt open: " + filename)
        return
    line = f.readline()
    line = f.readline()
    line = line.split()
    tau   = float(line[0])
    bi    = float(line[1])
    kpar  = float(line[2]) #kpar rho
    kperp = float(line[3]) #kperp rho
    vts   = float(line[4])
    mu    = float(line[5])
    omega = complex(float(line[6]),float(line[7]))

    line = f.readline()
    line = f.readline()
    line = line.split()
    vperpmin = float(line[0])
    vperpmax = float(line[1])
    vparmin  = float(line[2])
    vparmax  = float(line[3])
    delv     = float(line[4])
    resonant_int = omega/math.sqrt(bi) #calc resonant interval
    species = ''
    if(filename.find('specie02')>=0):
        resonant_int = resonant_int*tau**(.5)*mu**(-.5)
        if(verbose):print("Calculated resonant interval (elec): " + str(resonant_int))
        species = 'elec'
    else:
        if(verbose):print("Calculated resonant interval (ion): " + str(resonant_int))
        species = 'ion'

    resonant_int = resonant_int.real

    vpar = []
    vperp = []
    vparindex = vparmin
    vperpindex = vperpmin
    while(vparindex <= vparmax):
        vpar.append(float(vparindex))
        vparindex += delv
    while(vperpindex <= vperpmax):
        vperp.append(float(vperpindex))
        vperpindex += delv

    line = f.readline()
    C = []
    line = f.readline()
    while (line != ''):
        line = line.split()
        row = []
        for k in line:
            if(math.isnan(float(k))):
                row.append(0.)
            else:
                row.append(float(k))
        C.append(row)
        line = f.readline()

    C = np.asarray(C)

    #index rounding bug fix
    if(len(vpar) < len(C[0])):
        vpar.append(float(vparindex))
    elif(len(vpar) > len(C[0])):
        vpar.pop()
    if(len(vperp) < len(C)):
        vperp.append(float(vperpindex))
    elif(len(vperp) > len(C)):
        vperp.pop()

    linfpcdata = {'CEpar':C,'vpar':vpar,'vperp':vperp,'resonant_int':resonant_int,'vperpmin':vperpmin,'vperpmax':vperpmax,'vparmin':vparmin,'vparmax':vparmax,'delv':delv,'species':species,'omega':omega}

    return linfpcdata

def loadlinfpcceperp(filename,verbose=False):
    """
    Loads fpc data for Eperp in gyrotropic coordinates

    Paramters
    ---------
    filename : string
        filename to be loaded
    verbose : bool
        if true, write print statements

    Returns
    -------
    linfpcdata : dict
        dict containing fpc data
    """

    import copy
    linfpcdata = loadlinfpccepar(filename,verbose=verbose)

    linfpcdata['CEperp'] = copy.deepcopy(linfpcdata['CEpar'])
    del linfpcdata['CEpar']

    return linfpcdata

def loadlinfpcgyro_dist(filenamereal,filenameimag):
    """
    Loads fpc data for fs1 in gyrotropic coordinates

    Paramters
    ---------
    filename : string
        filename to be loaded
    verbose : bool
        if true, write print statements

    Returns
    -------
    outdict : dict
        dict containing fs1 data
    """

    realpartdict = loadlinfpccepar(filenamereal)
    imagpartdict = loadlinfpccepar(filenameimag)

    outdict = {}

    for key in imagpartdict.keys():
        outdict[key] = imagpartdict[key]
    outdict['im_f'] = outdict['CEpar'] 

    for key in realpartdict.keys():
        outdict[key] = realpartdict[key]
    outdict['re_f'] = outdict['CEpar'] 
    del outdict['CEpar']

    return outdict

def loadlinfpccart_dist(filenamereal,filenameimag,idxoffset=0):
    """
    Loads fpc data for fs1 in cartesian coordintates

    Paramters
    ---------
    filename : string
        filename to be loaded
    verbose : bool
        if true, write print statements
    idxoffset : int
        integer used to attempt to fix off by one errors when generating grid. Not to be used by user.

    Returns
    -------
    outdict : dict
        dict containing fpc data
    """

    realpartdict = loadlinfpccart(filenamereal)
    imagpartdict = loadlinfpccart(filenameimag)

    outdict = {}

    for key in imagpartdict.keys():
        outdict[key] = imagpartdict[key]

    outdict['im_fvxvy'] = outdict['CEparvxvy'] 
    outdict['im_fvxvz'] = outdict['CEparvxvz'] 
    outdict['im_fvyvz'] = outdict['CEparvyvz']

    for key in realpartdict.keys():
        outdict[key] = realpartdict[key]
    outdict['re_fvxvy'] = outdict['CEparvxvy'] 
    outdict['re_fvxvz'] = outdict['CEparvxvz'] 
    outdict['re_fvyvz'] = outdict['CEparvyvz']

    del outdict['CEparvxvy']
    del outdict['CEparvxvz']
    del outdict['CEparvyvz']

    return outdict

def loadlinfpccart(filename,idxoffset=0,verbose=False):
    """
    Loads fpc data for fs1 in cartesian coordintates
    assumes equal bounds in vx vy and vz direction

    Paramters
    ---------
    filename : string
        filename to be loaded
    verbose : bool
        if true, write print statements
    idxoffset : int
        integer used to attempt to fix off by one errors when generating grid. Not to be used by user.

    Returns
    -------
    linfpcdata : dict
        dict containing fpc data
    """

    if(abs(idxoffset) > 2):
        print("Error! Couldn't load  ",filename)
        return
    try:
        if(verbose):print("Opening " + filename)
        try:
            f = open(filename)
        except:
            print("Couldnt open: " + filename)
            return
        line = f.readline()
        line = f.readline()
        line = line.split()
        tau   = float(line[0])
        bi    = float(line[1])
        kpar  = float(line[2]) #kpar rho
        kperp = float(line[3]) #kperp rho
        vts   = float(line[4])
        mu    = float(line[5])
        omega = complex(float(line[6]),float(line[7]))

        line = f.readline()
        line = f.readline()
        line = line.split()
        vxmin = float(line[0])
        vxmax = float(line[1])
        vymin  = float(line[2])
        vymax  = float(line[3])
        vzmin = float(line[4])
        vzmax = float(line[5])
        delv     = float(line[6])
        resonant_int = omega/math.sqrt(bi) #calc resonant interval
        species = ''
        if(filename.find('specie02')>=0):
            resonant_int = resonant_int*tau**(.5)*mu**(-.5)
            if(verbose): print("Calculated resonant interval (elec): " + str(resonant_int))
            species = 'elec'
        else:
            if(verbose): print("Calculated resonant interval (ion): " + str(resonant_int))
            species = 'ion'

        resonant_int = resonant_int.real

        ivxmin=int(round(vxmin/delv))
        ivxmax=int(round(vxmax/delv))
        ivymin=int(round(vymin/delv))
        ivymax=int(round(vymax/delv))
        ivzmin=int(round(vzmin/delv))
        ivzmax=int(round(vzmax/delv))

        numstepvx = ivxmax-ivxmin+1-idxoffset
        numstepvy = ivymax-ivymin+1-idxoffset
        numstepvz = ivzmax-ivzmin+1-idxoffset

        vx_xy = np.zeros((numstepvx,numstepvy))
        vy_xy = np.zeros((numstepvx,numstepvy))

        vx_xz = np.zeros((numstepvx,numstepvz))
        vz_xz = np.zeros((numstepvx,numstepvz))

        vy_yz = np.zeros((numstepvy,numstepvz))
        vz_yz = np.zeros((numstepvy,numstepvz))

        vx = []
        vy = []
        vz = []
        vxindex = vxmin
        vyindex = vymin
        vzindex = vzmin

        _i = 0
        vx.append(float(vxindex))
        while(_i < numstepvx):
            vxindex += delv
            vx.append(float(vxindex))
            _i += 1

        _i = 0
        vy.append(float(vyindex))
        while(_i < numstepvy):
            vyindex += delv
            vy.append(float(vyindex))
            _i += 1

        vz.append(float(vzindex))
        _i = 0
        while(_i < numstepvz):
            vzindex += delv
            vz.append(float(vzindex))
            _i += 1
            
        line = f.readline()
        Cvxvy = []
        line = f.readline()
        linecounter = 1
        while (linecounter < len(vx)):
            line = line.split()
            row = []
            rowcounter = 0
            for k in line:
                if(math.isnan(float(k))):
                    row.append(0.)
                else:
                    row.append(float(k))
                vx_xy[linecounter-1,rowcounter] = vx[linecounter-1]
                vy_xy[linecounter-1,rowcounter] = vy[rowcounter]
                rowcounter += 1
            Cvxvy.append(row)
            line = f.readline()
            linecounter += 1
        Cvxvy = np.asarray(Cvxvy)

        Cvxvz = []
        linecounter = 1
        line = f.readline()
        while (linecounter < len(vx)):
            line = line.split()
            row = []
            rowcounter = 0
            for k in line:
                if(math.isnan(float(k))):
                    row.append(0.)
                else:
                    row.append(float(k))
                vx_xz[linecounter-1,rowcounter] = vx[linecounter-1]
                vz_xz[linecounter-1,rowcounter] = vz[rowcounter]
                rowcounter += 1
            Cvxvz.append(row)
            line = f.readline()
            linecounter += 1
        Cvxvz = np.asarray(Cvxvz)

        Cvyvz = []
        linecounter = 1
        line = f.readline()
        while (linecounter < len(vy)):
            line = line.split()
            row = []
            rowcounter = 0
            for k in line:
                if(math.isnan(float(k))):
                    row.append(0.)
                else:
                    row.append(float(k))
                vy_yz[linecounter-1,rowcounter] = vy[linecounter-1]
                vz_yz[linecounter-1,rowcounter] = vz[rowcounter]
                rowcounter += 1
            Cvyvz.append(row)
            line = f.readline()
            linecounter += 1
        Cvyvz = np.asarray(Cvyvz)
        vy_yz = vy_yz.T
        vz_yz = vz_yz.T

        linfpcckeyname = 'CEpar'
        if('perp1' in filename):
            linfpcckeyname = 'CEperp1'
        elif('perp2' in filename):
            linfpcckeyname = 'CEperp2'

        linfpcdata = {linfpcckeyname+'vxvy':Cvxvy,linfpcckeyname+'vxvz':Cvxvz,linfpcckeyname+'vyvz':Cvyvz,
                      'vx':vx,'vy':vy,'vz':vz,'resonant_int':resonant_int,'vxmin':vxmin,'vxmax':vxmax,'vymin':vymin,'vymax':vymax,'vzmin':vzmin,'vzmax':vzmax,
                      'delv':delv,'species':species,'omega':omega,
                      'vx_xy':vx_xy,'vy_xy':vy_xy,'vx_xz':vx_xz,'vz_xz':vz_xz,'vy_yz':vy_yz,'vz_yz':vz_yz}

        return linfpcdata
    except:
        return loadlinfpccart(filename,idxoffset=idxoffset+1)


def load_plume_sweep(flnm, nspec = 0, heating = False, eigen = False, verbose = False, use_ps_split_new=True):

    #Attempts to figure out what is output based on the number of values in a line-> Technically not possible
    # for all values of nspec given there exists a least common multiple where nspec*(num_additional_elements_if_eigen_is_true_and_new_low_n_is_false) = nspec^prime(num_additional_elements_if_both_are_true)
    # where nspec != nspec^prime

    #verbose is currently not used....

    if(use_ps_split_new == False):
        print("Warning- the most current version of PLUME forces use_ps_split_new=True now!")
        return

    nspec = int(round(nspec)) #note, python truncates when casting to int, which could be an issue if due to floating division if we get something like *.9999999999999, so we round before casting as int

    if(nspec == 0):
        f = open(flnm)
        line = f.readline()
        nline = len(line.split())

        if(nline > 12):
            eigen = True

        skipforedgecase = False
        if(nline == 60): #Given that having nspec = 2 heating=True and Eigen = True is so common, I manually check \mu_r (mass_ref/mass_ref) in the conflicting case of Eigen = true, heating = False, nspec = 3 (in nspec = 2 case, element 44 (at index 43) is the power absoprtion by p2ld_zz (which is super unlikely to be exactly +1, and in the nspec=3 case, element 44 (at index 43) is mu_r=1 (unless the user creates a bad PLUME input....))
                         #in summary, unless the user has a bad input, it will either correctly load, or ask the user for more input, which is fine.
                         #      if the user has bad input (which we check for!)
            #Add elements where relevant power output should be if nspec=2 heating=True and Eigen=True (very very improbable this equation returns true if this setup is not used and nline=60) <- (fun note: computers assume statistically impossible but technically possible things all the time- for example checksums cryptographic hash collisiones, pseudo random number generators, hardware failure, quantum tunneling....)
            _tol = .000000001
            if(abs(float(line.split()[43])-1) > _tol): #_tol is for small potential floating point rounding error
                nspec = 2
                eigen = True
                heating = True
                skipforedgecase = True #If the above is true, we know what it is and should keep going
                noutperspec = 15

        if(not(skipforedgecase)):
            #first check if two of these are true, then we can't determine what the output is given just the number of elements per line (except by checking specific values like if the sum of species power subelements equals total species power (approximately since it's technically missing terms) TODO: see above example and write for more? <- lot of work to not be used often, just the one common edge case is fine for now.)
            if(int((nline - 18) % (8+6) == 0)+int((nline - 6) % (7+6) == 0)+int((nline - 18) % (15+6) == 0)+int((nline - 6) % (6)==0) > 1):
                print("Could not automatically determine nspec, heating=T/F, and eigen=T/F (as number line elements, which we use to determine these values, could be produced by more than one combination of nspec heating and eigen)")
                print("Please call this function again and pass as optional paramters nspec=*val*, heating=True/False, eigen=True/False")
                print('Number of elements in first line: ',nline)
                print("Possible if Eigen is true and heat is false (bool->1 or 0) | Possible if Eigen is false and heat is true (bool->1 or 0) | Possible if Eigen is true and heat is true (bool->1 or 0)")
                print("int((nline - 18) % (8+6) == 0):",int((nline - 18) % (8+6) == 0), "    |int((nline - 6) % (7+6) == 0)",int((nline - 6) % (7+6) == 0),"   |int((nline - 6) % (6)==0)",int((nline - 6) % (6)==0))
                f.close()

                return

            if((nline - 18) % (8+6) == 0): #some outputs at front, 6nspec at end, and the middle is determined by nspec and heating/eigen (if eigen = True then  +8, if heating = True then +7)
                eigen = True
                heating = False
                nspec = int(round((nline - 18)/(8+6))) #note, python truncates when casting to int, which could be an issue if due to floating division if we get something like *.9999999999999, so we round before casting as int
                
                noutperspec = 8 #note this is the *additional* number of output per spec accounting for the 6 that is always there

            elif((nline - 6) % (7+6) == 0): #some outputs at front, 6nspec at end, and the middle is determined by nspec and heating/eigen (if eigen = True then  +8, if heating = True then +7)
                eigen = False
                heating = True
                nspec = int(round(nline - 6)/(7+6)) #note, python truncates when casting to int, which could be an issue if due to floating division if we get something like *.9999999999999, so we round before casting as int
                noutperspec = 7

            elif((nline - 18) % (15+6) == 0): #some outputs at front, 6nspec at end, and the middle is determined by nspec and heating/eigen (if eigen = True then  +8, if heating = True then +7)
                eigen = True
                heating = True
                nspec = int(round(nline - 18)/(15+6)) #note, python truncates when casting to int, which could be an issue if due to floating division if we get something like *.9999999999999, so we round before casting as int
                noutperspec = 15

            elif((nline - 6) % (6) == 0):
                eigen = False
                heating = False
                nspec = int(round(nline - 6)/(6)) #note, python truncates when casting to int, which could be an issue if due to floating division if we get something like *.9999999999999, so we round before casting as int
                noutperspec = 0

            else:
                print("Error, could not determine input format from output! It seems there are not the correct number of elements in a line for any setup...")

                return

        f.close()

    else:

        f = open(flnm)
        line = f.readline()
        nline = len(line.split())
        f.close()

        if(heating == False and eigen == True):
            noutperspec = 8
        elif(heating == True and eigen == False):
            noutperspec = 7
        elif(heating == True and eigen == True):
            noutperspec = 15
        elif(heating == False and eigen == False):
            noutperspec = 0

    #make python dictionary
    plume_sweep = {
            "kperp": [],
            "kpar": [],
            "betap": [],
            "vtp": [],
            "w": [],
            "g": []
        }
    if(eigen):
        plume_sweep["bxr"] = []
        plume_sweep["bxi"] = []
        plume_sweep["byr"] = []
        plume_sweep["byi"] = []
        plume_sweep["bzr"] = []
        plume_sweep["bzi"] = []
        plume_sweep["exr"] = []
        plume_sweep["exi"] = []
        plume_sweep["eyr"] = []
        plume_sweep["eyi"] = []
        plume_sweep["ezr"] = []
        plume_sweep["ezi"] = []
        for _i in range(0,nspec):
            plume_sweep["ux"+str(_i+1)+"r"] = []
            plume_sweep["ux"+str(_i+1)+"i"] = []
            plume_sweep["uy"+str(_i+1)+"r"] = []
            plume_sweep["uy"+str(_i+1)+"i"] = []
            plume_sweep["uz"+str(_i+1)+"r"] = []
            plume_sweep["uz"+str(_i+1)+"i"] = []

        for _i in range(0,nspec):
            plume_sweep["n"+str(_i+1)+"r"] = []
            plume_sweep["n"+str(_i+1)+"i"] = []

    if(heating):
        for _i in range(0,nspec):
            plume_sweep["p"+str(_i+1)] = []
        for _i in range(0,nspec):
            plume_sweep["ps"+str(_i+1)+"ttd_yy"] = []
        for _i in range(0,nspec):
            plume_sweep["p"+str(_i+1)+"ttd_yz"] = []
        for _i in range(0,nspec):
            plume_sweep["p"+str(_i+1)+"ld_zy"] = []
        for _i in range(0,nspec):
            plume_sweep["p"+str(_i+1)+"ld_zz"] = []
        for _i in range(0,nspec):
            plume_sweep["p"+str(_i+1)+"n_eq_0"] = []
        for _i in range(0,nspec):
            plume_sweep["p"+str(_i+1)+"cd_n_pm"] = []

    for _i in range(0,nspec):
        plume_sweep["p"+str(_i+1)+"tau"] = []
        plume_sweep["p"+str(_i+1)+"mu"] = []
        plume_sweep["p"+str(_i+1)+"alph"] = []
        plume_sweep["p"+str(_i+1)+"q"] = []
        plume_sweep["p"+str(_i+1)+"D"] = []
        plume_sweep["p"+str(_i+1)+"vv"] = []

    f = open(flnm)

    line = f.readline()
    while (line != ''):
        line = line.split()
        if(len(line) > 0):
            #NOTE THIS README INDEXING STARTS AT 1 BUT PYTHON INDEX STARTS AT ZERO so the leftmost number in each equation is one less here than in the PLUME output readme!
            plume_sweep["kperp"].append(float(line[0]))
            plume_sweep["kpar"].append(float(line[1]))
            plume_sweep["betap"].append(float(line[2]))
            plume_sweep["vtp"].append(float(line[3]))
            plume_sweep["w"].append(float(line[4]))
            plume_sweep["g"].append(float(line[5]))

            if(eigen):
                plume_sweep["bxr"].append(float(line[6]))
                plume_sweep["bxi"].append(float(line[7]))
                plume_sweep["byr"].append(float(line[8]))
                plume_sweep["byi"].append(float(line[9]))
                plume_sweep["bzr"].append(float(line[10]))
                plume_sweep["bzi"].append(float(line[11]))
                plume_sweep["exr"].append(float(line[12]))
                plume_sweep["exi"].append(float(line[13]))
                plume_sweep["eyr"].append(float(line[14]))
                plume_sweep["eyi"].append(float(line[15]))
                plume_sweep["ezr"].append(float(line[16]))
                plume_sweep["ezi"].append(float(line[17]))
                for _i in range(1,nspec+1):
                    plume_sweep["ux"+str(_i)+"r"].append(float(line[18+6*(_i-1)]))
                    plume_sweep["ux"+str(_i)+"i"].append(float(line[19+6*(_i-1)]))
                    plume_sweep["uy"+str(_i)+"r"].append(float(line[20+6*(_i-1)]))
                    plume_sweep["uy"+str(_i)+"i"].append(float(line[21+6*(_i-1)]))
                    plume_sweep["uz"+str(_i)+"r"].append(float(line[22+6*(_i-1)]))
                    plume_sweep["uz"+str(_i)+"i"].append(float(line[23+6*(_i-1)]))

                for _i in range(1,nspec+1):
                    plume_sweep["n"+str(_i)+"r"].append(float(line[18+6*(nspec)+2*(_i-1)]))
                    plume_sweep["n"+str(_i)+"i"].append(float(line[19+6*(nspec)+2*(_i-1)]))

            if(heating):
                for _i in range(1,nspec+1):
                    plume_sweep["p"+str(_i)].append(float(line[17+(8-6*int(not(eigen)))*(nspec)+(_i)-12*int(not(eigen))]))
                for _i in range(1,nspec+1):
                    plume_sweep["ps"+str(_i)+"ttd_yy"].append(float(line[17+(9-6*int(not(eigen)))*(nspec)+(_i)-12*int(not(eigen))]))
                for _i in range(1,nspec+1):
                    plume_sweep["p"+str(_i)+"ttd_yz"].append(float(line[17+(10-6*int(not(eigen)))*(nspec)+(_i)-12*int(not(eigen))]))
                for _i in range(1,nspec+1):
                    plume_sweep["p"+str(_i)+"ld_zy"].append(float(line[17+(11-6*int(not(eigen)))*(nspec)+(_i)-12*int(not(eigen))]))
                for _i in range(1,nspec+1):
                    plume_sweep["p"+str(_i)+"ld_zz"].append(float(line[17+(12-6*int(not(eigen)))*(nspec)+(_i)-12*int(not(eigen))]))
                for _i in range(1,nspec+1):
                    plume_sweep["p"+str(_i)+"n_eq_0"].append(float(line[17+(13-6*int(not(eigen)))*(nspec)+(_i)-12*int(not(eigen))]))
                for _i in range(1,nspec+1):
                    plume_sweep["p"+str(_i)+"cd_n_pm"].append(float(line[17+(14-6*int(not(eigen)))*(nspec)+(_i)-12*int(not(eigen))]))

            for _i in range(1,nspec+1):
                plume_sweep["p"+str(_i)+"tau"].append(float(line[18+noutperspec*nspec+6*(_i-1)-12*int(not(eigen))]))
                plume_sweep["p"+str(_i)+"mu"].append(float(line[19+noutperspec*nspec+6*(_i-1)-12*int(not(eigen))]))
                plume_sweep["p"+str(_i)+"alph"].append(float(line[20+noutperspec*nspec+6*(_i-1)-12*int(not(eigen))]))
                plume_sweep["p"+str(_i)+"q"].append(float(line[21+noutperspec*nspec+6*(_i-1)-12*int(not(eigen))]))
                plume_sweep["p"+str(_i)+"D"].append(float(line[22+noutperspec*nspec+6*(_i-1)-12*int(not(eigen))]))
                plume_sweep["p"+str(_i)+"vv"].append(float(line[23+noutperspec*nspec+6*(_i-1)-12*int(not(eigen))]))

        line = f.readline()

    for key in plume_sweep.keys():
            plume_sweep[key] = np.asarray(plume_sweep[key])


    return plume_sweep



#TODO: REMOVE after making sure it is obsolete! 
def load_plume_sweep_orig(flnm,verbose=False,use_ps_split_new=True):
    """
    Load data from plume sweep

    Assumes 2 species (TODO: generalize...)

    Parameters
    ----------
    flnm : str
        path to sweep to be loaded
    verbose : bool
        if true, write print statements
    use_ps_split_new : bool
        use new or old power split (warning- must match .new_low_n. in vars.f90 (recompile if changed))

    Returns
    -------
    plume_sweep : dict
        dictionary of data related to plume
    """
    
    #TODO: check if number of columns to see if user is trying to load more than 2 species or is using the wrong ps split!!!

    if(verbose):
        print("WARNING: assuming 2 species...\n TODO: write load_plume_sweep_nspec...")
        if(not(use_ps_split_new)):
            print("WARNING: assuming low_n is true (rather than new_low_n) (pass use_ps_split_new=True to change this)")
            print("If unsure, please check vars.f90 (change requires recompile)")
            print("If both are true, new_low_n is used")
        else:
            print("WARNING: assuming new_low_n is true (rather than low_n)")
            print("If unsure, please check vars.f90 (change requires recompile by calling 'make clean' then 'make' in main dir')")
            print("If both are true, new_low_n is used")

    f = open(flnm)

    if(not(use_ps_split_new)):
        plume_sweep = {
            "kperp": [],
            "kpar": [],
            "betap": [],
            "vtp": [],
            "w": [],
            "g": [],
            "bxr": [],
            "bxi": [],
            "byr": [],
            "byi": [],
            "bzr": [],
            "bzi": [],
            "exr": [],
            "exi": [],
            "eyr": [],
            "eyi": [],
            "ezr": [],
            "ezi": [],
            "ux1r": [],
            "ux1i": [],
            "uy1r": [],
            "uy1i": [],
            "uz1r": [],
            "uz1i": [],
            "ux2r": [],
            "ux2i": [],
            "uy2r": [],
            "uy2i": [],
            "uz2r": [],
            "uz2i": [],
            
            "n1r": [],
            "n1i": [],
            "n2r": [],
            "n2i": [],
            "ps1": [],
            "ps2": [],
            "p1ld": [],
            "p1ttd": [],
            "p1n0": [],
            "p1cd": [],
            "p2ld": [],
            "p2ttd": [],
            "p2n0": [],
            "p2cd": [],

            "p1tau": [],
            "p1mu": [],
            "p1alph": [],
            "p1q": [],
            "p1D": [],
            "p1vv": [],

            "p2tau": [],
            "p2mu": [],
            "p2alph": [],
            "p2q": [],
            "p2D": [],
            "p2vv": []
        }

        line = f.readline()
        while (line != ''):
            line = line.split()
            if(len(line) > 0):
                plume_sweep['kperp'].append(float(line[0]))
                plume_sweep['kpar'].append(float(line[1]))
                plume_sweep['betap'].append(float(line[2]))
                plume_sweep['vtp'].append(float(line[3]))
                plume_sweep['w'].append(float(line[4]))
                plume_sweep['g'].append(float(line[5]))
                plume_sweep['bxr'].append(float(line[6]))
                plume_sweep['bxi'].append(float(line[7]))
                plume_sweep['byr'].append(float(line[8]))
                plume_sweep['byi'].append(float(line[9]))
                plume_sweep['bzr'].append(float(line[10]))
                plume_sweep['bzi'].append(float(line[11]))
                plume_sweep['exr'].append(float(line[12]))
                plume_sweep['exi'].append(float(line[13]))
                plume_sweep['eyr'].append(float(line[14]))
                plume_sweep['eyi'].append(float(line[15]))
                plume_sweep['ezr'].append(float(line[16]))
                plume_sweep['ezi'].append(float(line[17]))
                plume_sweep['ux1r'].append(float(line[18]))
                plume_sweep['ux1i'].append(float(line[19]))
                plume_sweep['uy1r'].append(float(line[20]))
                plume_sweep['uy1i'].append(float(line[21]))
                plume_sweep['uz1r'].append(float(line[22]))
                plume_sweep['uz1i'].append(float(line[23]))
                plume_sweep['ux2r'].append(float(line[24]))
                plume_sweep['ux2i'].append(float(line[25]))
                plume_sweep['uy2r'].append(float(line[26]))
                plume_sweep['uy2i'].append(float(line[27]))
                plume_sweep['uz2r'].append(float(line[28]))
                plume_sweep['uz2i'].append(float(line[29]))
                
                plume_sweep['n1r'].append(float(line[30]))
                plume_sweep['n1i'].append(float(line[31]))
                plume_sweep['n2r'].append(float(line[32]))
                plume_sweep['n2i'].append(float(line[33]))
                plume_sweep['ps1'].append(float(line[34]))
                plume_sweep['ps2'].append(float(line[35]))
                plume_sweep['p1ld'].append(float(line[36]))
                plume_sweep['p1ttd'].append(float(line[37]))
                plume_sweep['p1n0'].append(float(line[38]))
                plume_sweep['p1cd'].append(float(line[39]))
                plume_sweep['p2ld'].append(float(line[40]))
                plume_sweep['p2ttd'].append(float(line[41]))
                plume_sweep['p2n0'].append(float(line[42]))
                plume_sweep['p2cd'].append(float(line[43]))

                plume_sweep['p1tau'].append(float(line[44]))
                plume_sweep['p1mu'].append(float(line[45]))
                plume_sweep['p1alph'].append(float(line[46]))
                plume_sweep['p1q'].append(float(line[47]))
                plume_sweep['p1D'].append(float(line[48]))
                plume_sweep['p1vv'].append(float(line[49]))

                plume_sweep['p2tau'].append(float(line[50]))
                plume_sweep['p2mu'].append(float(line[51]))
                plume_sweep['p2alph'].append(float(line[52]))
                plume_sweep['p2q'].append(float(line[53]))
                plume_sweep['p2D'].append(float(line[54]))
                plume_sweep['p2vv'].append(float(line[55]))
            line = f.readline()

        for key in plume_sweep.keys():
            plume_sweep[key] = np.asarray(plume_sweep[key])

        #normalize B
        plume_sweep['bxr'] = plume_sweep['bxr']*plume_sweep['vtp']
        plume_sweep['byr'] = plume_sweep['byr']*plume_sweep['vtp']
        plume_sweep['bzr'] = plume_sweep['bzr']*plume_sweep['vtp']
        plume_sweep['bxi'] = plume_sweep['bxi']*plume_sweep['vtp']
        plume_sweep['byi'] = plume_sweep['byi']*plume_sweep['vtp']
        plume_sweep['bzi'] = plume_sweep['bzi']*plume_sweep['vtp']
        
        return plume_sweep
        
    else:
        plume_sweep = {
            "kperp": [],
            "kpar": [],
            "betap": [],
            "vtp": [],
            "w": [],
            "g": [],
            "bxr": [],
            "bxi": [],
            "byr": [],
            "byi": [],
            "bzr": [],
            "bzi": [],
            "exr": [],
            "exi": [],
            "eyr": [],
            "eyi": [],
            "ezr": [],
            "ezi": [],
            "ux1r": [],
            "ux1i": [],
            "uy1r": [],
            "uy1i": [],
            "uz1r": [],
            "uz1i": [],
            "ux2r": [],
            "ux2i": [],
            "uy2r": [],
            "uy2i": [],
            "uz2r": [],
            "uz2i": [],
            
            "n1r": [],
            "n1i": [],
            "n2r": [],
            "n2i": [],
            "ps1": [], #power into species 1
            "ps2": [],
            "p1ttd1": [], #diagnol term in tensor; main transit time damping damping term
            "p1ttd2": [], #cross term in tensor; typically small
            "p1ld1": [], #diagnol term in tensor; main landau damping term
            "p1ld2": [], #cross term in tensor; typically small
            "p1n0": [],
            "p1cd": [],
            "p2ld1": [],
            "p2ld2": [],
            "p2ttd1": [],
            "p2ttd2": [],
            "p2n0": [],
            "p2cd": [],

            "p1tau": [],
            "p1mu": [],
            "p1alph": [],
            "p1q": [],
            "p1D": [],
            "p1vv": [],

            "p2tau": [],
            "p2mu": [],
            "p2alph": [],
            "p2q": [],
            "p2D": [],
            "p2vv": []
        }

        line = f.readline()

        while (line != ''):
            line = line.split()
            if(len(line) > 0):
                plume_sweep['kperp'].append(float(line[0]))
                plume_sweep['kpar'].append(float(line[1]))
                plume_sweep['betap'].append(float(line[2]))
                plume_sweep['vtp'].append(float(line[3]))
                plume_sweep['w'].append(float(line[4]))
                plume_sweep['g'].append(float(line[5]))
                plume_sweep['bxr'].append(float(line[6]))
                plume_sweep['bxi'].append(float(line[7]))
                plume_sweep['byr'].append(float(line[8]))
                plume_sweep['byi'].append(float(line[9]))
                plume_sweep['bzr'].append(float(line[10]))
                plume_sweep['bzi'].append(float(line[11]))
                plume_sweep['exr'].append(float(line[12]))
                plume_sweep['exi'].append(float(line[13]))
                plume_sweep['eyr'].append(float(line[14]))
                plume_sweep['eyi'].append(float(line[15]))
                plume_sweep['ezr'].append(float(line[16]))
                plume_sweep['ezi'].append(float(line[17]))
                plume_sweep['ux1r'].append(float(line[18]))
                plume_sweep['ux1i'].append(float(line[19]))
                plume_sweep['uy1r'].append(float(line[20]))
                plume_sweep['uy1i'].append(float(line[21]))
                plume_sweep['uz1r'].append(float(line[22]))
                plume_sweep['uz1i'].append(float(line[23]))
                plume_sweep['ux2r'].append(float(line[24]))
                plume_sweep['ux2i'].append(float(line[25]))
                plume_sweep['uy2r'].append(float(line[26]))
                plume_sweep['uy2i'].append(float(line[27]))
                plume_sweep['uz2r'].append(float(line[28]))
                plume_sweep['uz2i'].append(float(line[29]))
                
                plume_sweep['n1r'].append(float(line[30]))
                plume_sweep['n1i'].append(float(line[31]))
                plume_sweep['n2r'].append(float(line[32]))
                plume_sweep['n2i'].append(float(line[33]))
                plume_sweep['ps1'].append(float(line[34]))
                plume_sweep['ps2'].append(float(line[35]))
                plume_sweep['p1ttd1'].append(float(line[36]))
                plume_sweep['p1ttd2'].append(float(line[37]))
                plume_sweep['p1ld1'].append(float(line[38]))
                plume_sweep['p1ld2'].append(float(line[39]))
                plume_sweep['p1n0'].append(float(line[40]))
                plume_sweep['p1cd'].append(float(line[41]))
                plume_sweep['p2ttd1'].append(float(line[42]))
                plume_sweep['p2ttd2'].append(float(line[43]))
                plume_sweep['p2ld1'].append(float(line[44]))
                plume_sweep['p2ld2'].append(float(line[45]))
                plume_sweep['p2n0'].append(float(line[46]))
                plume_sweep['p2cd'].append(float(line[47]))

                plume_sweep['p1tau'].append(float(line[48]))
                plume_sweep['p1mu'].append(float(line[49]))
                plume_sweep['p1alph'].append(float(line[50]))
                plume_sweep['p1q'].append(float(line[51]))
                plume_sweep['p1D'].append(float(line[52]))
                plume_sweep['p1vv'].append(float(line[53]))

                plume_sweep['p2tau'].append(float(line[54]))
                plume_sweep['p2mu'].append(float(line[55]))
                plume_sweep['p2alph'].append(float(line[56]))
                plume_sweep['p2q'].append(float(line[57]))
                plume_sweep['p2D'].append(float(line[58]))
                plume_sweep['p2vv'].append(float(line[59]))
            line = f.readline()

        for key in plume_sweep.keys():
            plume_sweep[key] = np.asarray(plume_sweep[key])

        #normalize B
        plume_sweep['bxr'] = plume_sweep['bxr']*plume_sweep['vtp']
        plume_sweep['byr'] = plume_sweep['byr']*plume_sweep['vtp']
        plume_sweep['bzr'] = plume_sweep['bzr']*plume_sweep['vtp']
        plume_sweep['bxi'] = plume_sweep['bxi']*plume_sweep['vtp']
        plume_sweep['byi'] = plume_sweep['byi']*plume_sweep['vtp']
        plume_sweep['bzi'] = plume_sweep['bzi']*plume_sweep['vtp']

        return plume_sweep

#TODO: make this (and other relevant functions) work for N species!!!
def load_plume_sweep_orig(flnm,verbose=False,use_ps_split_new=True):
    """
    Load data from plume sweep

    Assumes 2 species (TODO: generalize...)

    Parameters
    ----------
    flnm : str
        path to sweep to be loaded
    verbose : bool
        if true, write print statements
    use_ps_split_new : bool
        use new or old power split (warning- must match .new_low_n. in vars.f90 (recompile if changed))

    Returns
    -------
    plume_sweep : dict
        dictionary of data related to plume
    """
    
    #TODO: check if number of columns to see if user is trying to load more than 2 species or is using the wrong ps split!!!

    if(verbose):
        print("WARNING: assuming 2 species...\n TODO: write load_plume_sweep_nspec...")
        if(not(use_ps_split_new)):
            print("WARNING: assuming low_n is true (rather than new_low_n) (pass use_ps_split_new=True to change this)")
            print("If unsure, please check vars.f90 (change requires recompile)")
            print("If both are true, new_low_n is used")
        else:
            print("WARNING: assuming new_low_n is true (rather than low_n)")
            print("If unsure, please check vars.f90 (change requires recompile by calling 'make clean' then 'make' in main dir')")
            print("If both are true, new_low_n is used")

    f = open(flnm)

    if(not(use_ps_split_new)):
        plume_sweep = {
            "kperp": [],
            "kpar": [],
            "betap": [],
            "vtp": [],
            "w": [],
            "g": [],
            "bxr": [],
            "bxi": [],
            "byr": [],
            "byi": [],
            "bzr": [],
            "bzi": [],
            "exr": [],
            "exi": [],
            "eyr": [],
            "eyi": [],
            "ezr": [],
            "ezi": [],
            "ux1r": [],
            "ux1i": [],
            "uy1r": [],
            "uy1i": [],
            "uz1r": [],
            "uz1i": [],
            "ux2r": [],
            "ux2i": [],
            "uy2r": [],
            "uy2i": [],
            "uz2r": [],
            "uz2i": [],
            
            "n1r": [],
            "n1i": [],
            "n2r": [],
            "n2i": [],
            "ps1": [],
            "ps2": [],
            "p1ld": [],
            "p1ttd": [],
            "p1n0": [],
            "p1cd": [],
            "p2ld": [],
            "p2ttd": [],
            "p2n0": [],
            "p2cd": [],

            "p1tau": [],
            "p1mu": [],
            "p1alph": [],
            "p1q": [],
            "p1D": [],
            "p1vv": [],

            "p2tau": [],
            "p2mu": [],
            "p2alph": [],
            "p2q": [],
            "p2D": [],
            "p2vv": []
        }

        line = f.readline()
        while (line != ''):
            line = line.split()
            if(len(line) > 0):
                plume_sweep['kperp'].append(float(line[0]))
                plume_sweep['kpar'].append(float(line[1]))
                plume_sweep['betap'].append(float(line[2]))
                plume_sweep['vtp'].append(float(line[3]))
                plume_sweep['w'].append(float(line[4]))
                plume_sweep['g'].append(float(line[5]))
                plume_sweep['bxr'].append(float(line[6]))
                plume_sweep['bxi'].append(float(line[7]))
                plume_sweep['byr'].append(float(line[8]))
                plume_sweep['byi'].append(float(line[9]))
                plume_sweep['bzr'].append(float(line[10]))
                plume_sweep['bzi'].append(float(line[11]))
                plume_sweep['exr'].append(float(line[12]))
                plume_sweep['exi'].append(float(line[13]))
                plume_sweep['eyr'].append(float(line[14]))
                plume_sweep['eyi'].append(float(line[15]))
                plume_sweep['ezr'].append(float(line[16]))
                plume_sweep['ezi'].append(float(line[17]))
                plume_sweep['ux1r'].append(float(line[18]))
                plume_sweep['ux1i'].append(float(line[19]))
                plume_sweep['uy1r'].append(float(line[20]))
                plume_sweep['uy1i'].append(float(line[21]))
                plume_sweep['uz1r'].append(float(line[22]))
                plume_sweep['uz1i'].append(float(line[23]))
                plume_sweep['ux2r'].append(float(line[24]))
                plume_sweep['ux2i'].append(float(line[25]))
                plume_sweep['uy2r'].append(float(line[26]))
                plume_sweep['uy2i'].append(float(line[27]))
                plume_sweep['uz2r'].append(float(line[28]))
                plume_sweep['uz2i'].append(float(line[29]))
                
                plume_sweep['n1r'].append(float(line[30]))
                plume_sweep['n1i'].append(float(line[31]))
                plume_sweep['n2r'].append(float(line[32]))
                plume_sweep['n2i'].append(float(line[33]))
                plume_sweep['ps1'].append(float(line[34]))
                plume_sweep['ps2'].append(float(line[35]))
                plume_sweep['p1ld'].append(float(line[36]))
                plume_sweep['p1ttd'].append(float(line[37]))
                plume_sweep['p1n0'].append(float(line[38]))
                plume_sweep['p1cd'].append(float(line[39]))
                plume_sweep['p2ld'].append(float(line[40]))
                plume_sweep['p2ttd'].append(float(line[41]))
                plume_sweep['p2n0'].append(float(line[42]))
                plume_sweep['p2cd'].append(float(line[43]))

                plume_sweep['p1tau'].append(float(line[44]))
                plume_sweep['p1mu'].append(float(line[45]))
                plume_sweep['p1alph'].append(float(line[46]))
                plume_sweep['p1q'].append(float(line[47]))
                plume_sweep['p1D'].append(float(line[48]))
                plume_sweep['p1vv'].append(float(line[49]))

                plume_sweep['p2tau'].append(float(line[50]))
                plume_sweep['p2mu'].append(float(line[51]))
                plume_sweep['p2alph'].append(float(line[52]))
                plume_sweep['p2q'].append(float(line[53]))
                plume_sweep['p2D'].append(float(line[54]))
                plume_sweep['p2vv'].append(float(line[55]))
            line = f.readline()

        for key in plume_sweep.keys():
            plume_sweep[key] = np.asarray(plume_sweep[key])

        #normalize B
        plume_sweep['bxr'] = plume_sweep['bxr']*plume_sweep['vtp']
        plume_sweep['byr'] = plume_sweep['byr']*plume_sweep['vtp']
        plume_sweep['bzr'] = plume_sweep['bzr']*plume_sweep['vtp']
        plume_sweep['bxi'] = plume_sweep['bxi']*plume_sweep['vtp']
        plume_sweep['byi'] = plume_sweep['byi']*plume_sweep['vtp']
        plume_sweep['bzi'] = plume_sweep['bzi']*plume_sweep['vtp']
        
        return plume_sweep
        
    else:
        plume_sweep = {
            "kperp": [],
            "kpar": [],
            "betap": [],
            "vtp": [],
            "w": [],
            "g": [],
            "bxr": [],
            "bxi": [],
            "byr": [],
            "byi": [],
            "bzr": [],
            "bzi": [],
            "exr": [],
            "exi": [],
            "eyr": [],
            "eyi": [],
            "ezr": [],
            "ezi": [],
            "ux1r": [],
            "ux1i": [],
            "uy1r": [],
            "uy1i": [],
            "uz1r": [],
            "uz1i": [],
            "ux2r": [],
            "ux2i": [],
            "uy2r": [],
            "uy2i": [],
            "uz2r": [],
            "uz2i": [],
            
            "n1r": [],
            "n1i": [],
            "n2r": [],
            "n2i": [],
            "ps1": [], #power into species 1
            "ps2": [],
            "p1ttd1": [], #diagnol term in tensor; main transit time damping damping term
            "p1ttd2": [], #cross term in tensor; typically small
            "p1ld1": [], #diagnol term in tensor; main landau damping term
            "p1ld2": [], #cross term in tensor; typically small
            "p1n0": [],
            "p1cd": [],
            "p2ld1": [],
            "p2ld2": [],
            "p2ttd1": [],
            "p2ttd2": [],
            "p2n0": [],
            "p2cd": [],

            "p1tau": [],
            "p1mu": [],
            "p1alph": [],
            "p1q": [],
            "p1D": [],
            "p1vv": [],

            "p2tau": [],
            "p2mu": [],
            "p2alph": [],
            "p2q": [],
            "p2D": [],
            "p2vv": []
        }

        line = f.readline()

        while (line != ''):
            line = line.split()
            if(len(line) > 0):
                plume_sweep['kperp'].append(float(line[0]))
                plume_sweep['kpar'].append(float(line[1]))
                plume_sweep['betap'].append(float(line[2]))
                plume_sweep['vtp'].append(float(line[3]))
                plume_sweep['w'].append(float(line[4]))
                plume_sweep['g'].append(float(line[5]))
                plume_sweep['bxr'].append(float(line[6]))
                plume_sweep['bxi'].append(float(line[7]))
                plume_sweep['byr'].append(float(line[8]))
                plume_sweep['byi'].append(float(line[9]))
                plume_sweep['bzr'].append(float(line[10]))
                plume_sweep['bzi'].append(float(line[11]))
                plume_sweep['exr'].append(float(line[12]))
                plume_sweep['exi'].append(float(line[13]))
                plume_sweep['eyr'].append(float(line[14]))
                plume_sweep['eyi'].append(float(line[15]))
                plume_sweep['ezr'].append(float(line[16]))
                plume_sweep['ezi'].append(float(line[17]))
                plume_sweep['ux1r'].append(float(line[18]))
                plume_sweep['ux1i'].append(float(line[19]))
                plume_sweep['uy1r'].append(float(line[20]))
                plume_sweep['uy1i'].append(float(line[21]))
                plume_sweep['uz1r'].append(float(line[22]))
                plume_sweep['uz1i'].append(float(line[23]))
                plume_sweep['ux2r'].append(float(line[24]))
                plume_sweep['ux2i'].append(float(line[25]))
                plume_sweep['uy2r'].append(float(line[26]))
                plume_sweep['uy2i'].append(float(line[27]))
                plume_sweep['uz2r'].append(float(line[28]))
                plume_sweep['uz2i'].append(float(line[29]))
                
                plume_sweep['n1r'].append(float(line[30]))
                plume_sweep['n1i'].append(float(line[31]))
                plume_sweep['n2r'].append(float(line[32]))
                plume_sweep['n2i'].append(float(line[33]))
                plume_sweep['ps1'].append(float(line[34]))
                plume_sweep['ps2'].append(float(line[35]))
                plume_sweep['p1ttd1'].append(float(line[36]))
                plume_sweep['p1ttd2'].append(float(line[37]))
                plume_sweep['p1ld1'].append(float(line[38]))
                plume_sweep['p1ld2'].append(float(line[39]))
                plume_sweep['p1n0'].append(float(line[40]))
                plume_sweep['p1cd'].append(float(line[41]))
                plume_sweep['p2ttd1'].append(float(line[42]))
                plume_sweep['p2ttd2'].append(float(line[43]))
                plume_sweep['p2ld1'].append(float(line[44]))
                plume_sweep['p2ld2'].append(float(line[45]))
                plume_sweep['p2n0'].append(float(line[46]))
                plume_sweep['p2cd'].append(float(line[47]))

                plume_sweep['p1tau'].append(float(line[48]))
                plume_sweep['p1mu'].append(float(line[49]))
                plume_sweep['p1alph'].append(float(line[50]))
                plume_sweep['p1q'].append(float(line[51]))
                plume_sweep['p1D'].append(float(line[52]))
                plume_sweep['p1vv'].append(float(line[53]))

                plume_sweep['p2tau'].append(float(line[54]))
                plume_sweep['p2mu'].append(float(line[55]))
                plume_sweep['p2alph'].append(float(line[56]))
                plume_sweep['p2q'].append(float(line[57]))
                plume_sweep['p2D'].append(float(line[58]))
                plume_sweep['p2vv'].append(float(line[59]))
            line = f.readline()

        for key in plume_sweep.keys():
            plume_sweep[key] = np.asarray(plume_sweep[key])

        #normalize B
        plume_sweep['bxr'] = plume_sweep['bxr']*plume_sweep['vtp']
        plume_sweep['byr'] = plume_sweep['byr']*plume_sweep['vtp']
        plume_sweep['bzr'] = plume_sweep['bzr']*plume_sweep['vtp']
        plume_sweep['bxi'] = plume_sweep['bxi']*plume_sweep['vtp']
        plume_sweep['byi'] = plume_sweep['byi']*plume_sweep['vtp']
        plume_sweep['bzi'] = plume_sweep['bzi']*plume_sweep['vtp']

        return plume_sweep

def loadeigen(flnm):
    """
    Used to load plume eigenfunction output made by FPC routine for debugging.

    Parameters
    ----------
    flnm : string
        filename to be loaded

    Returns
    -------
    eigendict : dict
        dict containing eigenfunction data.
    """

    eigendict = {}

    f = open(flnm)
    line = f.readline()
    line = line.split()

    #TODO: this assumes nspec = 2. TODO: generalize for any num of nspec

    eigendict['kperp'] = line[0]
    eigendict['kpar'] = line[1]
    eigendict['betap'] = line[2]
    eigendict['vtp'] = line[3]
    eigendict['omega'] = complex(float(line[4]),float(line[5]))
    eigendict['bx'] = complex(float(line[6]),float(line[7]))
    eigendict['by'] = complex(float(line[8]),float(line[9]))
    eigendict['bz'] = complex(float(line[10]),float(line[11]))
    eigendict['ex'] = complex(float(line[12]),float(line[13]))
    eigendict['ey'] = complex(float(line[14]),float(line[15]))
    eigendict['ez'] = complex(float(line[16]),float(line[17]))
    eigendict['uxi'] = complex(float(line[18]),float(line[19]))
    eigendict['uyi'] = complex(float(line[20]),float(line[21]))
    eigendict['uzi'] = complex(float(line[22]),float(line[23]))
    eigendict['uxe'] = complex(float(line[24]),float(line[25]))
    eigendict['uye'] = complex(float(line[26]),float(line[27]))
    eigendict['uze'] = complex(float(line[28]),float(line[29]))
    eigendict['ni'] = complex(float(line[30]),float(line[31]))
    eigendict['ne'] = complex(float(line[32]),float(line[33]))
    eigendict['Pi'] = complex(float(line[34]),float(line[35]))
    eigendict['Pe'] = complex(0,0)#complex(float(line[36]),float(line[37]))

    return eigendict

def loadmoms(flnm):
    """
    Used to load jet-plume fs1-moment output made by FPC routine for debugging.

    Parameters
    ----------
    flnm : string
        filename to be loaded

    Returns
    -------
    momsdict : dict
        dict containing fs1 moments data.
    """

    momsdict = {}

    f = open(flnm)
    line = f.readline()
    line = line.split()

    momsdict['ni'] = complex(float(line[0]),float(line[1]))
    momsdict['ne'] = complex(float(line[2]),float(line[3]))
    momsdict['uxi'] = complex(float(line[4]),float(line[5]))
    momsdict['uyi'] = complex(float(line[6]),float(line[7]))
    momsdict['uzi'] = complex(float(line[8]),float(line[9]))
    momsdict['uxe'] = complex(float(line[10]),float(line[11]))
    momsdict['uye'] = complex(float(line[12]),float(line[13]))
    momsdict['uze'] = complex(float(line[14]),float(line[15]))
    momsdict['jxiex'] = float(line[16])
    momsdict['jyiey'] = float(line[17])
    momsdict['jziez'] = float(line[18])
    momsdict['jxeex'] = float(line[19])
    momsdict['jyeey'] = float(line[20])
    momsdict['jzeez'] = float(line[21])

    return momsdict

def branch_2var_scan_from_root(plume_input,stylenum1,stylenum2,var1key,var1min,var1max,var2key,var2min,var2max,root,inputflnm,outputname,swlog=False,outlog='outlog',nsamps=25,verbose=False,use_ps_split_new=True):
    """
    Makes four 2d sweeps that all start at point specified by params namelist and root and branches out in opposite directions.

    Please make sure that your point falls within the sweepmin and sweepmax range!

    Paramters
    ---------
    plume_input : class
        class containing input parameters
    sytlenum : int
        integer of sweep style (see main PLUME readme)
    sweepvarkey : string
        key of variable to be swept over
    var_n_min/var_n_max : float
        min and max value of variable to be swept over
    root : complex
        root to refine
    inputflnm : string
        name to call input file name
    outputname : string
        name to call output data name
    outlog : string
        filename to write plume output to (print statements- not data)
    verbose : bool
        if true, write print statements
    nsamps : int
        number of samples per sweep
    use_ps_split_new : bool
        use new or old power split (warning- must match .new_low_n. in vars.f90 (recompile if changed))
    swlog : bool
        use log spacing between sweep points

    Returns
    -------
    sweep : dict
        dict containing parallel arrays of sweep values
    """


    if(verbose):
        print("OVERWRITING OPTION AND NUM GUESS AND USE_MAP; TODO CHECK THAT plume_input IS CORRECT INSTEAD...")
    plume_input.params['option'] = 2
    plume_input.params['use_map'] = '.false.'
    plume_input.params['nroot_max'] = 1
    plume_input.params['nscan'] = 2

    sweepvar = -1
    midsweepval = 0.

    #Style -1 is not supported here.
    if(stylenum1 == 0):
        if(var1key == 'kperp'):
                sweepvar1 = '0'
                midsweepval1 = plume_input.params['kperp']
        elif(var1key == 'kpar'):
                sweepvar1 = '1'
                midsweepval1 = plume_input.params['kpar']
        elif(var1key == 'betap'):
                sweepvar1 = '2'
                midsweepval1 = plume_input.params['betap']
        elif(var1key == 'vtp'):
                sweepvar1 = '3'
                midsweepval1 = plume_input.params['vtp']
        else:
                print("********************************************************************************")
                print("Please input a valid var1key for this style!************************************")
                print("********************************************************************************")
    elif(stylenum1 >= 1):
        if(var2key == 'tauS'):
                sweepvar1 = '0'
                midsweepval1 = plume_input.species[stylenum1-1]['tauS']
        elif(var2key == 'muS'):
                sweepvar1 = '1'
                midsweepval1 = plume_input.species[stylenum1-1]['muS']
        elif(var2key == 'alphS'):
                sweepvar1 = '2'
                midsweepval1 = plume_input.species[stylenum1-1]['alphS']
        elif(var2key == 'Qs'):
                sweepvar1 = '3'
                midsweepval1 = plume_input.species[stylenum1-1]['Qs']
        elif(var2key == 'Ds'):
                sweepvar1 = '4'
                midsweepval1 = plume_input.species[stylenum1-1]['Ds']
        elif(var2key == 'vvS'):
                sweepvar1 = '5'
                midsweepval1 = plume_input.species[stylenum1-1]['vvS']
        else:
                print("********************************************************************************")
                print("Please input a valid var2key for this style!************************************")
                print("********************************************************************************")
    else:
        print("Please enter a valid stylenum1...")

    if(stylenum2 == 0):
        if(var2key == 'kperp'):
                sweepvar2 = '0'
                midsweepval2 = plume_input.params['kperp']
        elif(var2key == 'kpar'):
                sweepvar2 = '1'
                midsweepval2 = plume_input.params['kpar']
        elif(var2key == 'betap'):
                sweepvar2 = '2'
                midsweepval2 = plume_input.params['betap']
        elif(var2key == 'vtp'):
                sweepvar2 = '3'
                midsweepval2 = plume_input.params['vtp']
        else:
                print("********************************************************************************")
                print("Please input a valid var2key for this style!************************************")
                print("********************************************************************************")
    elif(stylenum2 >= 1):
        if(var2key == 'tauS'):
                sweepvar2 = '0'
                midsweepval2 = plume_input.species[stylenum1-1]['tauS']
        elif(var2key == 'muS'):
                sweepvar2 = '1'
                midsweepval2 = plume_input.species[stylenum1-1]['muS']
        elif(var2key == 'alphS'):
                sweepvar2 = '2'
                midsweepval2 = plume_input.species[stylenum1-1]['alphS']
        elif(var2key == 'Qs'):
                sweepvar2 = '3'
                midsweepval2 = plume_input.species[stylenum1-1]['Qs']
        elif(var2key == 'Ds'):
                sweepvar2 = '4'
                midsweepval2 = plume_input.species[stylenum1-1]['Ds']
        elif(var2key == 'vvS'):
                sweepvar2 = '5'
                midsweepval2 = plume_input.species[stylenum1-1]['vvS']
        else:
                print("********************************************************************************")
                print("Please input a valid var2key for this style!************************************")
                print("********************************************************************************")
    else:
        print("Please enter a valid stylenum1...")


    plume_input.guesses = [{}]
    g_om = root.real
    g_gam = root.imag
    plume_input.make_guess(g_om,g_gam)

    #generate quadrant 1 sweep
    outputnametemp1 = outputname+'sweep1'
    inputflnmtemp1 = inputflnm+'sweep1'
    scan_type=sweepvar1
    scan_style=stylenum1
    swi=midsweepval1
    swf=var1max
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)

    scan_type=sweepvar2
    scan_style=stylenum2
    swi=midsweepval2
    swf=var2max
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)

    plume_input.write_input(inputflnmtemp1,outputnametemp1,verbose=verbose)
    cmd = './plume.e ' + inputflnmtemp1 + '.in'
    cmd += ' >> ' + outlog
    if(verbose):
        print(cmd)
    os.system(cmd)

    #generate quadrant 2 sweep
    outputnametemp2 = outputname+'sweep2'
    inputflnmtemp2 = inputflnm+'sweep2'
    scan_type=sweepvar1
    scan_style=stylenum1
    swi=midsweepval1
    swf=var1min
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)

    scan_type=sweepvar2
    scan_style=stylenum2
    swi=midsweepval2
    swf=var2max
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)

    plume_input.write_input(inputflnmtemp2,outputnametemp2,verbose=verbose)
    cmd = './plume.e ' + inputflnmtemp2 + '.in'
    cmd += ' >> ' + outlog
    if(verbose):
        print(cmd)
    os.system(cmd)

    #generate quadrant 3 sweep
    outputnametemp3 = outputname+'sweep3'
    inputflnmtemp3 = inputflnm+'sweep3'
    scan_type=sweepvar1
    scan_style=stylenum1
    swi=midsweepval1
    swf=var1max
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)

    scan_type=sweepvar2
    scan_style=stylenum2
    swi=midsweepval2
    swf=var2min
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)

    plume_input.write_input(inputflnmtemp3,outputnametemp3,verbose=verbose)
    cmd = './plume.e ' + inputflnmtemp3 + '.in'
    cmd += ' >> ' + outlog
    if(verbose):
        print(cmd)
    os.system(cmd)

    #generate quadrant 4 sweep
    outputnametemp4 = outputname+'sweep4'
    inputflnmtemp4 = inputflnm+'sweep4'
    scan_type=sweepvar1
    scan_style=stylenum1
    swi=midsweepval1
    swf=var1min
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)

    scan_type=sweepvar2
    scan_style=stylenum2
    swi=midsweepval2
    swf=var2min
    swlog=swlog
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)

    plume_input.write_input(inputflnmtemp4,outputnametemp4,verbose=verbose)
    cmd = './plume.e ' + inputflnmtemp4 + '.in'
    cmd += ' >> ' + outlog
    if(verbose):
        print(cmd)
    os.system(cmd)

    #load sweeps
    flnmsweep1 = 'data/'+plume_input.dataname+'/'+outputnametemp1+'_'+var1key+'_'+var2key+'.mode1'
    if(verbose):
        print("Loading ",flnmsweep1,"...")
    sweep1 = load_plume_sweep(flnmsweep1,verbose=verbose,use_ps_split_new=use_ps_split_new)

    flnmsweep2 = 'data/'+plume_input.dataname+'/'+outputnametemp2+'_'+var1key+'_'+var2key+'.mode1'
    if(verbose):
        print("Loading ",flnmsweep2,"...")
    sweep2 = load_plume_sweep(flnmsweep2,verbose=verbose,use_ps_split_new=use_ps_split_new)

    flnmsweep3 = 'data/'+plume_input.dataname+'/'+outputnametemp3+'_'+var1key+'_'+var2key+'.mode1'
    if(verbose):
        print("Loading ",flnmsweep3,"...")
    sweep3 = load_plume_sweep(flnmsweep3,verbose=verbose,use_ps_split_new=use_ps_split_new)

    flnmsweep4 = 'data/'+plume_input.dataname+'/'+outputnametemp4+'_'+var1key+'_'+var2key+'.mode1'
    if(verbose):
        print("Loading ",flnmsweep4,"...")
    sweep4 = load_plume_sweep(flnmsweep4,verbose=verbose,use_ps_split_new=use_ps_split_new)

    if(verbose):
        print("Combining data and returning as 1 sweep...")
        print("Data has no particular order.")

    sweep = {}
    for _key in sweep1.keys():
        sweep[_key] = np.concatenate((sweep1[_key],sweep2[_key],sweep3[_key],sweep4[_key]))

    return sweep