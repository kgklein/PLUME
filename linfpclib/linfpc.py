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

    def set_fpc(self,vperpmin,vperpmax,vparmin,vparmax,delv,vxmin=None,vxmax=None,vymin=None,vymax=None,vzmin=None,vzmax=None):
        self.fpc = {'vperpmin':vperpmin,
                    'vperpmax':vperpmax,
                    'vparmin':vparmin,
                    'vparmax':vparmax,
                    'delv':delv
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

    def write_input(self,flnm,outputname,desc='',verbose=True):
        #TODO: option branches to check that all relevant info is there
        #calls write_namelist
        #TODO: write out date and time at top of file
        #TODO: write out description at top

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

    def load_from_file(self,flnm): #TODO: consider moving outside of class and have it make and return class
        paramkeys = ['betap','kperp','kpar','vtp','nspec','nscan','option','nroot_max','use_map','writeOut','dataName','outputName']
        fpckeys = ['vperpmin','vperpmax','vparmin','vparmax','delv']
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
                            tempparamdict[parse[0]] = int(parse[1].replace('\n',''))
                        else:
                            try:
                                tempparamdict[parse[0]] = float(parse[1].replace('\n',''))
                            except:
                                tempparamdict[parse[0]] = parse[1].replace('\n','')
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
                            tempfpcdict[parse[0]] = float(parse[1].replace('\n',''))
                        except:
                            tempfpcdict[parse[0]] = parse[1].replace('\n','')
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
                            tempmapdict[parse[0]] = float(parse[1].replace('\n',''))
                        except:
                            tempmapdict[parse[0]] = parse[1].replace('\n','')
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
                            tempspecdict[parse[0]] = float(parse[1].replace('\n',''))
                        except:
                            tempspecdict[parse[0]] = parse[1].replace('\n','')
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
                            tempscandict[parse[0]] = int(parse[1].replace('\n',''))
                        else:
                            try:
                                tempscandict[parse[0]] = float(parse[1].replace('\n',''))
                            except:
                                tempscandict[parse[0]] = parse[1].replace('\n','')
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
                            tempguessdict[parse[0]] = float(parse[1].replace('\n',''))
                        except:
                            tempguessdict[parse[0]] = parse[1].replace('\n','')
                    _i += 1
                if(len(self.guesses[0]) == 0):
                    self.guesses = [tempguessdict]
                else:
                    self.guesses.append(tempguessdict)


            _i += 1

        #TODO: handle dataName and outputName carefully


    #main dict
    namelists = {}
    params = {}
    fpc = {}
    species = [{}]
    maps = {}
    scan_inputs = [{}]
    guesses = [{}]

    dataname = 'default'

def _replace_input_aux(inputflnm,verbose=True):
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

#TODO: write output of plume somewhere automatically when python calls plume!!!! (instead of using '>> outlog')
def compute_roots(plume_input,inputflnm,outputname,outlog='outlog'):
    #use option 0

    print("OVERWRITING OPTION; TODO CHECK THAT plume_input IS CORRECT INSTEAD...")
    plume_input.params['option'] = 0
    plume_input.params['use_map'] = '.true.'

    #TODO: checks that input file is setup coorect and then call plume.e after making inputfile
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
    print("Reading roots from ",rootflnm)

    roots = []
    tempfile = open(rootflnm, "r")
    for line in tempfile:
        parse = line.split()
        temp_om = float(parse[4])
        temp_gam = float(parse[5])
        roots.append(temp_om+temp_gam*1j)
    tempfile.close()

    return np.asarray(roots)
    
def refine_root_and_calc_eigen_guess(plume_input,root,inputflnm,outputname,outlog='outlog',verbose=True):
    outputnametemp1 = outputname+'eigenatpoint'
    inputflnmtemp1 = inputflnm+'eigenatpoint'
    
    if(verbose):
        print("OVERWRITING OPTION; TODO CHECK THAT plume_input IS CORRECT INSTEAD...")   #TODO: or consider making a temporary copy
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
    
    #TODO: checks that input file is set up correctly and then call plume.e after making inputfile
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
    
    

def load_roots():
    pass

def compute_fpc_from_root(plume_input,root,inputflnm,outputname,outlog='outlog',cart=False):
    #use guess to compute fpc
    print("OVERWRITING OPTION AND NUM GUESS AND USE_MAP; TODO CHECK THAT plume_input IS CORRECT INSTEAD...")
    if(not(cart)):
        plume_input.params['option'] = 6
    else:
        plume_input.params['option'] = 7
    plume_input.params['use_map'] = '.false.'
    plume_input.params['nroot_max'] = 1

    print("Forcing input file to use root as initial guess")
    plume_input.guesses = [{}]
    g_om = root.real
    g_gam = root.imag
    plume_input.make_guess(g_om,g_gam)

    #TODO: checks that input file is set up correctly and then call plume.e after making inputfile
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
        cdatafilenames.append('data/'+str(plume_input.dataname)+'/'+outputname+'.specie0'+str(_nspec)+'.fs0') #warning: assumes num of species < 10
        cdatafilenames.append('data/'+str(plume_input.dataname)+'/'+outputname+'.specie0'+str(_nspec)+'.fs0') #warning: assumes num of species < 10

    return cdatafilenames

def make_sweeps_that_branch_from_params(plume_input,sweepvarkey,sweepmin,sweepmax,root,inputflnm,outputname,outlog='outlog',nsamps=200,verbose=True,use_ps_split_new=False):
    #makes sweeps that start at params namelist, and branches out

    if(verbose):
        print("OVERWRITING OPTION AND NUM GUESS AND USE_MAP; TODO CHECK THAT plume_input IS CORRECT INSTEAD...")
    plume_input.params['option'] = 1
    plume_input.params['use_map'] = '.false.'
    plume_input.params['nroot_max'] = 1

    sweepvar = -1
    midsweepval = 0.
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
            print("Please input a valid sweepvarkey!***********************************************")
            print("********************************************************************************")

    plume_input.guesses = [{}]
    g_om = root.real
    g_gam = root.imag
    plume_input.make_guess(g_om,g_gam)

    #do first sweep
    outputnametemp1 = outputname+'sweep1'
    inputflnmtemp1 = inputflnm+'sweep1'
    scan_type=sweepvar
    scan_style=0
    swi=midsweepval
    swf=sweepmax
    swlog=True
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)
    #TODO: checks that input file is set up correctly and then call plume.e after making inputfile
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
    scan_style=0
    swi=midsweepval
    swf=sweepmin
    swlog=True
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)
    #TODO: checks that input file is set up correctly and then call plume.e after making inputfile
    plume_input.write_input(inputflnmtemp2,outputnametemp2,verbose=verbose)
    cmd = './plume.e ' + inputflnmtemp2 + '.in'
    cmd += ' >> ' + outlog
    if(verbose):
        print(cmd)
    os.system(cmd)

    try:
        #load sweeps
        flnmsweep1 = 'data/'+plume_input.dataname+'/'+outputnametemp1+'_'+sweepvarkey+'_'+str(int(midsweepval*1000))+'_'+str(int(sweepmax*1000))+'.mode1'
        if(verbose):
            print("Loading ",flnmsweep1,"...")
        sweephigh = load_plume_sweep(flnmsweep1,verbose=verbose,use_ps_split_new=use_ps_split_new)

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


        def __combine_out_sweeps(sweep1,sweep2):
            pass

    except:
        #TODO: automatically handle this
        print("Could not automatically load and combine sweeps.")
        print("Probable cause of error is described below:")
        print("If sweepmax and sweepmin and `middle of sweep' val (variable `midsweepval' in function) are not multiples of .001, PLUME will round these numbers")
        print("and thus the output file name will be rounded as well. ")
        print("For now, it is recommended that one rounds sweepmax and sweepmin to multiples of .001 manually but this should be fixed...")

def loadlinfpcdist(filename,vparmin,vperpmin,vparmax,vperpmax,delv): #TODO: write these vel vars to file and automatically load

    print("Opening " + filename)
    try:
        f = open(filename)
    except:
        print("Couldnt open: " + filename)
        return

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
    dist = []
    line = f.readline()
    while (line != ''):
        line = line.split()
        row = []
        for k in line:
            if(math.isnan(float(k))):
                row.append(0.)
            else:
                row.append(float(k))
        dist.append(row)
        line = f.readline()

    dist=np.flip(dist,axis=1) #flip coordinates as they are backwards due to inconsistency of basis between routines

        #weird index rounding bug fix
    if(len(vpar) < len(dist[0])):
        vpar.append(float(vparindex))
    elif(len(vpar) > len(dist[0])):
        vpar.pop()
    if(len(vperp) < len(dist)):
        vperp.append(float(vperpindex))
    elif(len(vperp) > len(dist)):
        vperp.pop()

    linfpcdata = {'dist':dist,'vpar':vpar,'vperp':vperp,'vperpmin':vperpmin,'vperpmax':vperpmax,'vparmin':vparmin,'vparmax':vparmax,'delv':delv}

    return linfpcdata

def loadlinfpccepar(filename):

    print("Opening " + filename)
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
        print("Calculated resonant interval (elec): " + str(resonant_int))
        species = 'elec'
    else:
        print("Calculated resonant interval (ion): " + str(resonant_int))
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

    C = np.flip(C,axis=1) #flip coordinates as they are backwards due to inconsistency of basis between routines

        #weird index rounding bug fix
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

#TODO: make one load function
def loadlinfpccart(filename):
    print("Opening " + filename)
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
        print("Calculated resonant interval (elec): " + str(resonant_int))
        species = 'elec'
    else:
        print("Calculated resonant interval (ion): " + str(resonant_int))
        species = 'ion'

    resonant_int = resonant_int.real


    numstepvx = int(math.floor((vxmax-vxmin)/delv)) #TODO: implement this fix into all load functions
    numstepvy = int(math.floor((vymax-vymin)/delv))
    numstepvz = int(math.floor((vzmax-vzmin)/delv))

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
    while (linecounter <= len(vx)):
        line = line.split()
        row = []
        for k in line:
            if(math.isnan(float(k))):
                row.append(0.)
            else:
                row.append(float(k))
        Cvxvy.append(row)
        line = f.readline()
        linecounter += 1
    Cvxvy = np.flip(Cvxvy,axis=1) #We must adjust for coordinate transform that is a result of the subtly different coordinates used to calculate the eigenfunctions and the FPC related quantities

    Cvxvz = []
    linecounter = 1
    line = f.readline()
    while (linecounter <= len(vx)):
        line = line.split()
        row = []
        for k in line:
            if(math.isnan(float(k))):
                row.append(0.)
            else:
                row.append(float(k))
        Cvxvz.append(row)
        line = f.readline()
        linecounter += 1
    Cvxvz = np.asarray(Cvxvz)

    Cvyvz = []
    linecounter = 1
    line = f.readline()
    while (linecounter <= len(vy)):
        line = line.split()
        row = []
        for k in line:
            if(math.isnan(float(k))):
                row.append(0.)
            else:
                row.append(float(k))
        Cvyvz.append(row)
        line = f.readline()
        linecounter += 1
    Cvyvz = np.flip(Cvyvz,axis=0) #We must adjust for coordinate transform that is a result of the subtly different coordinates used to calculate the eigenfunctions and the FPC related quantities

    linfpcckeyname = 'CEpar'
    if('perp1' in filename):
        linfpcckeyname = 'CEperp1'
    elif('perp2' in filename):
        linfpcckeyname = 'CEperp2'

    linfpcdata = {linfpcckeyname+'vxvy':Cvxvy,linfpcckeyname+'vxvz':Cvxvz,linfpcckeyname+'vyvz':Cvyvz,'vx':vx,'vy':vy,'vz':vz,'resonant_int':resonant_int,'vxmin':vxmin,'vxmax':vxmax,'vymin':vymin,'vymax':vymax,'vzmin':vzmin,'vzmax':vzmax,'delv':delv,'species':species,'omega':omega}

    return linfpcdata

def loadlinfpcceperp(filename):
    import copy
    linfpcdata = loadlinfpccepar(filename)

    linfpcdata['CEperp'] = copy.deepcopy(linfpcdata['CEpar'])
    del linfpcdata['CEpar']

    return linfpcdata

def load_plume_sweep(flnm,verbose=True,use_ps_split_new=False): #TODO: check which ps split is used automatically or eliminate old one entirely

    """
    Load data from plume sweep

    Assumes 2 species

    Parameters
    ----------
    flnm : str
        path to sweep to be loaded

    Returns
    -------
    plume_sweep : dict
        dictionary of data related to plume
    """
    
    if(verbose):
        print("WARNING: assuming 2 species...\n TODO: write load_plume_sweep_nspec...")
        if(not(use_ps_split_new)):
            print("WARNING: assuming low_n is true (rather than new_low_n) (pass use_ps_split_new=True to change this)")
            print("If unsure, please check vars.f90 (change requires recompile)")
            print("If both are true, new_low_n is used")
        else:
            print("WARNING: assuming new_low_n is true (rather than low_n)")
            print("If unsure, please check vars.f90 (change requires recompile)")
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
        }

        line = f.readline()
        while (line != ''):
            line = line.split()
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
            "p1ld1": [], #diagnol term in tensor; main landau damping term
            "p1ld2": [], #cross term in tensor; typically small
            "p1ttd1": [], #diagnol term in tensor; main transit time damping damping term
            "p1ttd2": [], #cross term in tensor; typically small
            "p1n0": [],
            "p1cd": [],
            "p2ld1": [],
            "p2ld2": [],
            "p2ttd1": [],
            "p2ttd2": [],
            "p2n0": [],
            "p2cd": [],
        }

        line = f.readline()
        while (line != ''):
            line = line.split()
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
    
def double_k_scan_from_root(plume_input,ktotsweepmin,ktotsweepmax,root,inputflnm,outputname,outlog='outlog',nsamps=200):
    #makes sweep over kperp and kpar from k0=np.sqrt(kperp^2+kpar^2) to k1=np.sqrt(kperp^2+kpar^2)

    print("OVERWRITING OPTION AND NUM GUESS AND USE_MAP; TODO CHECK THAT plume_input IS CORRECT INSTEAD...")
    plume_input.params['option'] = 1
    plume_input.params['use_map'] = '.false.'
    plume_input.params['nroot_max'] = 1

    plume_input.guesses = [{}]
    g_om = root.real
    g_gam = root.imag
    plume_input.make_guess(g_om,g_gam)

    #do double sweep
    outputnametemp1 = outputname+'doublesweep1'
    inputflnmtemp1 = inputflnm+'doublesweep1'
    scan_type=0
    scan_style=-1
    swi=ktotsweepmin
    swf=ktotsweepmax
    swlog=True
    ns=nsamps
    nres=1
    heating=True
    eigen=True
    plume_input.scan_inputs = [{}]
    plume_input.make_scan(scan_type,scan_style,swi,swf,swlog,ns,nres,heating,eigen)
    #TODO: checks that input file is set up correctly and then call plume.e after making inputfile
    plume_input.write_input(inputflnmtemp1,outputnametemp1)
    cmd = './plume.e ' + inputflnmtemp1 + '.in'
    cmd += ' >> ' + outlog
    print(cmd)
    os.system(cmd)

    #load sweep
    flnmsweep1 = 'data/'+plume_input.dataname+'/'+outputnametemp1+'_k_'+str(int(plume_input.params['kperp']*1000))+'_'+str(int(plume_input.params['kpar']*1000))+'_'+str(int(ktotsweepmin*1000))+'_'+str(int(ktotsweepmax*1000))+'.mode1'
    print("Loading ",flnmsweep1,"...")
    dblsweep = load_plume_sweep(flnmsweep1)
    
    return dblsweep 
