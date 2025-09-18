# linfpcplot.py>

import numpy as np
import matplotlib.pyplot as plt
import os

def plotlinfpc_gyro(linfpcdata,filename='',zoomin=False,vlim=None,plotlog=False,setequal=False,xlim=None,ylim=None,plotresonant=False,clim=None,computeEner=True,extraText=None):
    """
    Plots FPC velocity space signatures in gyrotropic space.

    Parameters
    ----------
    linfpcdata : dict
        dict returned by loadlinfpccepar/loadlinfpcceperp
    filename : string
        filename to save figure too
    """

    from matplotlib import ticker, cm
    import matplotlib.colors as colors

    if('CEpar' in linfpcdata.keys()):
        C = linfpcdata['CEpar']
        plttitle = r'$C_{E_{||}}(v_{||},v_\perp)$'
    elif('CEperp' in linfpcdata.keys()):
        C = linfpcdata['CEperp']
        plttitle = r'$C_{E_{\perp}}(v_{||},v_\perp)$'
    vpar = linfpcdata['vpar']
    vperp = linfpcdata['vperp']
    resonant_int = linfpcdata['resonant_int']
    vperpmin = linfpcdata['vperpmin']
    vperpmax = linfpcdata['vperpmax']
    vparmin = linfpcdata['vparmin']
    vparmax = linfpcdata['vparmax']
    delv = linfpcdata['delv']

    #restrict data to zoom interval if requested
    if(zoomin):
        vparplot = [] #temp array used to grab requested data
        cplot = [] #temp array used to grab requested data
        for a in range(0,len(C)):
            crow = []
            for b in range(0,len(C[a])):
                if((vpar[b] >= zoomcenter-zoomwidth and vpar[b] <= zoomcenter+zoomwidth)):
                    crow.append(C[a][b])
                    if(not(vpar[b] in vparplot)):
                        vparplot.append(vpar[b])
            cplot.append(crow)
        vpar = vparplot
        C = cplot
        del vparplot
        del cplot

    #force contour range to be centered around 0
    maxval = max(map(max, C))
    minval = min(map(min, C))
    rng = maxval
    if(abs(minval) > maxval):
        rng = abs(minval)

    if(clim != None):
        rng = abs(clim)

    #make results directory for figures
    if not os.path.exists('figures'):
        print("Making figures folder...")
        os.makedirs('figures')

    plt.figure()

    if(plotresonant):
        if(not(zoomin)):
            plt.plot([resonant_int,resonant_int],[vperpmin,vperpmax],c="black",lw=.30)
        elif(resonant_int > zoomcenter-zoomwidth and resonant_int < zoomcenter+zoomwidth):
            plt.plot([resonant_int,resonant_int],[vperpmin,vperpmax],c="black",lw=.30)
        else:
            print("Resonant interval occurs outside of requested plot range and will not be plotted.")

    plt.set_cmap('bwr')
    if(not(plotlog)):
        plt.pcolormesh(vpar, vperp, C, vmax=rng,vmin=-rng, cmap="seismic", shading="gouraud")
    else:
        plt.pcolormesh(vpar, vperp, C, cmap="seismic", shading="gouraud",norm=colors.SymLogNorm(linthresh=1., linscale=1., vmin=-rng, vmax=rng))

    if(vlim != None):
        print("Changing xlim and ylim of plot!!!")
        plt.gca().set_xlim(-1.*vlim,vlim)
        plt.gca().set_ylim(0,vlim)

    if(xlim != None):
        print("Changing xlim of plot!!!")
        plt.gca().set_xlim(xlim[0],xlim[1])
    if(ylim != None):
        print("Changing ylim of plot!!!")
        plt.gca().set_ylim(ylim[0],zlim[1])
    plt.title(plttitle)
    plt.xlabel(r'$v_{||}/v_{ts}$')
    plt.ylabel(r'$v_{\perp}/v_{ts}$')
    plt.colorbar()
    plt.grid()

    if(setequal):
        plt.gca().set_aspect('equal')

    if(computeEner == True):
        ener = np.sum(C)*delv*delv
        xscale = np.max(plt.gca().get_xlim())-np.min(plt.gca().get_xlim())
        xtext = np.min(plt.gca().get_xlim())+xscale*.1
        ytext = np.max(plt.gca().get_ylim())*.75

        if('CEpar' in linfpcdata.keys()):
            plt.text(xtext,ytext,r"$\int C_{E_{||}}(v_{||},v_\perp) \, \, d\mathbf{v}$ = "+"{:.2e}".format(ener))
        elif('CEperp' in linfpcdata.keys()):
            plt.text(xtext,ytext,r"$\int C_{E_{\perp}}(v_{||},v_\perp) \, \, d\mathbf{v}$ = "+"{:.2e}".format(ener))

    if(extraText != None):
        xtext = .5
        ytext = 2.15
        plt.text(xtext,ytext,extraText)

    print("Saving figure to figures folder!")
    if(filename != ''):
        plt.savefig('figures/'+filename+'cmap.png',format='png',dpi=1000,facecolor='white', transparent=False)
   
    plt.show()

def plotlinfpc_gyro_dist(linfpcdata,filename,plotkey,zoomin=False,vlim=None,plotlog=False,setequal=False,xlim=None,ylim=None,plotresonant=False,clim=None,extraText=None):
    """
    Plots FPC velocity space signatures in gyrotropic space.

    Parameters
    ----------
    linfpcdata : dict
        dict returned by loadlinfpccepar/loadlinfpcceperp
    filename : string
        filename to save figure too
    plotkey : 're_f' or 'im_f'
        picks if real or imag part is plotted
    """

    from matplotlib import ticker, cm
    import matplotlib.colors as colors

    if(plotkey == 'im_f'):
        C = linfpcdata['im_f']
        plttitle = r'Im{$f(v_{||},v_\perp)$}'
    elif(plotkey == 're_f'):
        C = linfpcdata['re_f']
        plttitle = r'Re{$f(v_{||},v_\perp)$}'
    vpar = linfpcdata['vpar']
    vperp = linfpcdata['vperp']
    resonant_int = linfpcdata['resonant_int']
    vperpmin = linfpcdata['vperpmin']
    vperpmax = linfpcdata['vperpmax']
    vparmin = linfpcdata['vparmin']
    vparmax = linfpcdata['vparmax']
    delv = linfpcdata['delv']

    #restrict data to zoom interval if requested
    if(zoomin):
        vparplot = [] #temp array used to grab requested data
        cplot = [] #temp array used to grab requested data
        for a in range(0,len(C)):
            crow = []
            for b in range(0,len(C[a])):
                if((vpar[b] >= zoomcenter-zoomwidth and vpar[b] <= zoomcenter+zoomwidth)):
                    crow.append(C[a][b])
                    if(not(vpar[b] in vparplot)):
                        vparplot.append(vpar[b])
            cplot.append(crow)
        vpar = vparplot
        C = cplot
        del vparplot
        del cplot

    #force contour range to be centered around 0
    maxval = max(map(max, C))
    minval = min(map(min, C))
    rng = maxval
    if(abs(minval) > maxval):
        rng = abs(minval)

    if(clim != None):
        rng = abs(clim)

    #make results directory for figures
    if not os.path.exists('figures'):
        print("Making figures folder...")
        os.makedirs('figures')

    plt.figure()

    if(plotresonant):
        if(not(zoomin)):
            print("Plotting interval...")
            plt.plot([resonant_int,resonant_int],[vperpmin,vperpmax],c="black",lw=.30)
        elif(resonant_int > zoomcenter-zoomwidth and resonant_int < zoomcenter+zoomwidth):
            print("Plotting interval...")
            plt.plot([resonant_int,resonant_int],[vperpmin,vperpmax],c="black",lw=.30)
        else:
            print("Resonant interval occurs outside of requested plot range and will not be plotted.")

    plt.set_cmap('bwr')
    if(not(plotlog)):
        print("plotting linear scale!")
        plt.pcolormesh(vpar, vperp, C, vmax=rng,vmin=-rng, cmap="PuOr", shading="gouraud")
    else:
        print("plotting log scale!")
        plt.pcolormesh(vpar, vperp, C, cmap="PuOr", shading="gouraud",norm=colors.SymLogNorm(linthresh=1., linscale=1., vmin=-rng, vmax=rng))

    if(vlim != None):
        print("Changing xlim and ylim of plot!!!")
        plt.gca().set_xlim(-1.*vlim,vlim)
        plt.gca().set_ylim(0,vlim)

    if(xlim != None):
        print("Changing xlim of plot!!!")
        plt.gca().set_xlim(xlim[0],xlim[1])
    if(ylim != None):
        print("Changing ylim of plot!!!")
        plt.gca().set_ylim(ylim[0],zlim[1])
    plt.title(plttitle)
    plt.xlabel(r'$v_{||}/v_{ts}$')
    plt.ylabel(r'$v_{\perp}/v_{ts}$')
    plt.colorbar()
    _cmap = plt.get_cmap("PuOr")
    zero_color = _cmap(0.5)   # mid-point of the colormap = value 0
    plt.gca().set_facecolor(zero_color)
    plt.grid()

    if(setequal):
        plt.gca().set_aspect('equal')

    if(extraText != None):
        xtext = .5
        ytext = 2.15
        plt.text(xtext,ytext,extraText)

    print("Saving figure to figures folder!")
    if(filename != ''):
        plt.savefig('figures/'+filename+'cmap.png',format='png',dpi=1000,facecolor='white', transparent=False)
    
    plt.show()

def plot_9pan_cart(foldername,filenametag,dataoverwritepar=None,dataoverwriteperp1=None,dataoverwriteperp2=None,flnm='',specnum='01',computeEner=False, scalevelocity=1):
    """
    Makes 3x3 plot of projections of FPC vel signature in cartesian coordinates

    Parameters
    ----------
    foldername : string
        location of data
    filenametag : string
        name of data
    flnm : string
        filename to save figure as
    specnum : '01' or '02' etc.
        species number to be plotted
    """

    from linfpclib.linfpc import loadlinfpccart

    if(dataoverwritepar is None):
        flnmread = foldername + filenametag +'.cparcart.specie'+specnum+'.mode01'
        print('Reading: ',flnmread)
        cartpar = loadlinfpccart(flnmread)

        flnmread = foldername + filenametag +'.cperp1.specie'+specnum+'.mode01'
        cartperp1 = loadlinfpccart(flnmread)

        flnmread = foldername + filenametag +'.cperp2.specie'+specnum+'.mode01'
        cartperp2 = loadlinfpccart(flnmread)
    else:
        cartpar = dataoverwritepar
        cartperp1 = dataoverwriteperp1
        cartperp2 = dataoverwriteperp2

    fig, axs = plt.subplots(3,3,figsize=(3*5,3*5),sharex=True)

    _hspace = .2
    _wspace = .2
    fig.subplots_adjust(hspace=_hspace,wspace=_wspace)

    dirkeys = ['vxvy','vxvz','vyvz']
    cartcdatas = [cartpar,cartperp1,cartperp2]
    ckeyprefixes = ['CEpar','CEperp1','CEperp2']
    titlesprefixes = [r'$C_{E_{||}}',r'$C_{E_{\perp,1}}',r'$C_{E_{\perp,2}}']
    titlesuffixes = [r'(v_{\perp,1},v_{\perp,2})$',r'(v_{\perp,1},v_{||})$',r'(v_{\perp,2},v_{||})$']
    xaxlabels = [r'$v_{\perp,1}$',r'$v_{\perp,1}$',r'$v_{\perp,2}$']
    yaxlabels = [r'$v_{\perp,2}$',r'$v_{||}$',r'$v_{||}$']
    _i = 0
    _j = 0
    for dirkey in dirkeys:
        for cartcdata in cartcdatas:
            ckey = ckeyprefixes[_j]+dirkey
            absmax = np.max(np.abs(cartcdata[ckey]))
            # if(dirkey == 'vyvz'):
            #     _tempim = axs[_j,_i].pcolormesh(np.asarray(cartcdata[dirkey[2:4]])*scalevelocity,np.asarray(cartcdata[dirkey[0:2]])*scalevelocity,np.asarray(cartcdata[ckey])[:,:],vmax=absmax,vmin=-absmax,cmap="seismic")
            # else:
            #     _tempim = axs[_j,_i].pcolormesh(np.asarray(cartcdata[dirkey[2:4]])*scalevelocity,np.asarray(cartcdata[dirkey[0:2]])*scalevelocity,np.asarray(cartcdata[ckey]).T[:,:],vmax=absmax,vmin=-absmax,cmap="seismic")
            _plot1 = cartcdata['v'+dirkey[1]+'_'+dirkey[1]+dirkey[3]]
            _plot2 = cartcdata['v'+dirkey[3]+'_'+dirkey[1]+dirkey[3]]
            _tempim = axs[_j,_i].pcolormesh(np.asarray(_plot1)*scalevelocity,np.asarray(_plot2)*scalevelocity,np.asarray(cartcdata[ckey])[:,:],vmax=absmax,vmin=-absmax,cmap="seismic")
            axs[_j,_i].grid()
            axs[_j,_i].set_title(titlesprefixes[_j]+titlesuffixes[_i])
            axs[_j,_i].set_xlabel(xaxlabels[_i])
            axs[_j,_i].set_ylabel(yaxlabels[_i])

            axs[_j,_i].axis('equal')

            fig.colorbar(_tempim, ax=axs[_j,_i])
            
            if(computeEner == True):
                delv = cartcdata[dirkey[0:2]][1]-cartcdata[dirkey[0:2]][0] #WARNING: ASSUMES SQUARE GRID
                ener = np.sum(np.asarray(cartcdata[ckey]).T[:,:])*delv*delv*delv
                xscale = np.max(axs[_j,_i].get_xlim())-np.min(axs[_j,_i].get_xlim())
                xtext = np.min(axs[_j,_i].get_xlim())+xscale*.1
                ytext = np.max(axs[_j,_i].get_ylim())*.75

                if('par' in ckey):
                    axs[_j,_i].text(xtext,ytext,r"$\int C_{E_{||}} \, \, d\mathbf{v}$ = "+"{:.2e}".format(ener))
                elif('perp1' in ckey):
                    axs[_j,_i].text(xtext,ytext,r"$\int C_{E_{\perp,1}}(v_{||},v_\perp) \, \, d\mathbf{v}$ = "+"{:.2e}".format(ener))
                elif('perp2' in ckey):
                    axs[_j,_i].text(xtext,ytext,r"$\int C_{E_{\perp,2}}(v_{||},v_\perp) \, \, d\mathbf{v}$ = "+"{:.2e}".format(ener))
            
            _j = _j + 1
        _j = 0
        _i = _i + 1

    if(flnm != ''):
        plt.savefig(flnm,format='png',dpi=300)

    plt.show()

def plot_fs1_re_im_cart(foldername,filenametag,dataoverwritedistfunccart=None,flnm='',specnum='01', scalevelocity=1):
    """
    Makes 2x3 plot of projections of fs1 in cartesian coordinates

    Parameters
    ----------
    foldername : string
        location of data
    filenametag : string
        name of data
    flnm : string
        filename to save figure as
    specnum : '01' or '02' etc.
        species number to be plotted
    """

    from linfpclib.linfpc import loadlinfpccart_dist

    if(dataoverwritedistfunccart is None):
        flnmreadimag = foldername + filenametag +'.dfs.imag.specie'+specnum+'.mode01'
        flnmreadreal = foldername + filenametag +'.dfs.real.specie'+specnum+'.mode01'

        distfunccart = loadlinfpccart_dist(flnmreadreal,flnmreadimag)

    else:
        distfunccart = dataoverwritedistfunccart

    fig, axs = plt.subplots(2,3,figsize=(3*5,2*5),sharex=True)

    _hspace = .2
    _wspace = .2
    fig.subplots_adjust(hspace=_hspace,wspace=_wspace)

    dirkeys = ['vxvy','vxvz','vyvz']
    cartcdatas = [distfunccart,distfunccart]
    ckeyprefixes = ['re_f','im_f']
    titlesprefixes = [r'Real{$f',r'Imag{$f']
    titlesuffixes = [r'(v_{\perp,1},v_{\perp,2})$}',r'(v_{\perp,1},v_{||})$}',r'(v_{\perp,2},v_{||})$}']
    xaxlabels = [r'$v_{\perp,1}$',r'$v_{\perp,1}$',r'$v_{\perp,2}$']
    yaxlabels = [r'$v_{\perp,2}$',r'$v_{||}$',r'$v_{||}$']
    _i = 0
    _j = 0
    for dirkey in dirkeys:
        for cartcdata in cartcdatas:
            ckey = ckeyprefixes[_j]+dirkey
            absmax = np.max(np.abs(cartcdata[ckey]))
            _plot1 = cartcdata['v'+dirkey[1]+'_'+dirkey[1]+dirkey[3]]
            _plot2 = cartcdata['v'+dirkey[3]+'_'+dirkey[1]+dirkey[3]]
            _tempim = axs[_j,_i].pcolormesh(np.asarray(_plot1)*scalevelocity,np.asarray(_plot2)*scalevelocity,np.asarray(cartcdata[ckey])[:,:],vmax=absmax,vmin=-absmax,cmap="PuOr")
            axs[_j,_i].grid()
            axs[_j,_i].set_title(titlesprefixes[_j]+titlesuffixes[_i])
            axs[_j,_i].set_xlabel(xaxlabels[_i])
            axs[_j,_i].set_ylabel(yaxlabels[_i])

            axs[_j,_i].axis('equal')

            _cmap = plt.get_cmap("PuOr")
            zero_color = _cmap(0.5)   # mid-point of the colormap = value 0
            axs[_j,_i].set_facecolor(zero_color)

            fig.colorbar(_tempim, ax=axs[_j,_i])
            
            _j = _j + 1
        _j = 0
        _i = _i + 1

    if(flnm != ''):
        plt.savefig(flnm,format='png',dpi=300)

    plt.show()


def plot_roots(roots,flnm='',xlim=[],ylim=[]):
    """
    Plots roots in complex freq space

    Parameters
    ----------
    roots : 1darray (complex)
        list of roots to plot
    flnm : string
        name of file to save figure to
    """

    oms = [rt.real for rt in roots]
    gams = [rt.imag for rt in roots]
    lbls = [_i for _i in range(0,len(oms))]

    plt.figure()
    plt.scatter(oms,gams)
    for _i, txt in enumerate(lbls):
        plt.gca().annotate(txt, (oms[_i], gams[_i]))
    plt.gca().set_aspect('equal')
    plt.grid()
    plt.xlabel(r"$\omega$") #normalized as in PLUME 
    plt.ylabel(r"$\gamma / \omega$")
    if(xlim != []):
        plt.xlim(xlim[0],xlim[1])
    if(ylim != []):
        plt.ylim(ylim[0],ylim[1])

    if(flnm != ''):
        plt.savefig(flnm,format='png',dpi=300)

    plt.show()

def plot_disp_rel(plumeinput, root, sweep, xkey, ykey, xlabel, ylabel, flnm ='', plot_root = True, xlim = [], ylim = [], semilogx = False):
    """
    Example plot of sweep to be used with 'make_sweeps_that_branch_from_params'

    Parameters
    ----------
    plumeinput : class
        plumeinput class
    root : complex
        start root of the sweep
    sweep : dict
        sweep data
    xkey : string
        sweep key that is to be plotted as x axis
    ykey : string
        sweep key that is to be plotted as y axis
    x/ylabel : string
        label to put on x/y axis
    flnm : string
        name of file to save figure to
    """
    plt.figure()
    xplot = sweep[xkey]
    yplot = sweep[ykey]

    #make positive
    if(not(semilogx) and sum(1 for y in yplot if y < 0) > len(yplot) / 2): #check if using full log plot and if there are more negative than positive values...
        yplot = -yplot
        ylabel = '-'+ylabel

    if(not(semilogx)):        
        plt.loglog(xplot,yplot)
    else:
        plt.semilogx(xplot,yplot)

    if(plot_root):
        plt.scatter(plumeinput.params[xkey],root.real)


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if(xlim != []):
        plt.xlim(xlim[0],ylim[0])
    if(ylim != []):
        plt.ylim(ylim[0],ylim[1])

    plt.grid()

    if(flnm != ''):
        plt.savefig(flnm,format='png',dpi=300)

    plt.show()

def plot_disp_power_2spec(sweep,flnm='',xlim=[],ylim=[]):
    """
    Quick example plot routine that plots power vs kperp

    Parameters
    ----------
    sweep : dict
        sweep data
    flnm : string
        filename of figure
    """

    #plot power rate due to different mechanisms/ susc tensor 'splits' 
    #INTO the particles in units of gamma/omega

    #note: here, a positive power means particles are gaining energy, 
    #i.e. the wave is losing energy
    #This is the opposite convention used above for the dispersion relation

    #note: as values can be positive and negative, 
    #we separate into positive and negative ararys to plot using a log scale
    #Here, we only plot positive values of power

    #p1 is species 1, ions. p2 is species 2, electrons

    ildposh = sweep['p1ld_zz']+sweep['p1ld_zy']
    ildnegh = -1*(sweep['p1ld_zz']+sweep['p1ld_zy'])
    eldposh = (sweep['p2ld_zz']+sweep['p2ld_zy'])
    eldnegh = -1*(sweep['p2ld_zz']+sweep['p2ld_zy'])

    ittdposh = sweep['p1ttd_yy']+sweep['p1ttd_yz']
    ittdnegh = -1*(sweep['p1ttd_yy']+sweep['p1ttd_yz'])
    ettdposh = (sweep['p2ttd_yy']+sweep['p2ttd_yz'])
    ettdnegh = -1*(sweep['p2ttd_yy']+sweep['p2ttd_yz'])

    itoth = sweep['p1']
    etoth = sweep['p2']

    gamma_over_omega = []
    for _i in range(0,len(sweep['kperp'])):
        gamma_over_omega.append(-sweep['g'][_i]/sweep['w'][_i])
    gamma_over_omega = np.asarray(gamma_over_omega)

    fig = plt.plot(figsize=(8,8))
    plt.loglog(sweep['kperp'],itoth,color='red',label=r'$\gamma_i/\omega$',ls='-',lw=2)
    plt.loglog(sweep['kperp'],etoth,color='blue',label=r'$\gamma_e/\omega$',ls='-',lw=2)
    plt.loglog(sweep['kperp'],ildposh,color='red',label=r'$\gamma_{i,ld}/\omega$',ls='-.',lw=2)
    plt.loglog(sweep['kperp'],eldposh,color='blue',label=r'$\gamma_{e,ld}/\omega$',ls='-.',lw=2)
    plt.loglog(sweep['kperp'],ittdposh,color='red',label=r'$\gamma_{i,ttd}/\omega$',ls='--',lw=2)
    plt.loglog(sweep['kperp'],ettdposh,color='blue',label=r'$\gamma_{e,ttd}/\omega$',ls='--',lw=2)
    plt.loglog(sweep['kperp'],gamma_over_omega,color='black',label=r'$\gamma/\omega$',lw=.75)

    plt.legend(fontsize=12)

    plt.tick_params(axis='both', which='both', direction='in',top=True, bottom=True, right=True, left=True)
    plt.grid()

    plt.ylabel(r'$-\gamma/\omega$')
    plt.xlabel(r'$k_\perp \rho_i$')

    if(ylim != []):
        plt.ylim(ylim[0],ylim[1])
    if(xlim != []):
        plt.xlim(xlim[0],xlim[1])

    if(flnm != ''):
        plt.savefig(flnm,format='png',dpi=300)

    plt.show()

def sweep2dplot(sweep2d,xkey,ykey,zkey,xlabel,ylabel,zlabel,flnm = '', xlim=[],ylim=[],vmin=None,vmax=None):
    """
    Plots 2d sweep

    Parameters
    ----------
    sweep2d : dict
        2d sweep data
    xkey : string
        sweep key that is to be plotted as x axis
    ykey : string
        sweep key that is to be plotted as y axis
    zkey : string
        sweep key that is to be plotted as z axis (i.e. color)
    x/y/zlabel : string
        label to put on x/y/z axis
    flnm : string
        name of file to save figure to
    """

    from matplotlib.tri import Triangulation
    from matplotlib.colors import LogNorm

    x = sweep2d[xkey]
    y = sweep2d[ykey]
    z = sweep2d[zkey]

    triang = Triangulation(x, y)
    plt.tripcolor(triang, z, cmap='viridis', edgecolors='none', vmin=vmin,vmax=vmax)
    cbar = plt.colorbar(label=zlabel, norm=LogNorm())

    if(xlim != []):
        plt.xlim(xlim[0],xlim[1])
    if(ylim != []):
        plt.ylim(ylim[0],ylim[1])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if(flnm != ''):
        plt.savefig(flnm,format='png',dpi=300)

    plt.show()

def compute_2d_from_3d_hist(out):
    #out is data from 3d load for hist/fs1 data

    out2df1ronly = lfpc.reduce_3d_to_projections(out['fs1_r'], out['vx'], out['vy'], out['vz'], 're_f')
    out2dfironly = lfpc.reduce_3d_to_projections(out['fs1_i'], out['vx'], out['vy'], out['vz'], 'im_f')
    out2d = out2df1ronly.copy()
    out2d['im_fvxvy'] = out2dfironly['im_fvxvy'].copy()
    out2d['im_fvxvz'] = out2dfironly['im_fvxvz'].copy()
    out2d['im_fvyvz'] = out2dfironly['im_fvyvz'].copy()
    
    return out2d
