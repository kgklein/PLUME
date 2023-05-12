import numpy as np
import matplotlib.pyplot as plt
import os

def plotlinfpcv2(linfpcdata,filename,zoomin=False,vlim=None,plotlog=False,setequal=False,xlim=None,ylim=None,plotresonant=False,clim=None,computeEner=True,extraText=None):
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

    if(plotlog):
        print('************************')
        print("WARNING: CHANGED C NORM AS TEMP FIX FOR LOG PLOT")
        print('************************')
        C = np.asarray(C)*100000


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
    print(plotlog)
    if(not(plotlog)):
        print("plotting linear scale!")
        plt.pcolormesh(vpar, vperp, C, vmax=rng,vmin=-rng, cmap="seismic", shading="gouraud")
    else:
        print("plotting log scale!")
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
    plt.xlabel('$v_{||}/v_{ts}$')
    plt.ylabel('$v_{\perp}/v_{ts}$')
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
    plt.savefig('figures/'+filename+'cmap.png',format='png',dpi=1000,facecolor='white', transparent=False)
    # plt.show() #TODO: option to save or show
    plt.close()

def plot_dist_func(linfpcdata,filename,plotlog=True): #TODO: stack multiple species
    from matplotlib import ticker, cm

    dist = linfpcdata['dist']
    vpar = linfpcdata['vpar']
    vperp = linfpcdata['vperp']
    vperpmin = linfpcdata['vperpmin']
    vperpmax = linfpcdata['vperpmax']
    vparmin = linfpcdata['vparmin']
    vparmax = linfpcdata['vparmax']
    delv = linfpcdata['delv']

    maxval = max(map(max, dist))
    minval = min(map(min, dist))
    rng = maxval

    #make results directory for figures
    if not os.path.exists('figures'):
        print("Making figures folder...")
        os.makedirs('figures')

    plt.figure()
    clev = np.arange(minval,rng+rng/100.,rng/50.) #makes contour 'smoother' (more levels)

    plt.set_cmap('plasma')
    if(plotlog):
        plt.contourf(vpar, vperp, dist, locator=ticker.LogLocator())
    else:
        plt.contourf(vpar, vperp, dist)
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim(vparmin,vparmax)
    plt.gca().set_ylim(vperpmin,vperpmax)
    print(vparmin,vparmax,vperpmin,vperpmax)
    # if(vlim != None):
    #     print("Changing xlim and ylim of plot!!!")
    #     plt.gca().set_xlim(-1.*vlim,vlim)
    #     plt.gca().set_ylim(0,vlim)

    # if(xlim != None):
    #     print("Changing xlim of plot!!!")
    #     plt.gca().set_xlim(xlim[0],xlim[1])
    plt.title('$f(v_{||},v_\perp)$')
    plt.xlabel('$v_{||}/v_{ts}$')
    plt.ylabel('$v_{\perp}/v_{ts}$')
    plt.colorbar()
    plt.grid()

    print("Saving figure to figures folder!")
    plt.savefig('figures/'+filename+'dist.png',format='png',dpi=1000,facecolor='white', transparent=False)
    plt.show()
    plt.close()

def plot_9pan_cart(foldername,flnm='',specnum='01',computeEner=False):
    from linfpclib.linfpc import loadlinfpccart

    print("Loading files...")
    print("Warning: assuming folder does not contain FPC data in cartesian coordinates for multiple predictions.")

    #TODO: this assumes that outputname was fpc when writing output. TODO: generalize this
    flnmread = foldername + 'fpc.cparcart.specie'+specnum+'.mode01'
    print('Reading: ',flnmread)
    cartpar = loadlinfpccart(flnmread)

    flnmread = foldername + 'fpc.cperp1.specie'+specnum+'.mode01'
    cartperp1 = loadlinfpccart(flnmread)

    flnmread = foldername + 'fpc.cperp2.specie'+specnum+'.mode01'
    cartperp2 = loadlinfpccart(flnmread)


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
            _tempim = axs[_j,_i].pcolormesh(cartcdata[dirkey[2:4]],cartcdata[dirkey[0:2]],np.asarray(cartcdata[ckey]).T[:,::-1],vmax=absmax,vmin=-absmax,cmap="seismic")
            axs[_j,_i].grid()
            axs[_j,_i].set_title(titlesprefixes[_j]+titlesuffixes[_i])
            axs[_j,_i].set_xlabel(xaxlabels[_i])
            axs[_j,_i].set_ylabel(yaxlabels[_i])

            fig.colorbar(_tempim, ax=axs[_j,_i])
            
            if(computeEner == True):
                delv = cartcdata[dirkey[0:2]][1]-cartcdata[dirkey[0:2]][0] #WARNING: ASSUMES SQUARE GRID
                ener = np.sum(np.asarray(cartcdata[ckey]).T[:,:])*delv*delv
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
        plt.close()
    else:
        plt.show()


#
