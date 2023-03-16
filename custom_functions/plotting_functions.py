import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#from basic_functions import *
import basic_functions as bf
import extinction_fitting_functions as eff
quite_val = False

##################################MAKE_REVISED_PPXF_FIGURE##################################

#######################CUSTOM_CREATION_OF_COLORBAR_IN_SUBPLOTS_OF_MATPLOTLIB#######################
#This function is called for creating custom colorbar for subplots
#See https://stackoverflow.com/questions/23876588/matplotlib-colorbar-in-each-subplot for more details

def add_colorbar(mappable):
	last_axes = plt.gca()
	ax = mappable.axes
	fig = ax.figure
	divider = make_axes_locatable(ax)
	cax1 = divider.append_axes("right", size="5%", pad=0.15)
	cbar = plt.colorbar(mappable, cax=cax1)
	cbar.set_ticks(ticker.LogLocator(), update_ticks=True)
	cbar.ax.tick_params(size=0)
	return cbar

def add_colorbar_lin(mappable):
	last_axes = plt.gca()
	ax = mappable.axes
	fig = ax.figure
	divider = make_axes_locatable(ax)
	cax1 = divider.append_axes("right", size="5%", pad=0.15)
	cbar = plt.colorbar(mappable, cax=cax1)
	cbar.set_ticks(ticker.LinearLocator(), update_ticks=True)
	cbar.ax.tick_params(size=0)
	return cbar

#######################CUSTOM_CREATION_OF_COLORBAR_IN_SUBPLOTS_OF_MATPLOTLIB#######################



#----------------------------------------------------------------------------
def plot_kinemetry_profiles_velocity(k, fitcentre=False, name=None, figname='test.pdf'):
    """
    Based on the kinemetry results (passed in k), this routine plots radial
    profiles of the position angle (PA), flattening (Q), k1 and k5 terms.
    Last two plots are for X0,Y0 and systemic velocity
    """

    k0 = k.cf[:,0]
    k1 = np.sqrt(k.cf[:,1]**2 + k.cf[:,2]**2)
    k5 = np.sqrt(k.cf[:,5]**2 + k.cf[:,6]**2)
    k51 = k5/k1
    erk1 = (np.sqrt( (k.cf[:,1]*k.er_cf[:,1])**2 + (k.cf[:,2]*k.er_cf[:,2])**2 ))/k1
    erk5 = (np.sqrt( (k.cf[:,5]*k.er_cf[:,5])**2 + (k.cf[:,6]*k.er_cf[:,6])**2 ))/k5
    erk51 = ( np.sqrt( ((k5/k1) * erk1)**2 + erk5**2  ) )/k1
    

    fig,ax =plt.subplots(figsize=(7,8))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1,1,1])
    

    ax1 = plt.subplot(gs[0])
    ax1.errorbar(k.rad, k.pa, yerr=[k.er_pa, k.er_pa], fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
    ax1.set_ylabel('PA [deg]', fontweight='bold')
    if name:
        ax1.set_title(name, fontweight='bold')

    ax1.tick_params(axis='both', which='both', top=True, right=True)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.yaxis.set_tick_params(length=6)
    ax1.xaxis.set_tick_params(width=2)
    ax1.yaxis.set_tick_params(width=2)
    ax1.xaxis.set_tick_params(length=6)
    ax1.tick_params(which='minor', length=3)
    ax1.tick_params(which='minor', width=1)
    ax1.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)
    ax1.get_xaxis().set_ticklabels([])

    ax2 = plt.subplot(gs[1])
    ax2.errorbar(k.rad, k.q, yerr=[k.er_q, k.er_q], fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
    ax2.set_ylabel('Q ', fontweight='bold')
    #ax2.set_xlabel('R [arsces]')
    ax2.set_ylim(0,1)
    if fitcentre:
        ax2.set_title('Velocity, fit centre', fontweight='bold')
    else:
        ax2.set_title('Velocity, fixed centre', fontweight='bold')


    ax2.tick_params(axis='both', which='both', top=True, right=True)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax2.yaxis.set_tick_params(length=6)
    ax2.xaxis.set_tick_params(width=2)
    ax2.yaxis.set_tick_params(width=2)
    ax2.xaxis.set_tick_params(length=6)
    ax2.tick_params(which='minor', length=3)
    ax2.tick_params(which='minor', width=1)
    ax2.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax2.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax2.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax2.spines[axis].set_linewidth(2)
    ax2.get_xaxis().set_ticklabels([])


    
    ax3 = plt.subplot(gs[2])
    ax3.errorbar(k.rad, k1, yerr=[erk1, erk1], fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
    ax3.set_ylabel(r'$k_1$ [km/s]', fontweight='bold')

    ax3.tick_params(axis='both', which='both', top=True, right=True)
    ax3.tick_params(axis='both', which='major', labelsize=10)
    ax3.yaxis.set_tick_params(length=6)
    ax3.xaxis.set_tick_params(width=2)
    ax3.yaxis.set_tick_params(width=2)
    ax3.xaxis.set_tick_params(length=6)
    ax3.tick_params(which='minor', length=3)
    ax3.tick_params(which='minor', width=1)
    ax3.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax3.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax3.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax3.spines[axis].set_linewidth(2)
    ax3.get_xaxis().set_ticklabels([])


    ax4 = plt.subplot(gs[3])
    ax4.errorbar(k.rad, k51, yerr=[erk51, erk51], fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
    ax4.set_ylabel(r'$k_{51}$', fontweight='bold')
    ax4.set_xlabel('R [arsces]', fontweight='bold')
    
    ax4.tick_params(axis='both', which='both', top=True, right=True)
    ax4.tick_params(axis='both', which='major', labelsize=10)
    ax4.yaxis.set_tick_params(length=6)
    ax4.xaxis.set_tick_params(width=2)
    ax4.yaxis.set_tick_params(width=2)
    ax4.xaxis.set_tick_params(length=6)
    ax4.tick_params(which='minor', length=3)
    ax4.tick_params(which='minor', width=1)
    ax4.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax4.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax4.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax4.spines[axis].set_linewidth(2)
    ax4.get_xaxis().set_ticklabels([])
        
    ax5 = plt.subplot(gs[4])
    ax5.errorbar(k.rad, k.xc, yerr=k.er_xc, fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3, label='Xc')
    ax5.errorbar(k.rad, k.yc, yerr=k.er_yc, fmt='--o', mec='k', mew=1.2, color='salmon', mfc='salmon', capsize=3, label='Yc')
    ax5.set_ylabel(r'$X_c, Y_c$ [arsces]', fontweight='bold')
    ax5.set_xlabel('R [arsces]', fontweight='bold')
    ax5.legend()

    ax5.tick_params(axis='both', which='both', top=True, right=True)
    ax5.tick_params(axis='both', which='major', labelsize=10)
    ax5.yaxis.set_tick_params(length=6)
    ax5.xaxis.set_tick_params(width=2)
    ax5.yaxis.set_tick_params(width=2)
    ax5.xaxis.set_tick_params(length=6)
    ax5.tick_params(which='minor', length=3)
    ax5.tick_params(which='minor', width=1)
    ax5.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax5.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax5.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax5.spines[axis].set_linewidth(2)
    
    
    ax6 = plt.subplot(gs[5])
    ax6.errorbar(k.rad, k0, yerr=k.er_cf[:,0], fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
    ax6.hlines(np.median(k0),5,20,linestyles='dashed', colors='skyblue', label=r'median $V_{sys}$')
    ax6.set_ylabel(r'V$_{sys}$ [km/s]', fontweight='bold')
    ax6.set_xlabel('R [arsces]', fontweight='bold')
    ax6.legend()
    
    ax6.tick_params(axis='both', which='both', top=True, right=True)
    ax6.tick_params(axis='both', which='major', labelsize=10)
    ax6.yaxis.set_tick_params(length=6)
    ax6.xaxis.set_tick_params(width=2)
    ax6.yaxis.set_tick_params(width=2)
    ax6.xaxis.set_tick_params(length=6)
    ax6.tick_params(which='minor', length=3)
    ax6.tick_params(which='minor', width=1)
    ax6.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax6.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax6.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax6.spines[axis].set_linewidth(2)
        
    fig.tight_layout()
    #plt.show()
    plt.savefig(figname)


#----------------------------------------------------------------------------
def plot_kinemetry_profiles_sigma(k, fitcentre=False, name=None, photo=False, figname='test.pdf'):

    k0 = k.cf[:,0]
    er_k0 = k.er_cf[:,0]
    b2 = k.cf[:,4]
    er_b2 = k.er_cf[:,4]
    b4 = k.cf[:,8]
    er_b4 = k.er_cf[:,8]

    k20 = b2/k0
    er_k20 = np.sqrt( (er_b2/k0)**2 + (er_k0*(b2/k0**2))**2  )
    k40 = b4/k0
    er_k40 = np.sqrt( (er_b4/k0)**2 + (er_k0*(b4/k0**2))**2  )

    fig,ax =plt.subplots(figsize=(7,9))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1,1,1])

    ax1 = plt.subplot(gs[0])

    ax1.errorbar(k.rad, k.pa, yerr=[k.er_pa, k.er_pa], fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3, label='fit cen')
    ax1.set_ylabel('PA [deg]', fontweight='bold')
    if name:
        ax1.set_title(name, fontweight='bold')
    ax1.legend()
    if photo:
        ax1.set_xscale('log')

    ax1.tick_params(axis='both', which='both', top=True, right=True)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.yaxis.set_tick_params(length=6)
    ax1.xaxis.set_tick_params(width=2)
    ax1.yaxis.set_tick_params(width=2)
    ax1.xaxis.set_tick_params(length=6)
    ax1.tick_params(which='minor', length=3)
    ax1.tick_params(which='minor', width=1)
    ax1.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)
    ax1.get_xaxis().set_ticklabels([])


    ax2 = plt.subplot(gs[1])
    ax2.errorbar(k.rad, k.q, yerr=[k.er_q, k.er_q], fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
    ax2.set_ylabel('Q ', fontweight='bold')
    ax2.set_ylim(0,1)
    if fitcentre:
        if photo:
            ax2.set_title('photometry, fit centre', fontweight='bold')
            ax2.set_xscale('log')
        else:
            ax2.set_title(r'$\sigma$, fit centre', fontweight='bold')
    else:
        if photo:
            ax2.set_title('photometry, fixed centre', fontweight='bold')
            ax2.set_xscale('log')
        else:
            ax2.set_title(r'$\sigma$, fixed centre', fontweight='bold')



    ax2.tick_params(axis='both', which='both', top=True, right=True)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax2.yaxis.set_tick_params(length=6)
    ax2.xaxis.set_tick_params(width=2)
    ax2.yaxis.set_tick_params(width=2)
    ax2.xaxis.set_tick_params(length=6)
    ax2.tick_params(which='minor', length=3)
    ax2.tick_params(which='minor', width=1)
    ax2.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax2.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax2.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax2.spines[axis].set_linewidth(2)
    ax2.get_xaxis().set_ticklabels([])


    ax3 = plt.subplot(gs[2])
    if photo:
        ax3.errorbar(k.rad, k0, yerr=k.er_cf[:,0], fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
        ax3.set_ylabel(r'log$_{10} a_0$', fontweight='bold')
        ax3.set_xscale('log')
        ax3.set_yscale('log')
    else:
        ax3.errorbar(k.rad, k0, yerr=k.er_cf[:,0], fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
        ax3.set_ylabel(r'$\sigma_0$ [km/s]', fontweight='bold')

    ax3.tick_params(axis='both', which='both', top=True, right=True)
    ax3.tick_params(axis='both', which='major', labelsize=10)
    ax3.yaxis.set_tick_params(length=6)
    ax3.xaxis.set_tick_params(width=2)
    ax3.yaxis.set_tick_params(width=2)
    ax3.xaxis.set_tick_params(length=6)
    ax3.tick_params(which='minor', length=3)
    ax3.tick_params(which='minor', width=1)
    ax3.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax3.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax3.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax3.spines[axis].set_linewidth(2)
    ax3.get_xaxis().set_ticklabels([])
  
    ax4 = plt.subplot(gs[3])
    ax4.errorbar(k.rad, k20, yerr=er_k20, fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
    ax4.set_ylabel(r'$b_2/b_0$ [km/s]', fontweight='bold')
    if photo:
        ax4.set_xscale('log')

    ax4.tick_params(axis='both', which='both', top=True, right=True)
    ax4.tick_params(axis='both', which='major', labelsize=10)
    ax4.yaxis.set_tick_params(length=6)
    ax4.xaxis.set_tick_params(width=2)
    ax4.yaxis.set_tick_params(width=2)
    ax4.xaxis.set_tick_params(length=6)
    ax4.tick_params(which='minor', length=3)
    ax4.tick_params(which='minor', width=1)
    ax4.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax4.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax4.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax4.spines[axis].set_linewidth(2)
    ax4.get_xaxis().set_ticklabels([])


    ax5 = plt.subplot(gs[4])
    ax5.errorbar(k.rad, k40, yerr=er_k40, fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3)
    ax5.set_ylabel(r'$b_4/\sigma_0$', fontweight='bold')
    if photo:
        ax5.set_xscale('log')
        ax5.set_xlabel(r'log$_{10}$ R [arsces]', fontweight='bold')
    else:
        ax5.set_xlabel('R [arsces]', fontweight='bold')


    ax5.tick_params(axis='both', which='both', top=True, right=True)
    ax5.tick_params(axis='both', which='major', labelsize=10)
    ax5.yaxis.set_tick_params(length=6)
    ax5.xaxis.set_tick_params(width=2)
    ax5.yaxis.set_tick_params(width=2)
    ax5.xaxis.set_tick_params(length=6)
    ax5.tick_params(which='minor', length=3)
    ax5.tick_params(which='minor', width=1)
    ax5.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax5.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax5.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax5.spines[axis].set_linewidth(2)


    ax6 = plt.subplot(gs[5])
    ax6.errorbar(k.rad, k.xc, yerr=k.er_xc, fmt='--o', mec='k', mew=1.2, color='skyblue', mfc='skyblue', capsize=3, label='Xc')
    ax6.errorbar(k.rad, k.yc, yerr=k.er_yc, fmt='--o', mec='k', mew=1.2, color='salmon', mfc='salmon', capsize=3, label='Yc')
    ax6.set_ylabel('$X_c, Y_c$ [arsces]', fontweight='bold')
    ax6.legend()
    if photo:
        ax6.set_xscale('log')
        ax6.set_xlabel(r'log)$_{10}$ R [arsces]', fontweight='bold')
    else:
        ax6.set_xlabel('R [arsces]', fontweight='bold')
    
    ax6.tick_params(axis='both', which='both', top=True, right=True)
    ax6.tick_params(axis='both', which='major', labelsize=10)
    ax6.yaxis.set_tick_params(length=6)
    ax6.xaxis.set_tick_params(width=2)
    ax6.yaxis.set_tick_params(width=2)
    ax6.xaxis.set_tick_params(length=6)
    ax5.tick_params(which='minor', length=3)
    ax6.tick_params(which='minor', width=1)
    ax6.tick_params(axis='both', which='both', top=True, right=True)
    for tick in ax6.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax6.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax6.spines[axis].set_linewidth(2)
        
    fig.tight_layout()
    #plt.show()
    plt.savefig(figname)


#----------------------------------------------------------------------------
def plot_kinemetry_maps_rev(xbin, ybin, velbin, k, sigma=False, figname='test.pdf'):
    """
    Based on the kinemetry results (k) and original coordinates (xbin,ybin) and
    the analysed moment (i.e. velocity), this routine plots the original moment
    (i.e. velocity) map with overplotted best fitted ellispes, reconstructed
    (rotation) map and a map based on the full Fourier analysis.
    
    """

    
    k0 = k.cf[:,0]
    k1 = np.sqrt(k.cf[:,1]**2 + k.cf[:,2]**2)

    vsys=np.median(k0)
    if sigma:
        mx=np.max(k0)
        mn=np.min(k0)
        vsys=0
    else:
        mx=np.max(k1)
        mn=-mx
        vsys=np.median(k0)

    
    tmp=np.where(k.velcirc < 123456789)

    fig,ax =plt.subplots(figsize=(10,8))

    gs = gridspec.GridSpec(2, 2, width_ratios=[1,1], height_ratios=[1,1])

    ax1 = plt.subplot(gs[0, 0])
    im1 = plot_velfield(xbin, ybin, velbin - vsys, colorbar=False, label='km/s', nodots=True, vmin=mn, vmax=mx)
    ax1.plot(k.Xellip, k.Yellip, ',', label ='ellipse locations')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb=fig.colorbar(im1, cax=cax)
    cb.ax.tick_params(labelsize=12, width=2)
    cb.set_ticks([mn,0,mx],update_ticks=True)
    for axis in ['top','bottom','left','right']:
        cb.ax.spines[axis].set_linewidth(5)
    ax1.set_xlabel('arcsec', fontweight='bold')
    ax1.set_ylabel('arcsec', fontweight='bold')
    if sigma:
        ax1.set_title(r'$\sigma$', fontweight='bold')
        cb.set_label(r'$\sigma$ [km/s]', fontweight='bold')
    else:
        ax1.set_title('V', fontweight='bold')
        cb.set_label(r'V [km/s]', fontweight='bold')

    ax1.plot(k.xc, k.yc, '+', label='(Xc,Yc)')
    ax1.legend()

    ax1.tick_params(axis='both', which='both', top=True, right=True)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.yaxis.set_tick_params(length=6)
    ax1.xaxis.set_tick_params(width=2)
    ax1.yaxis.set_tick_params(width=2)
    ax1.xaxis.set_tick_params(length=6)
    ax1.tick_params(which='minor', length=3)
    ax1.tick_params(which='minor', width=1)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)



    ax1 = plt.subplot(gs[0, 1])
    im1 = plot_velfield(xbin[tmp], ybin[tmp], k.velcirc[tmp]- vsys, colorbar=False, label='km/s', nodots=True, vmin=mn, vmax=mx)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb=fig.colorbar(im1, cax=cax)
    cb.ax.tick_params(labelsize=12, width=2)
    cb.set_ticks([mn,0,mx],update_ticks=True)
    for axis in ['top','bottom','left','right']:
        cb.ax.spines[axis].set_linewidth(5)
    ax1.set_xlabel('arcsec', fontweight='bold')
    if sigma:
        ax1.set_title(r'$\sigma_0$', fontweight='bold')
        cb.set_label(r'$\sigma_0$ [km/s]', fontweight='bold')
    else:
        ax1.set_title(r'V$_{disk}$', fontweight='bold')
        cb.set_label(r'V$_{disk}$ [km/s]', fontweight='bold')

    ax1.plot(k.xc, k.yc, '+', label='(Xc,Yc)')
    ax1.legend()

    ax1.tick_params(axis='both', which='both', top=True, right=True)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.yaxis.set_tick_params(length=6)
    ax1.xaxis.set_tick_params(width=2)
    ax1.yaxis.set_tick_params(width=2)
    ax1.xaxis.set_tick_params(length=6)
    ax1.tick_params(which='minor', length=3)
    ax1.tick_params(which='minor', width=1)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)


    ax1 = plt.subplot(gs[1, 0])
    im1 = plot_velfield(xbin[tmp], ybin[tmp], k.velkin[tmp]-vsys, colorbar=False, label='km/s', nodots=True, vmin=mn, vmax=mx)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb=fig.colorbar(im1, cax=cax)
    cb.ax.tick_params(labelsize=12, width=2)
    cb.set_ticks([mn,0,mx],update_ticks=True)
    
    for axis in ['top','bottom','left','right']:
        cb.ax.spines[axis].set_linewidth(5)
    ax1.set_xlabel('arcsec', fontweight='bold')
    if sigma:
        ax1.set_title(r'$\sigma_{kin}$', fontweight='bold')
        cb.set_label(r'$\sigma$ [km/s]', fontweight='bold')
    else:
        ax1.set_title(r'V$_{kin}$', fontweight='bold')
        cb.set_label('V [km/s]', fontweight='bold')


    ax1.tick_params(axis='both', which='both', top=True, right=True)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.yaxis.set_tick_params(length=6)
    ax1.xaxis.set_tick_params(width=2)
    ax1.yaxis.set_tick_params(width=2)
    ax1.xaxis.set_tick_params(length=6)
    ax1.tick_params(which='minor', length=3)
    ax1.tick_params(which='minor', width=1)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)

    ax1 = plt.subplot(gs[1, 1])
    im1 = plot_velfield(xbin[tmp], ybin[tmp], k.velkin[tmp]-k.velcirc[tmp], colorbar=False, label='km/s', nodots=True, vmin=mn/10, vmax=mx/10)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb=fig.colorbar(im1, cax=cax)
    cb.ax.tick_params(labelsize=12, width=2)
    cb.set_ticks([mn/10,0,mx/10],update_ticks=True)
    
    for axis in ['top','bottom','left','right']:
        cb.ax.spines[axis].set_linewidth(5)
    ax1.set_xlabel('arcsec', fontweight='bold')
    if sigma:
        ax1.set_title(r'$\sigma_{residual}$', fontweight='bold')
        cb.set_label(r'$\sigma$ [km/s]', fontweight='bold')
    else:
        ax1.set_title(r'V$_{residual}$', fontweight='bold')
        cb.set_label('V [km/s]', fontweight='bold')


    ax1.tick_params(axis='both', which='both', top=True, right=True)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.yaxis.set_tick_params(length=6)
    ax1.xaxis.set_tick_params(width=2)
    ax1.yaxis.set_tick_params(width=2)
    ax1.xaxis.set_tick_params(length=6)
    ax1.tick_params(which='minor', length=3)
    ax1.tick_params(which='minor', width=1)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)

    fig.tight_layout()
    #plt.show()
    plt.savefig(figname)




def ppxf_figure(lam_gal, galaxy, noise, residuals, bestfit_solution_array, **kwargs):
	savefig = kwargs.get('savefig', 'save') # Save figure?
	good_index = kwargs.get('good_index', np.ones([len(lam_gal)], dtype=bool)) # Goodpixels in fit
	str_fit_val = kwargs.get('str_fit_val', str('Fit')) # str_fit_val
	quiet = kwargs.get('quiet', False) # quiet
	fig_x_size = kwargs.get('fig_x_size', 20) # fig_x_size
	fig_y_size = kwargs.get('fig_y_size', 16) # fig_y_size
	dpi_val = kwargs.get('dpi_val', 100) # dpi_val
	st_age_unique = kwargs.get('st_age_unique', np.array([0.1, 0.5, 1.0, 5.0, 10.0, 12.0])) # st_age_unique
	st_mass_unique = kwargs.get('st_mass_unique', np.zeros_like(st_age_unique)) # st_mass_unique
	st_lum_unique = kwargs.get('st_lum_unique', np.zeros_like(st_age_unique)) # st_lum_unique
	size_of_font = kwargs.get('size_of_font', 20) # size_of_font
	vel_window = kwargs.get('vel_window', 2000) # vel_window
	lick_index_species = kwargs.get('lick_index_species', np.array(['Hbeta_o', 'Mg1', 'NaD', 'Fe6189', 'Halpha'])) # lick_index_species
	mean_wave_list = kwargs.get('mean_wave_list', np.array([4864.8735, 5120.541667, 5898.791667, 6170.666667, 6540.166667])) # mean_wave_list
	em_index_species = kwargs.get('em_index_species', np.array(['Hbeta', '[OIII]4958', '[OIII]5007', '[NII]6547', 'Halpha', '[NII]6583', '[SII]6716', '[SII]6730'])) # em_index_species
	em_mean_wave_list = kwargs.get('em_mean_wave_list', np.array([4861.32, 4958.83, 5006.77, 6547.96, 6562.8, 6583.34, 6716.31, 6730.68])) # em_mean_wave_list
	fig_name_cust = kwargs.get('figname', 'test.pdf') # figname
	redshift = kwargs.get('redshift', 0.0) # redshift

	lam_gal = lam_gal / (1.+redshift)
	
	if not (savefig=='display'):
		fig = plt.figure(figsize=[fig_x_size, fig_y_size], dpi=dpi_val)
	else:
		fig = plt.figure()
		size_of_font = size_of_font/2

	if (np.any(st_mass_unique)!=0.):
		gs = gridspec.GridSpec(4, len(mean_wave_list), height_ratios=[len(mean_wave_list), 2, 2, 2])
	else:
		gs = gridspec.GridSpec(4, int(len(em_mean_wave_list)/2.), height_ratios=[int(len(em_mean_wave_list)/2.), 2, 2, 2])

	ax0 = plt.subplot(gs[0, :])
	ax1 = plt.subplot(gs[1, :])
	ax2 = []
	ax3 = []

	if (np.any(st_mass_unique)!=0.):
		ax2 = plt.subplot(gs[3, :])
		for i in range(len(mean_wave_list)):
			ax3 = np.append(ax3, plt.subplot(gs[2, i]))
	else:
		for i in range(0, int(len(em_mean_wave_list)/2.)):
			ax3 = np.append(ax3, plt.subplot(gs[2, i]))
		for j in range(0, int(len(em_mean_wave_list))-int(len(em_mean_wave_list)/2.)):
			ax2 = np.append(ax2, plt.subplot(gs[3, j]))

	ax0.plot(lam_gal, galaxy, color='tab:blue', drawstyle='steps-mid', alpha=0.4, label=r'data', zorder=1)
	ax0.plot(lam_gal[good_index], galaxy[good_index], color='tab:blue', drawstyle='steps-mid', alpha=0.8, label=r'goodpixel', zorder=2)
	ax0.plot(lam_gal, bestfit_solution_array, 'r--', label=str(str_fit_val), zorder=3)
	ax0.legend(loc='best', fontsize=size_of_font)
	ax0.set_xticks([])
	ax0.set_ylabel(r"Relative Flux", fontsize=size_of_font)
	ax1.plot(lam_gal, residuals, 'k.', markersize=size_of_font/4, alpha=0.2, zorder=4)
	ax1.plot(lam_gal[good_index], residuals[good_index], 'k.', markersize=size_of_font/4, alpha=0.8, zorder=5)
	ax1.axhline(0.0, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(-3, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(3, color='green', ls='dashed', alpha=0.5)
	ax1.text(6500, -4.0, s=r'-3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.text(6500, 3.2, s=r'3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.set_ylim(-5,5)
	ax1.set_xlabel(r"Rest Wavelength ($\rm \AA$)", fontsize=size_of_font)
	ax1.set_ylabel(r"Residuals", fontsize=size_of_font)

	if (np.any(st_mass_unique)!=0.):
		label_x_ticks = np.array(["<0.1", "0.1-0.5", "0.5-1.0", "1.0-5.0", "5.0-10.0", ">10.0"])
		ax2.plot(label_x_ticks, np.log10(st_mass_unique), 'bo', markersize=size_of_font, label='Stellar Mass')
		ax2.plot(label_x_ticks, np.log10(st_lum_unique), 'r*', markersize=size_of_font, label='Stellar Luminosity')
		ax2.set_xlabel(r"Age (Gyr)", fontsize=size_of_font)
		ax2.set_ylabel(r"Relative Weights (log)", fontsize=size_of_font)
		ax2.legend(loc='best', fontsize=size_of_font)
		ax2.tick_params(size=size_of_font)
		ax2.tick_params(axis='both', labelsize=size_of_font)

		for i in range(len(mean_wave_list)):
			vel31 = bf.vel_prof(lam_gal, mean_wave_list[i])
			idx_start = np.searchsorted(vel31, -vel_window)
			idx_end = np.searchsorted(vel31, vel_window)
			ax3[i].errorbar(vel31[idx_start:idx_end], galaxy[idx_start:idx_end], yerr=noise[idx_start:idx_end], ds='steps-mid', color='tab:blue')
			ax3[i].plot(vel31[idx_start:idx_end], bestfit_solution_array[idx_start:idx_end], 'r--')
			ax3[i].set_xlabel(r"Rel. Velocity (kms$^{-1}$)", fontsize=size_of_font)
			ax3[i].set_title(lick_index_species[i], fontsize=size_of_font)
			ax3[i].tick_params(axis='both', labelsize=size_of_font)
		ax3[0].set_ylabel(r"Rel. Flux", fontsize=size_of_font)

	else:
		for i in range(0, int(len(em_mean_wave_list)/2.)):
			vel31 = bf.vel_prof(lam_gal, em_mean_wave_list[i])
			idx_start = np.searchsorted(vel31, -vel_window)
			idx_end = np.searchsorted(vel31, vel_window)
			ax3[i].errorbar(vel31[idx_start:idx_end], galaxy[idx_start:idx_end], yerr=noise[idx_start:idx_end], ds='steps-mid', color='tab:blue')
			ax3[i].plot(vel31[idx_start:idx_end], bestfit_solution_array[idx_start:idx_end], 'r--')
			ax3[i].set_xlabel(r"Rel. Velocity (kms$^{-1}$)", fontsize=size_of_font)
			ax3[i].set_title(em_index_species[i], fontsize=size_of_font)
			ax3[i].tick_params(axis='both', labelsize=size_of_font)

		for j in range(int(len(em_mean_wave_list)/2.), int(len(em_mean_wave_list))):
			k = int(j - int(len(em_mean_wave_list)/2.))
			vel31 = bf.vel_prof(lam_gal, em_mean_wave_list[j])
			idx_start = np.searchsorted(vel31, -vel_window)
			idx_end = np.searchsorted(vel31, vel_window)
			ax2[k].errorbar(vel31[idx_start:idx_end], galaxy[idx_start:idx_end], yerr=noise[idx_start:idx_end], ds='steps-mid', color='tab:blue')
			ax2[k].plot(vel31[idx_start:idx_end], bestfit_solution_array[idx_start:idx_end], 'r--')
			ax2[k].set_xlabel(r"Rel. Velocity (kms$^{-1}$)", fontsize=size_of_font)
			ax2[k].set_title(em_index_species[j], fontsize=size_of_font)
			ax2[k].tick_params(axis='both', labelsize=size_of_font)

		ax3[0].set_ylabel(r"Rel. Flux", fontsize=size_of_font)
		ax2[0].set_ylabel(r"Rel. Flux", fontsize=size_of_font)

	ax0.tick_params(axis='both', labelsize=size_of_font)
	ax1.tick_params(axis='both', labelsize=size_of_font)
	#ax2.tick_params(axis='both', labelsize=size_of_font)
	ax0.tick_params(size=size_of_font)
	ax1.tick_params(size=size_of_font)
	#ax2.tick_params(size=size_of_font)
	#fig_name_cust = str("./all_pixel_images/Pos_") + str(int(x_pos)) + str("_") + str(int(y_pos)) + str("_fit_figure.pdf")
	if not (savefig=='display'):
		fig.tight_layout()
		plt.savefig(fig_name_cust)
		bf.print_cust(f'{fig_name_cust}, saved', quiet_val=quiet)
		plt.close()
	else:
		plt.show()
		





def make_figure_ppxf(fig_x_size, fig_y_size, dpi_val, lam_gal, galaxy, noise, residuals, bestfit_solution_array, st_age_unique, st_mass_unique, st_lum_unique, size_of_font, vel_window, mean_wave_list, lick_index_species, x_pos, y_pos, **kwargs):

	savefig = kwargs.get('savefig', 'display') # Save figure?
	good_index = kwargs.get('good_index', np.ones([len(lam_gal)], dtype=np.bool)) # Goodpixels in fit
	str_fit_val = kwargs.get('str_fit_val', str('Fit')) # Goodpixels in fit
	quiet = kwargs.get('quiet', False) # Goodpixels in fit

	if not (savefig=='display'):
		fig = plt.figure(figsize=[fig_x_size, fig_y_size], dpi=dpi_val)
	else:
		fig = plt.figure()
		size_of_font = size_of_font/2

	#gs = gridspec.GridSpec(4, 5, height_ratios=[5, 2, 2, 2])
	gs = gridspec.GridSpec(4, len(mean_wave_list), height_ratios=[5, 2, 2, 2])
	ax0 = plt.subplot(gs[0, :])
	ax1 = plt.subplot(gs[1, :])
	ax2 = plt.subplot(gs[3, :])
	ax3 = []
	for i in range(len(mean_wave_list)):
		ax3 = np.append(ax3, plt.subplot(gs[2, i]))

	#lam_gal = lam_gal / (1.+0.07527116015236746)
	ax0.plot(lam_gal, galaxy, color='tab:blue', drawstyle='steps-mid', alpha=0.4, label=r'data', zorder=1)
	ax0.plot(lam_gal[good_index], galaxy[good_index], color='tab:blue', drawstyle='steps-mid', alpha=0.8, label=r'goodpixel', zorder=2)
	ax0.plot(lam_gal, bestfit_solution_array, 'r--', label=str(str_fit_val), zorder=3)
	ax0.legend(loc='best', fontsize=size_of_font)
	ax0.set_xticks([])
	ax0.set_ylabel(r"Relative Flux", fontsize=size_of_font)
	ax1.plot(lam_gal, residuals, 'k.', markersize=size_of_font/4, alpha=0.2, zorder=4)
	ax1.plot(lam_gal[good_index], residuals[good_index], 'k.', markersize=size_of_font/4, alpha=0.8, zorder=5)
	ax1.axhline(0.0, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(-3, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(3, color='green', ls='dashed', alpha=0.5)
	ax1.text(6500, -4.0, s=r'-3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.text(6500, 3.2, s=r'3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.set_ylim(-5,5)
	ax1.set_xlabel(r"Rest Wavelength ($\rm \AA$)", fontsize=size_of_font)
	ax1.set_ylabel(r"Residuals", fontsize=size_of_font)
	label_x_ticks = np.array(["<0.1", "0.1-0.5", "0.5-1.0", "1.0-5.0", "5.0-10.0", ">10.0"])
	ax2.plot(label_x_ticks, np.log10(st_mass_unique), 'bo', markersize=size_of_font, label='Stellar Mass')
	ax2.plot(label_x_ticks, np.log10(st_lum_unique), 'r*', markersize=size_of_font, label='Stellar Luminosity')
	ax2.set_xlabel(r"Age (Gyr)", fontsize=size_of_font)
	ax2.set_ylabel(r"Relative Weights (log)", fontsize=size_of_font)
	ax2.legend(loc='best', fontsize=size_of_font)
	
	for i in range(len(mean_wave_list)):
		vel31 = bf.vel_prof(lam_gal, mean_wave_list[i])
		idx_start = np.searchsorted(vel31, -vel_window)
		idx_end = np.searchsorted(vel31, vel_window)
		ax3[i].errorbar(vel31[idx_start:idx_end], galaxy[idx_start:idx_end], yerr=noise[idx_start:idx_end], ds='steps-mid', color='tab:blue')
		ax3[i].plot(vel31[idx_start:idx_end], bestfit_solution_array[idx_start:idx_end], 'r--')
		ax3[i].set_xlabel(r"Rel. Velocity (kms$^{-1}$)", fontsize=size_of_font)
		ax3[i].set_title(lick_index_species[i], fontsize=size_of_font)
		ax3[i].tick_params(axis='both', labelsize=size_of_font)

	ax3[0].set_ylabel(r"Rel. Flux", fontsize=size_of_font)
	ax0.tick_params(axis='both', labelsize=size_of_font)
	ax1.tick_params(axis='both', labelsize=size_of_font)
	ax2.tick_params(axis='both', labelsize=size_of_font)
	ax0.tick_params(size=size_of_font)
	ax1.tick_params(size=size_of_font)
	ax2.tick_params(size=size_of_font)
	fig_name_cust = str("./all_pixel_images/Pos_") + str(int(x_pos)) + str("_") + str(int(y_pos)) + str("_fit_figure.pdf")
	if not (savefig=='display'):
		fig.tight_layout()
		plt.savefig(fig_name_cust)
		bf.print_cust(f'{fig_name_cust}, saved', quiet_val=quiet)
		plt.close()
	else:
		plt.show()
##################################MAKE_REVISED_PPXF_FIGURE##################################


##################################MAKE_FIGURE##################################
def make_figure(x_slot, y_slot, fig_x_size, fig_y_size, dpi_val, size_of_font, qso_name, fig_name, name_list_sorted, center_list_sorted, wave_fit, flux_fit, err_fit, cont_array_fitted, result_fitted, reddening_array_fit, redshift_val, vel_window, center_list_init, comments_on_balmer, number_of_narrow_components_init, number_of_wide_components_init, amp_array_rev, center_array, sigma_array, quiet_val=False):
	f, axarr = plt.subplots(x_slot, y_slot, sharex=False, sharey=False, figsize=((fig_x_size), (fig_y_size)), dpi=dpi_val)
	f.subplots_adjust(wspace=0.4, hspace=0.4)
	for u in range(x_slot):
		for v in range(y_slot):
			if (len(center_list_sorted) > ((x_slot*v)+u)):
				test_x3 = bf.vel_prof(wave_fit/(1.+redshift_val), center_list_sorted[(x_slot*v)+u])
				idx1 = bf.find_nearest_idx(test_x3,-vel_window-100)
				idx2 = bf.find_nearest_idx(test_x3,vel_window+100)
				if (u>1):
					axarr[u,v].errorbar(wave_fit[idx1:idx2], flux_fit[idx1:idx2], yerr=err_fit[idx1:idx2], color='black', drawstyle='steps-mid', zorder=1)
					axarr[u,v].plot(wave_fit[idx1:idx2], result_fitted[idx1:idx2], color='red', ls='--', zorder=2)
				else:
					axarr[v].errorbar(wave_fit[idx1:idx2], flux_fit[idx1:idx2], yerr=err_fit[idx1:idx2], color='black', drawstyle='steps-mid', zorder=1)
					axarr[v].plot(wave_fit[idx1:idx2], result_fitted[idx1:idx2], color='red', ls='--', zorder=2)
				count_amp = 0
				for z in range(len(center_list_init)):
					centre = center_list_init[z]
					vel_array = bf.vel_prof(wave_fit/(1.+redshift_val), centre)
					count = 0
					if (comments_on_balmer[z]):
						for k in range(number_of_narrow_components_init):
							#group_prof += bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
							prof_fit = bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count]) + cont_array_fitted
							prof_fit_adv = eff.func_6(wave_fit/(1.+redshift_val), prof_fit, *reddening_array_fit)
							if (u>1):
								axarr[u,v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'g--')
								axarr[u,v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='green', alpha=0.5)
							else:
								axarr[v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'g--')
								axarr[v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='green', alpha=0.5)

							prof_fit = []
							prof_fit_adv = []
							count+=1
							count_amp+=1
						for l in range(number_of_wide_components_init):
							#group_prof += bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
							prof_fit = bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count]) + cont_array_fitted
							prof_fit_adv = eff.func_6(wave_fit/(1.+redshift_val), prof_fit, *reddening_array_fit)
							if (u>1):
								axarr[u,v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'm--')
								axarr[u,v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='magenta', alpha=0.5)
							else:
								axarr[v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'm--')
								axarr[v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='magenta', alpha=0.5)
							prof_fit = []
							prof_fit_adv = []
							count+=1
							count_amp+=1
					else:
						for m in range(number_of_narrow_components_init):
							#group_prof += bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
							prof_fit = bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count]) + cont_array_fitted
							prof_fit_adv = eff.func_6(wave_fit/(1.+redshift_val), prof_fit, *reddening_array_fit)
							if (u>1):
								axarr[u,v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'g--')
								axarr[u,v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='green', alpha=0.5)
							else:
								axarr[v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'g--')
								axarr[v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='green', alpha=0.5)
							prof_fit = []
							prof_fit_adv = []
							count+=1
							count_amp+=1
				if (u>1):
					axarr[u,v].set_title(name_list_sorted[(x_slot*v)+u], fontsize=size_of_font)
					axarr[u,v].set_xlim(wave_fit[idx1], wave_fit[idx2])
					axarr[u,v].tick_params(axis = 'both', which = 'major', direction='in', length=size_of_font/2, width=2, colors='k')
					axarr[u,v].tick_params(axis = 'both', which = 'minor', direction='in', length=size_of_font/4, width=1, colors='k')
					axarr[u,v].tick_params(axis='both', labelsize=size_of_font)
				else:
					axarr[v].set_title(name_list_sorted[(x_slot*v)+u], fontsize=size_of_font)
					axarr[v].set_xlim(wave_fit[idx1], wave_fit[idx2])
					axarr[v].tick_params(axis = 'both', which = 'major', direction='in', length=size_of_font/2, width=2, colors='k')
					axarr[v].tick_params(axis = 'both', which = 'minor', direction='in', length=size_of_font/4, width=1, colors='k')
					axarr[v].tick_params(axis='both', labelsize=size_of_font)
			else:
				if (u>1):
					axarr[u,v].set_visible(False)
				else:
					axarr[v].set_visible(False)
	f.text(0.48, 0.94, qso_name, fontsize=1.2*size_of_font)
	f.text(0.52, 0.06, r'Wavelength ($\rm \AA$)', ha='center', va='center', fontsize=1.2*size_of_font)
	f.text(0.06, 0.5, r'Flux (10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ AA$^{-1}$)', ha='center', va='center', rotation='vertical', fontsize=1.2*size_of_font)
	f.patch.set_linewidth(2)
	f.patch.set_edgecolor('black')
	plt.savefig(fig_name, edgecolor=f.get_edgecolor())
	bf.print_cust(f'{fig_name} saved...', quiet_val=quiet_val)

def plot_quick_figure(wave, flux, err, wave_fitted, flux_fitted, flux_err_fitted, result_fitted, continuum_fitted, center_list_init, redshift_val, comments_on_balmer, center_array, number_of_narrow_components_init, number_of_wide_components_init):
	fig_secondary_1d_data, (ax_secondary_1d_data) = plt.subplots()
	ax_secondary_1d_data.errorbar(wave, flux, yerr=err, color='tab:blue', alpha=0.5, label='Data', zorder=1)
	ax_secondary_1d_data.errorbar(wave_fitted, flux_fitted, yerr=flux_err_fitted, color='tab:blue', label='Fitted Region', zorder=2)
	ax_secondary_1d_data.plot(wave_fitted, result_fitted, 'r.-', label='Total Fit', zorder=4)
	ax_secondary_1d_data.plot(wave_fitted, continuum_fitted, 'g--', label='Continuum', zorder=6)
	count_amp = 0
	for j in range(len(center_list_init)):
		centre = center_list_init[j]*(1.+redshift_val)
		vel_array = bf.vel_prof(wave_fitted, centre)
		count = 0
		if (comments_on_balmer[j]):
			for k in range(number_of_narrow_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='green', alpha=0.5)
				count+=1
				count_amp+=1
			for l in range(number_of_wide_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='magenta', alpha=0.5)
				count+=1
				count_amp+=1
		else:
			for m in range(number_of_narrow_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='green', alpha=0.5)
				count+=1
				count_amp+=1
	ax_secondary_1d_data.set_xlabel(r"Wavelength ($\rm \AA$)")
	ax_secondary_1d_data.set_ylabel(r"Relative Flux")
	plt.legend()
	plt.show()

def save_quick_figure_rev(wave, flux, err, wave_fitted, flux_fitted, flux_err_fitted, result_fitted, continuum_fitted, residuals, center_list_init, redshift_val, comments_on_balmer, center_array, number_of_narrow_components_init, number_of_wide_components_init, fig_x_size, fig_y_size, dpi_val, x_pos, y_pos, size_of_font, str_fit, mean_wave_list, lick_index_species, vel_window, save_fig='display', **kwargs):
	quite_val = kwargs.get('quite_val', False)  # length of amplitude array to be fitted
	if not (save_fig=='display'):
		fig = plt.figure(figsize=[fig_x_size, fig_y_size], dpi=dpi_val)
	else:
		fig = plt.figure()

	gs = gridspec.GridSpec(3, len(mean_wave_list), height_ratios=[5, 2, 2])
	ax3 = []
	for i in range(len(mean_wave_list)):
		ax3 = np.append(ax3, plt.subplot(gs[2, i]))
	ax0 = plt.subplot(gs[0, :])
	ax1 = plt.subplot(gs[1, :], sharex=ax0)
	lam_gal = wave_fitted / (1.+redshift_val)
	wave = wave / (1.+redshift_val)
	galaxy = flux_fitted
	noise = flux_err_fitted
	bestfit_solution_array = result_fitted
	ax0.errorbar(wave, flux, yerr=err, color='tab:blue', drawstyle='steps-mid', alpha=0.5, label=r'data', zorder=1)
	ax0.plot(lam_gal, flux_fitted, color='tab:blue', drawstyle='steps-mid', alpha=0.8, label=r'to fit', zorder=2)
	ax0.plot(lam_gal, result_fitted, 'r--', label=str(str_fit), zorder=3)
	ax0.legend(loc='best', fontsize=size_of_font)
	ax0.set_ylabel(r"Relative Flux", fontsize=size_of_font)
	ax0.set_xlim(np.nanmin(lam_gal)-10., np.nanmax(lam_gal)+10.)
	ax0.set_ylim(np.nanmin(flux_fitted)-5., np.nanmax(flux_fitted)+5.)
	ax1.plot(lam_gal, residuals, 'k.', markersize=size_of_font/4, alpha=0.8, zorder=4)
	ax1.axhline(0.0, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(-3, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(3, color='green', ls='dashed', alpha=0.5)
	ax1.text(6500, -4.0, s=r'-3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.text(6500, 3.2, s=r'3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.set_ylim(-5,5)
	ax1.set_xlabel(r"Rest Wavelength ($\rm \AA$)", fontsize=size_of_font)
	ax1.set_ylabel(r"Residuals", fontsize=size_of_font)
	for i in range(len(mean_wave_list)):
		vel31 = bf.vel_prof(lam_gal, mean_wave_list[i])
		idx_start = np.searchsorted(vel31, -vel_window)
		idx_end = np.searchsorted(vel31, vel_window)
		ax3[i].errorbar(vel31[idx_start:idx_end], galaxy[idx_start:idx_end], yerr=noise[idx_start:idx_end], ds='steps-mid', color='tab:blue')
		ax3[i].plot(vel31[idx_start:idx_end], bestfit_solution_array[idx_start:idx_end], 'r--')
		ax3[i].set_xlabel(r"Rel. Velocity (kms$^{-1}$)", fontsize=size_of_font)
		ax3[i].set_title(lick_index_species[i], fontsize=size_of_font)
		ax3[i].tick_params(axis='both', labelsize=size_of_font)

	ax0.tick_params(axis='both', labelsize=size_of_font)
	ax1.tick_params(axis='both', labelsize=size_of_font)
	ax3[0].set_ylabel(r"Rel. Flux", fontsize=size_of_font)
	ax0.tick_params(size=size_of_font)
	ax1.tick_params(size=size_of_font)
	fig_name_extension = kwargs.get('fig_name_extension', '')  # length of amplitude array to be fitted
	fig_name_cust_orig = str("./all_pixel_images/Pos_") + str(int(x_pos)) + str("_") + str(int(y_pos)) + str("_fit_figure") + str(fig_name_extension) + str(".pdf")
	if not (save_fig=='display'):
		fig.tight_layout()
		plt.savefig(fig_name_cust_orig)
		bf.print_cust(f'{fig_name_cust_orig} saved', quiet_val=quite_val)
	else:
		plt.show()

	plt.close('all')

def save_quick_figure(wave, flux, err, wave_fitted, flux_fitted, flux_err_fitted, result_fitted, continuum_fitted, center_list_init, redshift_val, comments_on_balmer, center_array, number_of_narrow_components_init, number_of_wide_components_init, fig_x_size, fig_y_size, dpi_val, x_pos, y_pos, save_fig='display', quite_val=False):
	if not (save_fig=='display'):
		fig_secondary_1d_data, (ax_secondary_1d_data) = plt.subplots(figsize=[fig_x_size, fig_y_size], dpi=dpi_val)
	else:
		fig_secondary_1d_data, (ax_secondary_1d_data) = plt.subplots()

	ax_secondary_1d_data.errorbar(wave, flux, yerr=err, color='tab:blue', alpha=0.5, label='Data', zorder=1)
	ax_secondary_1d_data.errorbar(wave_fitted, flux_fitted, yerr=flux_err_fitted, color='tab:blue', label='Fitted Region', zorder=2)
	ax_secondary_1d_data.plot(wave_fitted, result_fitted, 'r.-', label='Total Fit', zorder=4)
	ax_secondary_1d_data.plot(wave_fitted, continuum_fitted, 'g--', label='Continuum', zorder=6)
	count_amp = 0
	for j in range(len(center_list_init)):
		centre = center_list_init[j]*(1.+redshift_val)
		vel_array = bf.vel_prof(wave_fitted, centre)
		count = 0
		if (comments_on_balmer[j]):
			for k in range(number_of_narrow_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='green', alpha=0.5)
				count+=1
				count_amp+=1
			for l in range(number_of_wide_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='magenta', alpha=0.5)
				count+=1
				count_amp+=1
		else:
			for m in range(number_of_narrow_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='green', alpha=0.5)
				count+=1
				count_amp+=1
	ax_secondary_1d_data.set_xlabel(r"Wavelength ($\rm \AA$)")
	ax_secondary_1d_data.set_ylabel(r"Relative Flux")
	fig_secondary_1d_data.legend()
	fig_name_cust = str("./all_pixel_images/Pos_") + str(int(x_pos)) + str("_") + str(int(y_pos)) + str("_fit_figure.pdf")
	if not (save_fig=='display'):
		plt.savefig(fig_name_cust)
		bf.print_cust(f'{fig_name_cust} saved', quiet_val=quite_val)
	else:
		plt.show()
	plt.close('all')
##################################MAKE_FIGURE##################################




####################GET_SNR_MAP_FOR_BINNED_DATA####################
def get_snr_map_revised(file_name_rev_linear, file_name_rev_binned, par_dict, physical_axes = False, snr_vmin_val = 0.1, snr_vmax_val = 20, quiet_val=False, return_map=False, data_type_requested='snr_map'):
	bf.print_cust('Plotting SNR Map...', quiet_val=quiet_val)
	assert ('voronoi_snr_type' in par_dict), f"voronoi_snr_type required..."
	y_array, x_array, ra_array, dec_array, signal1, noise1, signal2, noise2 = np.loadtxt(file_name_rev_linear).T
	x_array_binned, y_array_binned, bin_num = np.loadtxt(file_name_rev_binned).T
	bin_unique = np.unique(bin_num)
	signal_total = np.zeros_like(bin_num)
	signal_total_med = np.zeros_like(bin_num)
	noise_total = np.ones_like(bin_num)
	if 'ew' in par_dict['voronoi_snr_type']:
		signal = signal1
		noise = noise1
	else:
		signal = signal2
		noise = noise2
	snr = signal / noise
	for i in range(len(bin_unique)):
		mask = np.isin(bin_num, bin_unique[i])
		signal_tmp = np.nansum(signal[mask])
		signal_median_tmp = np.nanmedian(signal[mask])
		noise_tmp = np.sqrt(np.nansum(noise[mask]**2))
		signal_total[mask] = signal_tmp
		signal_total_med[mask] = signal_median_tmp
		noise_total[mask] = noise_tmp
		signal_tmp = 0.0
		noise_tmp = 1.0
		
	if (data_type_requested=='median'):
		snr_revised = signal_total_med
	else:
		snr_revised = signal_total / noise_total

	if (return_map):
		if (physical_axes):
			return (ra_array, dec_array, snr_revised)
		else:
			return (x_array_binned, y_array_binned, snr_revised)
	else:
		bf.print_cust('Plotting SNR Map...', quiet_val=quiet_val)
		fig2, axs2 = plt.subplots(2, figsize=(6,10), sharex=True, sharey=True)
		if (physical_axes):
			im2 = axs2[0].scatter(dec, ra, c = snr, cmap='viridis', vmin=snr_vmin_val, vmax=snr_vmax_val)
		else:
			im2 = axs2[0].scatter(x_array, y_array, c = snr, cmap='viridis', vmin=snr_vmin_val, vmax=snr_vmax_val)
		axs2[0].set_title("SNR")
		add_colorbar(im2)
		#snr_revised = np.log10(snr_revised)
		im22 = axs2[1].scatter(x_array_binned, y_array_binned, c = snr_revised, cmap='viridis', vmin=snr_vmin_val, vmax=snr_vmax_val*100)
		#im2 = axs[1].scatter(x_array_binned, y_array_binned, c = np.log10(snr_revised), cmap='viridis')
		axs2[1].set_title("Binning")
		add_colorbar(im22)
		#plt.pause(1)
		plt.show()
		bf.print_cust('SNR Map Plotted...', quiet_val=quiet_val)
		#return None
		#continue
####################GET_SNR_MAP_FOR_BINNED_DATA####################





