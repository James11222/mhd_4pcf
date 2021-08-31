import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from nbodykit.lab import *
import mpl_toolkits.axes_grid1 as axgrid
from mpl_toolkits.axes_grid1 import make_axes_locatable

class plot_util_4pcf(object):
    
    def __init__(self, data = None, data_source=None):
        """
        This class allows us to measure the 4pcf from some input data field
        """
        if data is None:
            raise AssertionError("""You need to input a 6D array of 4PCF coefficients
            such that it can be indexed as: zeta[l1,l2,l3,b1,b2,b3]""")
        elif type(data) == list:
            self.zeta_list = data
        elif len(np.shape(data)) == 6:
            self.zeta = data
            
        if data_source is None:
            raise AssertionError("""You need to specify whether the data is from the FFT code, the Model, or Encore""")
        elif data_source == 'Model':
            self.data_source = data_source
        elif data_source == 'FFT':
            self.data_source = data_source
        elif data_source == 'Encore':
            self.data_source = data_source
        else:
            raise AssertionError("""Your data_source argument must be a string either 'Model', 'FFT', or 'Encore' """)
            
    def plot_2D(self, ells='000', b_1=0, b_2=0, b_3=0,xylabel_fontsize = 22,
    sub_title_fontsize = 22, sup_title_fontsize = 30,
    text_fontsize = 30, upper_bound = 5e-10,
    lower_bound = -1e-10, colorbar_ax_percent = '5%', nbins=4, ax=None):

        def make_plot_model(self):
            """
            This function is a wrapper function for plotting a 2D plot of the
            4PCF coefficients predicted by Jiamin's Model GRF code.
            """
            if ells in list(self.zeta.zetas_dict.keys()):
                if ax is not None:
        #         f, (ax1) = plt.subplots(1,1, figsize=(9,9))
                    zeta = self.zeta.zetas_dict[ells]
                    ells_string = '$\ell_1, \ell_2, \ell_3 =$ ' + ells[0] + ',' + ells[1] + ',' + ells[2]
        #             f.suptitle('4PCF coefficients $\\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$', fontsize=sup_title_fontsize)
                    im1 = ax.imshow(zeta[:,b_2,:].real, origin='lower')
                    ax.text(0,nbins-1, ells_string, c='white', fontsize=text_fontsize)
                    ax.set_title('$b_2$ = ' + str(b_2), fontsize=sub_title_fontsize)
                    ax.set_xticks(range(nbins))
                    ax.set_yticks(range(nbins))
                    ax.set_xlabel("$b_1$",fontsize=xylabel_fontsize)
                    ax.set_ylabel("$b_3$", fontsize=xylabel_fontsize)
                    divider1 = axgrid.make_axes_locatable(ax)
                    cax1 = divider1.append_axes("right", size=colorbar_ax_percent, pad=0.05)
                    plt.colorbar(im1, cax=cax1)
                    #im1.set_clim(lower_bound, upper_bound)
                else:
                    f, (ax1) = plt.subplots(1,1, figsize=(9,9))
                    zeta = self.zeta.zetas_dict[ells]
                    ells_string = '$\ell_1, \ell_2, \ell_3 =$ ' + ells[0] + ',' + ells[1] + ',' + ells[2]
                    f.suptitle('4PCF coefficients $\\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$', fontsize=sup_title_fontsize)

                    im1 = ax1.imshow(zeta[:,b_2,:].real, origin='lower')
                    ax1.text(0,nbins-1, ells_string, c='white', fontsize=text_fontsize)
                    ax1.set_title('$b_2$ = ' + str(b_2), fontsize=sub_title_fontsize)
                    ax1.set_xticks(range(nbins))
                    ax1.set_yticks(range(nbins))
                    ax1.set_xlabel("$b_1$",fontsize=xylabel_fontsize)
                    ax1.set_ylabel("$b_3$", fontsize=xylabel_fontsize)
                    divider1 = axgrid.make_axes_locatable(ax1)
                    cax1 = divider1.append_axes("right", size=colorbar_ax_percent, pad=0.05)
                    plt.colorbar(im1, cax=cax1)
                    #im1.set_clim(lower_bound, upper_bound)

            else:
                if ax is not None:
        #         f, (ax1) = plt.subplots(1,1, figsize=(9,9))
                    zeta = np.zeros((4,4,4))
                    ells_string = '$\ell_1, \ell_2, \ell_3 =$ ' + "Invalid Combination"
        #             f.suptitle('4PCF coefficients $\\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$', fontsize=sup_title_fontsize)
                    im1 = ax.imshow(zeta[:,b_2,:].real, origin='lower')
                    ax.text(0,nbins-1, ells_string, c='white', fontsize=text_fontsize)
                    ax.set_title('$b_2$ = ' + str(b_2), fontsize=sub_title_fontsize)
                    ax.set_xticks(range(nbins))
                    ax.set_yticks(range(nbins))
                    ax.set_xlabel("$b_1$",fontsize=xylabel_fontsize)
                    ax.set_ylabel("$b_3$", fontsize=xylabel_fontsize)
                    divider1 = axgrid.make_axes_locatable(ax)
                    cax1 = divider1.append_axes("right", size=colorbar_ax_percent, pad=0.05)
                    plt.colorbar(im1, cax=cax1)
                    #im1.set_clim(lower_bound, upper_bound)
                else:
                    f, (ax1) = plt.subplots(1,1, figsize=(9,9))
                    zeta = np.zeros((4,4,4))
                    ells_string = '$\ell_1, \ell_2, \ell_3 =$ ' + "Invalid Combination"
                    f.suptitle('4PCF coefficients $\\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$', fontsize=sup_title_fontsize)

                    im1 = ax1.imshow(zeta[:,b_2,:].real, origin='lower')
                    ax1.text(0,nbins-1, ells_string, c='white', fontsize=text_fontsize)
                    ax1.set_title('$b_2$ = ' + str(b_2), fontsize=sub_title_fontsize)
                    ax1.set_xticks(range(nbins))
                    ax1.set_yticks(range(nbins))
                    ax1.set_xlabel("$b_1$",fontsize=xylabel_fontsize)
                    ax1.set_ylabel("$b_3$", fontsize=xylabel_fontsize)
                    divider1 = axgrid.make_axes_locatable(ax1)
                    cax1 = divider1.append_axes("right", size=colorbar_ax_percent, pad=0.05)
                    plt.colorbar(im1, cax=cax1)
                    #im1.set_clim(lower_bound, upper_bound)

        def make_plot_FFT(self):
            """
            This function is a wrapper function for plotting a 2D plot of the
            4PCF coefficients predicted by Jiamin's Model GRF code.
            """

            zeta = self.zeta.zeta_normed_values
            nbins = self.zeta.nbins
            ell_max = self.zeta.ell_max                            
            l_1,l_2,l_3 = int(ells[0]), int(ells[1]), int(ells[2]) 
            ells_string = '$\ell_1, \ell_2, \ell_3 =$ ' + ells[0] + ',' + ells[1] + ',' + ells[2]
            if ax is None:
                f, (ax1) = plt.subplots(1,1, figsize=(9,9))

                f.suptitle('4PCF coefficients $\\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$', fontsize=sup_title_fontsize)

                im1 = ax1.imshow(zeta[l_1,l_2,l_3,:,b_2,:].real, origin='lower')
                ax1.text(0,nbins-1, ells_string, c='white', fontsize=text_fontsize)
                ax1.set_title('$b_2$ = ' + str(b_2), fontsize=sub_title_fontsize)
                ax1.set_xticks(range(nbins))
                ax1.set_yticks(range(nbins))
                ax1.set_xlabel("$b_1$",fontsize=xylabel_fontsize)
                ax1.set_ylabel("$b_3$", fontsize=xylabel_fontsize)
                divider1 = axgrid.make_axes_locatable(ax1)
                cax1 = divider1.append_axes("right", size=colorbar_ax_percent, pad=0.05)
                plt.colorbar(im1, cax=cax1)
            else:
                im1 = ax.imshow(zeta[l_1,l_2,l_3,:,b_2,:].real, origin='lower')
                ax.text(0,nbins-1, ells_string, c='white', fontsize=text_fontsize)
                ax.set_title('$b_2$ = ' + str(b_2), fontsize=sub_title_fontsize)
                ax.set_xticks(range(nbins))
                ax.set_yticks(range(nbins))
                ax.set_xlabel("$b_1$",fontsize=xylabel_fontsize)
                ax.set_ylabel("$b_3$", fontsize=xylabel_fontsize)
                divider1 = axgrid.make_axes_locatable(ax)
                cax1 = divider1.append_axes("right", size=colorbar_ax_percent, pad=0.05)
                plt.colorbar(im1, cax=cax1)

        def make_plot_encore(self):

            l_1,l_2,l_3 = int(ells[0]), int(ells[1]), int(ells[2]) 
            ells_string = '$\ell_1, \ell_2, \ell_3 =$ ' + ells[0] + ',' + ells[1] + ',' + ells[2]
            if ax is None:
                f, (ax1) = plt.subplots(1,1, figsize=(9,9))
                f.suptitle('4PCF coefficients $\\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$', fontsize=sup_title_fontsize)

                im1 = ax1.imshow(self.zeta[l_1,l_2,l_3,:,b_2,:].real, origin='lower')
                ax1.text(0,nbins-1, ells_string, c='white', fontsize=text_fontsize)
                ax1.set_title('$b_2$ = ' + str(b_2), fontsize=sub_title_fontsize)
                ax1.set_xticks(range(nbins))
                ax1.set_yticks(range(nbins))
                ax1.set_xlabel("$b_1$",fontsize=xylabel_fontsize)
                ax1.set_ylabel("$b_3$", fontsize=xylabel_fontsize)
                divider1 = axgrid.make_axes_locatable(ax1)
                cax1 = divider1.append_axes("right", size=colorbar_ax_percent, pad=0.05)
                plt.colorbar(im1, cax=cax1)
            else:
                im1 = ax.imshow(self.zeta[l_1,l_2,l_3,:,b_2,:].real, origin='lower')
                ax.text(0,nbins-1, ells_string, c='white', fontsize=text_fontsize)
                ax.set_title('$b_2$ = ' + str(b_2), fontsize=sub_title_fontsize)
                ax.set_xticks(range(nbins))
                ax.set_yticks(range(nbins))
                ax.set_xlabel("$b_1$",fontsize=xylabel_fontsize)
                ax.set_ylabel("$b_3$", fontsize=xylabel_fontsize)
                divider1 = axgrid.make_axes_locatable(ax)
                cax1 = divider1.append_axes("right", size=colorbar_ax_percent, pad=0.05)
                plt.colorbar(im1, cax=cax1)
        
        if self.data_source == 'Model':
            make_plot_model(self)
        elif self.data_source == 'FFT':
            make_plot_FFT(self)
        elif self.data_source == 'Encore':
            make_plot_encore(self)
        else:
            raise AssertionError("You need to give a data_source so I know how to plot")
            
        
        
            
            
            
            
#####################################################################################
##################             Sawtooth Plots             ###########################
#####################################################################################
            
    
    
    def plot_sawtooth(self, average_bins=None, ells='000', ax = None,
                        mean_color = 'limegreen', 
                        error_color = 'steelblue',
                        nbins = 4, bin_overlap=True, factor=None, alpha_ =0.7):
        if factor is not None:
            factor_flag = True
        else:
            factor_flag = False
            
        if average_bins is None:
            raise AssertionError("you need to give a list of the centers of all radial bins")
        
        
        def plot_sawtooth_model(self):

            """
            This function is a wrapper function for plotting a sawtooth plot of the
            4PCF coefficients predicted by Jiamin's Model GRF code.
            """

            
            N = len(self.zeta_list)
            l_1 = int(ells[0])
            l_2 = int(ells[1])
            l_3 = int(ells[2])
            if bin_overlap==True:
                zeta_ell_all = []
                for j in range(N):
                    zeta_ell_j = []
                    zeta_j = self.zeta_list[j].zetas_dict[ells]
                    bin_indexes = []
                    b_i = 1
                    for b1 in range(0,nbins):
                        for b2 in range(b1,nbins):
                            for b3 in range(b2,nbins):
                                value = (zeta_j[b1,b2,b3].real * average_bins[b1] * 
                                average_bins[b2] * average_bins[b3])
                                zeta_ell_j.append(value)
                                bin_indexes.append(b_i)
                                b_i += 1
                    zeta_ell_all.append(zeta_ell_j)

                zeta_ell_all_array = np.array(zeta_ell_all)
                mean_zeta = np.mean(zeta_ell_all_array, axis=0)
                error_zeta = np.std(zeta_ell_all_array, axis=0)
            else:
                zeta_ell_all = []
                for j in range(N):
                    zeta_ell_j = []
                    zeta_j = self.zeta_list[j].zetas_dict[ells]
                    bin_indexes = []
                    b_i = 1
                    for b1 in range(0,nbins):
                        for b2 in range(b1+1,nbins):
                            for b3 in range(b2+1,nbins):
                                value = (zeta_j[b1,b2,b3].real * average_bins[b1] * 
                                average_bins[b2] * average_bins[b3])
                                zeta_ell_j.append(value)
                                bin_indexes.append(b_i)
                                b_i += 1
                    zeta_ell_all.append(zeta_ell_j)

                zeta_ell_all_array = np.array(zeta_ell_all)
                mean_zeta = np.mean(zeta_ell_all_array, axis=0)
                error_zeta = np.std(zeta_ell_all_array, axis=0)


            if factor_flag == True:
                mean_zeta *= factor
                error_zeta *= factor
                
            if ax is not None:
                ax.axhline(y=0, xmin=0, xmax=64, color='white', linestyle='--')
                ax.plot(bin_indexes, mean_zeta, color=mean_color, linewidth=3, 
                label="Model $\\Lambda" + " = (" + str(l_1) + "," + str(l_2) + "," + str(l_3)+")$")
                ax.fill_between(bin_indexes, y1 = mean_zeta - error_zeta, 
                                y2 = mean_zeta + error_zeta, color=error_color, alpha = alpha_)
                ax.set_xlabel("Bin Index")
                ax.set_ylabel("$\\bar{r}_1 \\bar{r}_2 \\bar{r}_3 \\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$")
                # ax.legend(fontsize=16)
            else:
                f, (ax1) = plt.subplots(1,1, figsize=(16,10), sharex=True)
                ax1.axhline(y=0, xmin=0, xmax=64, color='white', linestyle='--')
                ax1.plot(bin_indexes, mean_zeta, color=mean_color, linewidth=3, 
                label="Model $\\Lambda" + " = (" + str(l_1) + "," + str(l_2) + "," + str(l_3)+")$")
                ax1.fill_between(bin_indexes, y1 = mean_zeta - error_zeta, 
                                y2 = mean_zeta + error_zeta, color=error_color, alpha = alpha_)
                ax1.set_xlabel("Bin Index")
                ax1.set_ylabel("$\\bar{r}_1 \\bar{r}_2 \\bar{r}_3 \\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$")
                ax1.legend(fontsize=22)
                f.tight_layout()

        def plot_sawtooth_FFT(self):

            """
            This function is a wrapper function for plotting a sawtooth plot of the
            4PCF coefficients measured by the FFT Code.
            """
            N = len(self.zeta_list)
            l_1 = int(ells[0])
            l_2 = int(ells[1])
            l_3 = int(ells[2])
            nbins = len(average_bins)
            if bin_overlap == True:
                zeta_ell_all = []
                for j in range(N):
                    zeta_ell_j = []
                    zeta_j = self.zeta_list[j].zeta_normed_values
                    bin_indexes = []
                    b_i = 1
                    for b1 in range(nbins):
                        for b2 in range(b1,nbins):
                            for b3 in range(b2,nbins):
                                value = (zeta_j[l_1,l_2,l_3,b1,b2,b3].real * average_bins[b1] * 
                                average_bins[b2] * average_bins[b3])
                                zeta_ell_j.append(value)
                                bin_indexes.append(b_i)
                                b_i += 1
                    zeta_ell_all.append(zeta_ell_j)

                zeta_ell_all_array = np.array(zeta_ell_all)

                mean_zeta = np.mean(zeta_ell_all_array, axis=0)
                error_zeta = np.std(zeta_ell_all_array, axis=0)

            else:
                zeta_ell_all = []
                for j in range(N):
                    zeta_ell_j = []
                    zeta_j = self.zeta_list[j].zeta_normed_values
                    bin_indexes = []
                    b_i = 1
                    for b1 in range(0,nbins):
                        for b2 in range(b1+1,nbins):
                            for b3 in range(b2+1,nbins):
                                value = (zeta_j[l_1,l_2,l_3,b1,b2,b3].real * average_bins[b1] * 
                                average_bins[b2] * average_bins[b3])
                                zeta_ell_j.append(value)
                                bin_indexes.append(b_i)
                                b_i += 1
                    zeta_ell_all.append(zeta_ell_j)

                zeta_ell_all_array = np.array(zeta_ell_all)

                mean_zeta = np.mean(zeta_ell_all_array, axis=0)
                error_zeta = np.std(zeta_ell_all_array, axis=0)
                
            if factor_flag == True:
                mean_zeta *= factor
                error_zeta *= factor


            if ax is not None:
                ax.axhline(y=0, xmin=0, xmax=64, color='white', linestyle='--')
                ax.plot(bin_indexes, mean_zeta, color=mean_color, linewidth=3, 
                label="FFT $\\Lambda" + " = (" + str(l_1) + "," + str(l_2) + "," + str(l_3)+")$")
                ax.fill_between(bin_indexes, y1 = mean_zeta - error_zeta/2, 
                                y2 = mean_zeta + error_zeta/2, color= error_color, alpha = alpha_)
                ax.set_xlabel("Bin Index")
                ax.set_ylabel("$\\bar{r}_1 \\bar{r}_2 \\bar{r}_3 \\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$")
                #ax.legend(fontsize=16)
            else:
                f, (ax1) = plt.subplots(1,1, figsize=(16,10), sharex=True)
                # ax1.set_title('Unscaled')
                ax1.axhline(y=0, xmin=0, xmax=64, color='white', linestyle='--')
                ax1.plot(bin_indexes, mean_zeta, color=mean_color, linewidth=3, 
                label="FFT $\\Lambda" + " = (" + str(l_1) + "," + str(l_2) + "," + str(l_3)+")$")
                ax1.fill_between(bin_indexes, y1 = mean_zeta - error_zeta/2, 
                                y2 = mean_zeta + error_zeta/2, color=error_color, alpha = alpha_)
                ax1.set_xlabel("Bin Index")
                ax1.set_ylabel("$\\bar{r}_1 \\bar{r}_2 \\bar{r}_3 \\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$")
                ax1.legend(fontsize=22)
                f.tight_layout()


        def plot_sawtooth_encore(self):

            """
            This function is a wrapper function for plotting a sawtooth plot of the
            4PCF coefficients measured by the FFT Code.
            """

            N = len(self.zeta_list)
            nbins = len(average_bins)
            l_1,l_2,l_3 = int(ells[0]),int(ells[1]),int(ells[2])
            zeta_ell_all = []
            for j in range(N):
                zeta_ell_j = []
                zeta_j = self.zeta_list[j]
                bin_indexes = []
                b_i = 1
                for b1 in range(0,nbins):
                    for b2 in range(b1+1,nbins):
                        for b3 in range(b2+1,nbins):
                            value = (zeta_j[l_1,l_2,l_3,b1,b2,b3].real * average_bins[b1] * 
                            average_bins[b2] * average_bins[b3])
                            zeta_ell_j.append(value)
                            bin_indexes.append(b_i)
                            b_i += 1
                zeta_ell_all.append(zeta_ell_j)

            zeta_ell_all_array = np.array(zeta_ell_all)

            mean_zeta = np.mean(zeta_ell_all_array, axis=0)
            error_zeta = np.std(zeta_ell_all_array, axis=0)
            
            if factor_flag == True:
                mean_zeta *= factor
                error_zeta *= factor
            
            if ax is not None:
                ax.axhline(y=0, xmin=0, xmax=64, color='white', linestyle='--')
                ax.plot(bin_indexes, mean_zeta, color=mean_color, linewidth=3, 
                label="Encore $\\Lambda" + " = (" + str(l_1) + "," + str(l_2) + "," + str(l_3)+")$")
                ax.fill_between(bin_indexes, y1 = mean_zeta - error_zeta, 
                                y2 = mean_zeta + error_zeta, color= error_color, alpha = alpha_)
                ax.set_xlabel("Bin Index")
                ax.set_ylabel("$\\bar{r}_1 \\bar{r}_2 \\bar{r}_3 \\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$")
                #ax.legend(fontsize=16)
            else:
                f, (ax1) = plt.subplots(1,1, figsize=(16,10), sharex=True)
                # ax1.set_title('Unscaled')
                ax1.axhline(y=0, xmin=0, xmax=64, color='white', linestyle='--')
                ax1.plot(bin_indexes, mean_zeta, color=mean_color, linewidth=3, 
                label="Encore $\\Lambda" + " = (" + str(l_1) + "," + str(l_2) + "," + str(l_3)+")$")
                ax1.fill_between(bin_indexes, y1 = mean_zeta - error_zeta, 
                                y2 = mean_zeta + error_zeta, color=error_color, alpha = alpha_)
                ax1.set_xlabel("Bin Index")
                ax1.set_ylabel("$\\bar{r}_1 \\bar{r}_2 \\bar{r}_3 \\zeta^{b_1 b_2 b_3}_{\ell_1 \ell_2 \ell_3}$")
                ax1.legend(fontsize=22)
                f.tight_layout()
            
        if self.data_source == 'Model':
            plot_sawtooth_model(self)
        elif self.data_source == 'FFT':
            plot_sawtooth_FFT(self)
        elif self.data_source == 'Encore':
            plot_sawtooth_encore(self)
        else:
            raise AssertionError("You need to give a data_source so I know how to plot")




            



