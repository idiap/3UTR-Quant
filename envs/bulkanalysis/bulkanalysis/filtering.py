from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
import scipy
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from genomemanager.overlaps import get_overlap_file


def remove_genes_not_expressed_in_all_replicates(df_counts_genes):
    # remove genes if not expressed in all replicates
    reps = df_counts_genes.columns
    for rep in reps:
        df_counts_genes = df_counts_genes[(df_counts_genes[rep] > 0)]

    print(f"After removing genes that are not expressed in all replicates, we have {str(len(df_counts_genes))} genes.")
    return df_counts_genes

def gmm_threshold(data, 
                  quantile_threshold=0.9,
                    background_quantile=True, 
                    intersection=True, 
                    density=False, 
                    bins=50, 
                    baysian_gmm=True, 
                    scaling_to_max=True, 
                    show_plot=True, 
                    saving_folder=None, 
                    data_name="", 
                    title="", 
                    extension='pdf', 
                    return_int=True):
    """
    Fit two gaussian distributions to your data and calculated a treshold from the quantile of one of the two gaussians.

    Args:
        data (array): 1D numpy array from which you want to fit two guassian distributions 
        quantile_threshold (float, optional): Quantile threshold you want to apply on the foreground or background gaussian. Defaults to 0.9.
        background_quantile (bool, optional): Apply the quantile threshold on the background (the gaussian with the lowest mean). Defaults to True.
        intersection (bool, optional): To plot or not the intesersection between the two gaussians or not. Defaults to True.
        density (bool, optional): To use density instead of count for the histogram. Defaults to False.
        bins (int, optional): Number of bins for the histogram. Defaults to 50.
        baysian_gmm (bool, optional): To use BayesianGaussianMixture instead of GaussianMixture. Defaults to True.
        scaling_to_max (bool, optional): To scale the gaussian to their max inside [mean-std; mean+std], if False it is scaled to the mean inside the same itervall. Defaults to True.
        show_plot (bool, optional): To show the plot or not. Defaults to True.
        saving_folder (str, optional): Saving folder to save the figure. Defaults to None.
        data_name (str, optional): Name of your data, to label x axis. Defaults to "".
        title (str, optional): Title of your plot. Defaults to "".

    Returns:
        float: data threshold calculated from the guassian quantile
    """       
    if baysian_gmm:
        gmm = BayesianGaussianMixture(n_components=2, random_state=0, max_iter=1000)
    else:
        gmm = GaussianMixture(n_components=2, random_state=0, max_iter=500)

    gmm.fit(np.reshape(data, (data.shape[0],1)))

    fig, ax = plt.subplots(figsize=(6,5))

    y_hist, x_hist,_ = plt.hist(data, density=density, alpha=0.5, color="gray", bins=bins)
    x_hist=(x_hist[1:]+x_hist[:-1])/2
    hist_df = pd.DataFrame(np.array([y_hist, x_hist]).T, columns=["y", "x"]) 

    x = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)

    if gmm.means_[0][0] > gmm.means_[1][0]:
        ind_list = [1,0]
    else:
        ind_list = [0,1]

    pdfs = {}
    for i, ind in enumerate(ind_list):
        mean = gmm.means_[ind][0]
        covariance = gmm.covariances_[ind][0,0]
        dist = scipy.stats.norm(loc=mean, scale=np.sqrt(covariance))
        pdf = dist.pdf(x)

        if (background_quantile and i == 0) or (not background_quantile and i==1):
            threshold = dist.ppf(quantile_threshold)
            color="orange"
        else:
            color="green"

        gauss_max = np.max(pdf)
        hist_max = hist_df[np.logical_and(hist_df["x"]>mean-np.sqrt(covariance), hist_df["x"]<mean+np.sqrt(covariance))]["y"]
        if scaling_to_max:
            hist_max = hist_max.max()
        else:
            hist_max = np.nanmean(hist_max.apply(lambda y: float("nan") if y==0 else y))
        scaling_factor = hist_max / gauss_max

        plt.plot(x, pdf*scaling_factor, label=f'Component {i + 1}', color=color)
        pdfs[f"Component {i + 1}"] = np.reshape(pdf*scaling_factor, -1)

    plt.axvline(threshold, color="red", label=f'{quantile_threshold*100:.0f}th: {threshold:.3f}')
    print("Threshold based on the {}th quantile: {}".format(quantile_threshold*100, np.round(threshold, 3)))

    if intersection:
        intersection_ind = np.argwhere(np.diff(np.sign(pdfs["Component 1"] - pdfs["Component 2"]))).flatten()
        intersection_point = x[intersection_ind].flatten()
        if show_plot:

            plt.scatter(x[intersection_ind], pdfs["Component 1"][intersection_ind], c="red", label="Intersection: {}".format(np.round(intersection_point.max(),3)), marker="X")
        print("Intersection x point: {}".format(intersection_point[0]))

    plt.legend()
    plt.xlabel(data_name)
    if density:
        plt.ylabel("Density")
    else:
        plt.ylabel("Count")

    plt.title(title)
    
    if saving_folder is not None:
        plt.savefig(os.path.join(saving_folder, f"gmm_threshold_{data_name}.{extension}"), bbox_inches='tight')
    if show_plot:
        plt.show()
    else:
        plt.gcf().set_visible(False)
    if return_int:
        return intersection_point
    else:
        return threshold

def filter_samples_based_on_thresholds(df_genes: pd.DataFrame, 
                                       treatments: dict, 
                                       thresholds: dict,
                                       pct_in_treatment: float = 0.75,
                                       save_file=None):
    """Filter the samples according to their respective thresholds. 
    A gene is kept if reliably expressed in at least *pct_in_treatment* percent 
    of samples of the same treatment.
    Args:
        df_genes: raw counts
        treatments: dictionnary with treatment names in keys and list of corresponding sample names in values.
        thresholds: sample names in keys and thresholds to filter the log-transformed counts in values.
        pct_in_treatment (float, optional): Proportion of samples that should express a gene within a treatment
        for the gene to be kept. Defaults to 0.75.
        save_file: Path to save the figure of log-transformed counts after filtering. Defaults to None.

    Returns:
        pd.DataFrame: The filtered raw gene counts.
    """
        
    # df_genes_log = np.log2(df_genes + 1)
    df_genes_log = df_genes[thresholds.keys()].apply(lambda x: np.log2(x+1))
    df_bool = df_genes_log.copy()

    # for col in df_genes_log.columns:
    for col in thresholds.keys():
        df_bool[col] = df_bool[col] >= thresholds[col]

    all_genes_to_keep = []
    for treatment, samples in treatments.items():
        samples = [sample for sample in samples if sample in df_genes_log.columns]
        nb_samples_required = int(np.ceil(pct_in_treatment*len(samples)))
        print(f"For treatment {treatment} ({len(samples)} samples): the gene should be reliably expressed in {nb_samples_required} samples.")
        df_bool_treatment = df_bool[samples]
        print(f"Number of reliably expressed genes for this condition: {len(list(set(list(df_bool_treatment[df_bool_treatment.sum(axis=1) >= nb_samples_required].index))))}")
        all_genes_to_keep += list(df_bool_treatment[df_bool_treatment.sum(axis=1) >= nb_samples_required].index)

        df_filtered = df_genes.loc[list(set(all_genes_to_keep))]
        df_filtered_log = df_genes_log.loc[list(set(all_genes_to_keep))]

        print(f"Number of genes reliably expressed at least in one of the condition: {len(list(set(all_genes_to_keep)))}")

        plt.figure()
        # sns.histplot(df_filtered_log[thresholds.keys()], stat='proportion', alpha=0.05, element='step')
        sns.histplot(df_filtered_log, stat='proportion', alpha=0.05, element='step')
        plt.legend('',frameon=False)
        plt.title("Gene expression after filtering", weight='bold')
        plt.xlabel("log2(read counts + 1)")

        if save_file is not None:
            plt.savefig(save_file, bbox_inches='tight')

        print("Filter samples based on thresholds: done.")
        print(f"Columns not filtered: {[col for col in df_genes.columns if col not in thresholds.keys()]}")

    return df_filtered

def remove_genes_with_overlapping_annotations_on_both_strand(df_counts_genes, gtf_file_quantification):
    overlaps = get_overlap_file(gtf_file_quantification)
    df_counts_genes = df_counts_genes[df_counts_genes.index.isin(overlaps.index) == False]
    print(f"After further removing genes with overlapping annotations, we have {str(len(df_counts_genes))} reliably expressed genes.")

    return df_counts_genes

def remove_outliers_samples(df_genes: pd.DataFrame, thresholds: dict, save_path: str = None, extension: str = 'pdf'):
    """Remove outliers samples in terms of thresholds and number of reliably expressed genes per samples.

    
    """
    def get_outliers(serie):
        """Get outliers of a ps.Series.
        Args:
            serie (pd.Series): pd.Serie with numerical values

        Returns:
            ps.Series: Two series with the lower outliers and upper outliers.
        """
        Q1, Q3 = serie.quantile([0.25, 0.75]).values
        lower_bound = Q1 - 1.5*(Q3-Q1)
        upper_bound = Q3 + 1.5*(Q3-Q1)
        return serie[serie < lower_bound], serie[serie > upper_bound]
    
    #Boxplot of number of cutting thresholds between foreground and background

    fig, ax = plt.subplots(1, 1, figsize=(2, 5))
    sns.boxplot(list(thresholds.values())) 
    plt.ylabel("Threshold values")   
    plt.xlabel("All samples")
    plt.title("Distribution of the 95th quantile \n of the background for all samples", weight='bold')
    if save_path is not None:
        plt.savefig(os.path.join(save_path, f"boxplot_thresholds_dsistributions.{extension}"), bbox_inches='tight')


    l, u = get_outliers(pd.Series(thresholds))

    samples_to_drop = list(l.index) + list(u.index)

    df_genes.drop(samples_to_drop, axis=1, inplace=True)
    df_genes_log = np.log2(df_genes+1)


    #Boxplot of number of reliably expressed genes per sample
    nb_reliably_expressed = {}
    for col in df_genes_log.columns:
        nb_reliably_expressed[col] = len(df_genes_log[df_genes_log[col] > thresholds[col]])


    fig, ax = plt.subplots(1, 1, figsize=(2, 5))
    sns.boxplot(list(nb_reliably_expressed.values())) 
    plt.ylabel("Number of reliably expressed genes")   
    plt.xlabel("All samples")
    plt.title("Distribution of the number of \n reliably expressed genes for all samples", weight='bold')
    if save_path is not None:
        plt.savefig(os.path.join(save_path, f"boxplot_nb_reliably_expressed_genes_dsistributions.{extension}"), bbox_inches='tight')

    l, u = get_outliers(pd.Series(nb_reliably_expressed))

    samples_to_drop = list(l.index) + list(u.index)

    df_genes.drop(samples_to_drop, axis=1, inplace=True)

    return df_genes