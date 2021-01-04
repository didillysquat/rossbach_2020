#!/usr/bin/env python3
from skbio.stats.ordination import pcoa
import matplotlib as mpl
mpl.use('TKAgg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import itertools
from skbio.stats.distance import permanova, DistanceMatrix, permdisp
import ecopy as ep
from Bio import SeqIO
from collections import Counter
from matplotlib.patches import Circle
from collections import defaultdict

DEGREE_SIGN = u'\N{DEGREE SIGN}'

class Susann:
    def __init__(self):

        self.fig, self.ax = plt.subplots(1,3, figsize=(16,5))

        path_to_seq_rel_seq_abund = '109_20200609_2020-06-09_11-59-06.357078.seqs.relative.abund_only.txt'
        path_to_meta = '109_20200609_2020-06-09_11-59-06.357078.seqs.absolute.meta_only.txt'
        path_to_profile = '109_20200609_2020-06-09_11-59-06.357078.profiles.absolute.abund_and_meta.txt'
        profile_df = pd.read_table(path_to_profile, skiprows=6, index_col=1)
        profile_df = profile_df.iloc[:-2].astype(float).astype(int)
        a1_samples = list(profile_df[profile_df['A1'] != 0].index)

        self.dist_df, to_drop, self.id_to_name_dict = self._make_dist_df()

        # make an abundance dataframe that is sample by sequence (Symbiodinium sequences only)
        # for running the simper analysis
        self.abund_df = pd.read_csv(path_to_seq_rel_seq_abund, sep='\t', index_col=0)
        self.abund_df = self.abund_df.loc[:,
                        [_ for _ in list(self.abund_df) if (_.startswith('A') or _.endswith('_A'))]]
        # Drop the negative samples etc that were dropped from the dist matrix.
        self.abund_df = self.abund_df.loc[[uid for uid in self.abund_df.index if str(uid) in self.id_to_name_dict.keys()],:]
        self.abund_df.index = [self.id_to_name_dict[str(_)] for _ in self.abund_df.index.values]
        # The group names
        group_names = [str(_).split('_')[0] for _ in self.abund_df.index]
        simper_result = ep.simper(data=self.abund_df, factor=group_names)
        simper_result.to_csv('simper_output.csv', sep=',', header=True, index=True)
        self._make_meta_df(path_to_meta, to_drop)
        self.unique_loc = self.meta_df['loc'].unique()
        self.site_names = [
            f'29{DEGREE_SIGN} Gulf of Aqaba', f'27{DEGREE_SIGN} Duba', f'22{DEGREE_SIGN} Thuwal', f'20{DEGREE_SIGN} Al Lith', f'18{DEGREE_SIGN} Farasan Banks'
        ]
        self.loc_to_site_dict = {loc: site_name for loc, site_name in zip(self.unique_loc, self.site_names)}

        # compute the permanova
        print('running permdisp\n\n')
        print(permdisp(distance_matrix=DistanceMatrix(self.dist_df),
                        grouping=[_.split('_')[0] for _ in list(self.dist_df)], permutations=999))
        print('running permanova\n\n')
        print(permanova(distance_matrix=DistanceMatrix(self.dist_df), grouping=[_.split('_')[0] for _ in list(self.dist_df)], permutations=9999))

        # now run the pcoa
        self._plot_pcoa_plots(a1_samples)
        # I want to add an additional figure here which will be better for quantifying
        # the distances between and within the sites. For each site I want to have
        # one bar that is coloured the same colour of the site, one bar for each of the other
        # sites, and then finally a grey bar.
        # The bar that is the colour as the site in question will show the average within site
        # distance. For each of the other sites, it will give the average between sites
        # Finally it will give the average betweeen all other sites as a grey bar.

        # we can put the results into a df where each row is a site and the columns
        # are the sites plus one extra for the overall average distance
        cols = list(self.unique_loc)
        cols.append('between_all')
        self.result_df = pd.DataFrame(index=self.unique_loc, columns=cols)
        self.stdev_df = pd.DataFrame(index=self.unique_loc, columns=cols)
        self._calc_pw_distances()
        # for mat df
        self._plot_heat_map(ax=self.ax[0], loc_to_site_dict=self.loc_to_site_dict, result_df=self.result_df)
        # self._plot_heat_map(ax=self.ax[3], loc_to_site_dict=self.loc_to_site_dict, result_df=self.stdev_df)
        # result_df.plot.bar(ax=ax[2], yerr=stdev_df)

        plt.tight_layout()
        plt.savefig(f's_fig_braycurtis_sqrt.png', dpi=600)
        plt.savefig(f's_fig_braycurtis_sqrt.svg', dpi=600)

    def _calc_pw_distances(self):
        for site_self in self.result_df.index:
            for site_compare in list(self.result_df):
                # Here we do the computations
                if site_self == site_compare:
                    # Then we want to be doing the within site average distance
                    # Get a list of the sites
                    self._within_site_av_distance(site_compare, site_self)
                elif site_compare == 'between_all':
                    # then we want to be doing the overall distance to other sites
                    self._overall_dist_to_other_sites(site_compare, site_self)
                else:
                    # Then we want to be doing the specific between site average distance
                    self._between_site_average_dist(site_compare, site_self)

    def _within_site_av_distance(self, site_compare, site_self):
        samples = self.meta_df[self.meta_df['loc'] == site_self].index
        tot_dist = [self.dist_df.at[site_1, site_2] for (site_1, site_2) in itertools.combinations(samples, 2)]
        self.result_df.at[site_self, site_compare] = sum(tot_dist) / len(tot_dist)
        self.stdev_df.at[site_self, site_compare] = np.std(tot_dist)

    def _overall_dist_to_other_sites(self, site_compare, site_self):
        samples = self.meta_df[self.meta_df['loc'] == site_self].index
        other_samples = self.meta_df[self.meta_df['loc'] != site_self].index
        tot_dist = []
        for s in samples:
            for o_s in other_samples:
                tot_dist.append(self.dist_df.at[s, o_s])
        self.result_df.at[site_self, site_compare] = sum(tot_dist) / len(tot_dist)
        self.stdev_df.at[site_self, site_compare] = np.std(tot_dist)

    def _between_site_average_dist(self, site_compare, site_self):
        samples_in = self.meta_df[self.meta_df['loc'] == site_self].index
        samples_out = self.meta_df[self.meta_df['loc'] == site_compare].index
        tot_dist = []
        for s_in in samples_in:
            for s_out in samples_out:
                tot_dist.append(self.dist_df.at[s_in, s_out])
        self.result_df.at[site_self, site_compare] = sum(tot_dist) / len(tot_dist)
        self.stdev_df.at[site_self, site_compare] = np.std(tot_dist)

    def _plot_pcoa_plots(self, a1_samples):
        pcoa_df_raw = pcoa(self.dist_df)
        self.pcoa_df = pcoa_df_raw.samples
        self.pcoa_df.index = self.dist_df.index
        self.colours = []
        # We want to differentiate between the sites with symbols as well as colours
        # for colourblind.
        self.shape_dict = {loc: symbol for loc, symbol in zip(self.unique_loc, reversed(["o", "s", "+", "^", "1"]), )}
        cols = reversed([d['color'] for d in list(mpl.rcParams["axes.prop_cycle"][:len(self.site_names)])])
        self.col_dict = {loc: col for loc, col in zip(self.unique_loc, cols)}
        for loc in self.unique_loc:
            samples = self.meta_df[self.meta_df['loc'] == loc].index
            scat = self.ax[1].plot(self.pcoa_df.loc[samples, 'PC1'], self.pcoa_df.loc[samples, 'PC2'], linestyle='none',
                                   label=loc, marker=self.shape_dict[loc], fillstyle='none', c=self.col_dict[loc])
            self.colours.append(scat[0].get_color())
        self.ax[1].set_xlabel(f'PC1 {pcoa_df_raw.proportion_explained[0]:.2f}')
        self.ax[1].set_ylabel(f'PC2 {pcoa_df_raw.proportion_explained[1]:.2f}')
        # Plot an arrow next to the A1 samples
        a1_samples_pcoa_one = [_ for _ in self.pcoa_df.index if _ in a1_samples]
        scat = self.ax[1].plot(self.pcoa_df.loc[a1_samples_pcoa_one, 'PC1'] + 0.04,
                               self.pcoa_df.loc[a1_samples_pcoa_one, 'PC2'],
                               linestyle='none', marker='<', c='black')
        for loc, c in zip(self.unique_loc, self.colours):
            samples = self.meta_df[self.meta_df['loc'] == loc].index
            self.ax[2].plot(self.pcoa_df.loc[samples, 'PC1'], self.pcoa_df.loc[samples, 'PC3'], c=c,
                            label=self.loc_to_site_dict[loc],
                            linestyle='none', marker=self.shape_dict[loc], fillstyle='none')
        self.ax[2].set_xlabel(f'PC1 {pcoa_df_raw.proportion_explained[0]:.2f}')
        self.ax[2].set_ylabel(f'PC3 {pcoa_df_raw.proportion_explained[2]:.2f}')
        self.ax[2].legend(loc='upper right', fontsize='x-small')
        # Plot an arrow next to the A1 samples
        scat = self.ax[2].plot(self.pcoa_df.loc[a1_samples_pcoa_one, 'PC1'] + 0.04,
                               self.pcoa_df.loc[a1_samples_pcoa_one, 'PC3'],
                               linestyle='none', marker='<', c='black')
        self._set_lims(ax=self.ax[1])
        self._set_lims(ax=self.ax[2])
        self.ax[1].set_aspect('equal', 'box')
        self.ax[2].set_aspect('equal', 'box')

    def _make_dist_df(self):
        path_to_dist = '2020-06-09_11-59-06.357078_braycurtis_sample_distances_A_sqrt.dist'
        with open(path_to_dist, 'r') as f:
            d = [_.rstrip().split('\t') for _ in f]
        samples = [_[0] for _ in d]
        ids = [_[1] for _ in d]
        data = [_[2:] for _ in d]
        id_to_name_dict = {uid: name for uid, name in zip(ids, samples)}
        dist_df = pd.DataFrame(data, index=samples, columns=samples).astype(float)
        to_drop = ['S1_H2O', 'S3_H2O', 'extraction_neg', 'milliq_neg']
        dist_df = dist_df.drop(columns=to_drop, index=to_drop)
        id_to_name_dict = {uid: id_to_name_dict[uid] for uid, name in id_to_name_dict.items() if name in dist_df.index}
        return dist_df, to_drop, id_to_name_dict

    def _make_meta_df(self, path_to_meta, to_drop):
        self.meta_df = pd.read_csv(path_to_meta, sep='\t')
        self.meta_df = self.meta_df.iloc[:-3, ]
        self.meta_df.set_index('sample_name', drop=True, inplace=True)
        self.meta_df = self.meta_df.drop(index=to_drop)
        self.meta_df = self.meta_df.loc[:, ['collection_latitude', 'collection_longitude']]
        self.meta_df['loc'] = [f'{latitude},{longitude}' for longitude, latitude in
                               zip(self.meta_df['collection_longitude'], self.meta_df['collection_latitude'])]
        self.meta_df.sort_values(by='collection_latitude', inplace=True, ascending=False)

    def _plot_heat_map(self, ax, loc_to_site_dict, result_df):
        result_df = result_df.astype(float)
        result_df.columns = [loc_to_site_dict[_] if _ in loc_to_site_dict else 'between_all' for _ in list(result_df)]
        result_df.index = [loc_to_site_dict[_] if _ in loc_to_site_dict else 'between_all' for _ in
                           result_df.index.values]
        result_df_heat = result_df.iloc[:, :-1]
        result_df_heat = result_df_heat.where(np.tril(np.ones(result_df_heat.shape)).astype(np.bool))
        heat_map = ax.imshow(result_df_heat)
        ax.set_xticks(np.arange(len(list(result_df_heat))))
        ax.set_yticks(np.arange(len(result_df_heat.index)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(list(result_df_heat), rotation='vertical')
        ax.set_yticklabels(result_df_heat.index.values)
        plt.colorbar(mappable=heat_map, ax=ax)
        for i in range(len(list(result_df_heat))):
            for j in range(len(list(result_df_heat))):
                if result_df_heat.iat[i, j] > 0.4:
                    text = ax.text(j, i, f'{result_df_heat.iat[i, j]:.2f}\n+/- {self.stdev_df.iat[i,j]:.2f}',
                                      ha="center", va="center", color="black", fontsize='x-small')
                else:
                    text = ax.text(j, i, f'{result_df_heat.iat[i, j]:.2f}\n+/- {self.stdev_df.iat[i,j]:.2f}',
                                      ha="center", va="center", color="w", fontsize='x-small')

    def _set_lims(self, ax):
        # Get the longest side and then set the small side to be the same length
        x_len = ax.get_xlim()[1] - ax.get_xlim()[0]
        y_len = ax.get_ylim()[1] - ax.get_ylim()[0]
        if y_len > x_len:
            # Then the y is the longest and we should adust the x to be bigger
            x_mid = ax.get_xlim()[0] + (x_len / 2)
            x_max = x_mid + (y_len / 2)
            x_min = x_mid - (y_len / 2)
            ax.set_xlim(x_min, x_max)
            ax.set_aspect('equal', 'box')
            return
        else:
            # Then the y is the longest and we should adust the x to be bigger
            y_mid = ax.get_ylim()[0] + (y_len / 2)
            y_max = y_mid + (x_len / 2)
            y_min = y_mid - (x_len / 2)
            ax.set_ylim(y_min, y_max)
            ax.set_aspect('equal', 'box')
            return

class HaploPies:
    """
    Quick utility class for quantifying the haplotypes that were found in the Pappas and webber data
    """
    def __init__(self):
        self.path_to_msa = '/Users/benjaminhume/Documents/projects/susann_clams/mol_ecol/dryad_files/pappas_weber_alignment.fasta'
        self.fig, self.ax = plt.subplots(2,3, figsize=(20,20))
        with open(self.path_to_msa, "r") as handle:
            record_list = list(SeqIO.parse(handle, "fasta"))

        # One of the Webber sequences has an 'n' in the sequence right at the position of one of the
        # haplotype defining nucleotide positions of the pappas dataset.
        # For the sake of this analysis we will return this to consensus.
        for i in range(len(record_list)):
            if '69004' in record_list[i].name:
                foo = 'bar'
            record_list[i].seq._data = record_list[i].seq._data.replace('n', 'g')


        # NB there are many errors in the Webber dataset.
        # In the thesis the locations of the reefs are wrong in the table but right in the
        # M and M section.
        # The sequence accession names have the wrong abbreviations in them too for some of the samples
        # I have used the locations from the thesis table instead.
        # I create a dictionary here to do so.
        self.loc_dict = {
            "GU068991": 'Dahab', "GU069006": 'Dahab', "GU069005": "Dahab",
            "GU069004": "Dahab", "GU069003": "Dahab", "GU068995": "Ras Nasrani",
            "GU068984": "Ras Nasrani", "GU068994": "Ras Nasrani", "GU068993": "Ras Nasrani",
            "GU068992": "Ras Nasrani", "GU068990": "Hurghada", "GU068989": "Hurghada", "GU069002": "Hurghada",
            "GU069001": "Hurghada", "GU069000": "Hurghada", "GU068999": "El Qeseir", "GU068998": "El Qeseir",
            "GU068997": "El Qeseir", "GU068988": "El Qeseir", "GU068996": "El Qeseir"
        }

        seq_count = Counter([r.seq for r in record_list])
        hap_one, hap_two, hap_three = [tup for tup in sorted(seq_count.items(), key=lambda x: x[1], reverse=True)[:3]]
        self.g_dict = {'other': "white", str(hap_one[0]): "#E7E7E7", str(hap_two[0]): "#A5A5A5", str(hap_three[0]): "#5C5C5C"}
        # plot one is the circles to scale
        for i, hap in enumerate([hap_one, hap_two, hap_three]):
            circle = Circle((i,0), radius=hap[1]/250, color=self.g_dict[str(hap[0])])
            self.ax[0, 0].add_patch(circle)
        self.ax[0, 0].set_xlim((-0.75, 2.25))
        self.ax[0, 0].set_ylim((-1.5, 1.5))


        # Now plot up pie charts for the individual sites
        for ax, site in zip(
                [self.ax[0,1], self.ax[0,2], self.ax[1,0], self.ax[1,1], self.ax[1,2]],
                ['Dahab', 'Ras Nasrani', 'Hurghada', 'El Qeseir', 'Thuwal']):
            # collect the abundances of each haplotypes
            count_dd = defaultdict(int)
            for r in record_list:
                assession = r.name.split('.')[0]
                try:
                    if self.loc_dict[assession] == site:
                        if str(r.seq) == str(hap_one[0]):
                            count_dd['hap_one'] += 1
                        elif str(r.seq) == str(hap_two[0]):
                            count_dd['hap_two'] += 1
                        elif str(r.seq) == str(hap_three[0]):
                            count_dd['hap_three'] += 1
                        else:
                            count_dd['other'] += 1
                except KeyError:
                    # Then this is a Pappas sequence
                    assert(r.name.startswith('LT'))
                    if site == 'Thuwal':
                        if str(r.seq) == str(hap_one[0]):
                            count_dd['hap_one'] += 1
                        elif str(r.seq) == str(hap_two[0]):
                            count_dd['hap_two'] += 1
                        elif str(r.seq) == str(hap_three[0]):
                            count_dd['hap_three'] += 1
                        else:
                            count_dd['other'] += 1
            ax.pie(
                x=[count_dd[_] for _ in ['hap_one', 'hap_two', 'hap_three', 'other']],
                labels=None,
                colors=[self.g_dict[_] for _ in [str(hap_one[0]), str(hap_two[0]), str(hap_three[0]), 'other']]
            )
            ax.set_title(site)

        plt.savefig('haplotype.svg')
        foo = 'bar'

HaploPies()