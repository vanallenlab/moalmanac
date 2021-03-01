from config import COLNAMES
from config import CONFIG

import io
import base64
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams["font.family"] = "sans-serif"
plt.switch_backend('agg')


class Illustrator(object):
    tableau10 = {
        'blue': (78, 121, 167), 'orange': (242, 142, 43), 'red': (225, 87, 89),
        'cyan': (118, 183, 178), 'green': (89, 161, 79), 'yellow': (237, 201, 72),
        'purple': (176, 122, 161), 'pink': (225, 157, 167), 'brown': (156, 117, 95),
        'grey': (186, 176, 172)}

    for color in tableau10.keys():
        r, g, b = tableau10[color]
        tableau10[color] = (r / 255., g / 255., b / 255.)

    @staticmethod
    def save_fig(figure, patient_id, suffix):
        outname = '.'.join([patient_id, suffix])
        figure.savefig(outname, bbox_inches='tight')


class Signatures(Illustrator):
    colors = [Illustrator.tableau10['blue'], Illustrator.tableau10['orange'], Illustrator.tableau10['red'],
              Illustrator.tableau10['grey'], Illustrator.tableau10['cyan'], Illustrator.tableau10['pink']]

    @classmethod
    def plot_context(cls, data, title='', ylabel='', grid=False, hline=0.0):
        fig = plt.figure(figsize=(20, 4))

        ax = plt.subplot()
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        plt.tick_params(axis="both", which="both", bottom=False, top=False,
                        labelbottom=True, left=False, right=False, labelleft=True)

        plt.bar(range(0, 16), data[0:16], color=cls.colors[0])       # C > A
        plt.bar(range(16, 32), data[16:32], color=cls.colors[1])     # C > G
        plt.bar(range(32, 48), data[32:48], color=cls.colors[2])     # C > T
        plt.bar(range(48, 64), data[48:64], color=cls.colors[3])     # T > A
        plt.bar(range(64, 80), data[64:80], color=cls.colors[4])     # T > C
        plt.bar(range(80, 95), data[80:95], color=cls.colors[5])     # T > G

        if float(hline) != 0.0:
            plt.axhline(hline, linewidth=2, color='k', alpha=0.4, linestyle='dashed')

        label_height = ax.get_ylim()[1] * 1.10
        labels = ['C > A', 'C > G', 'C > T', 'T > A', 'T > C', 'T > G']
        for i in range(0, len(labels)):
            label_ = labels[i]
            plt.text(6.05 + 16 * i, label_height, label_,  fontsize=15,
                     bbox={'facecolor': cls.colors[i], 'alpha': 0.4, 'pad': 3})

        plt.xlim(0, 95)
        xlabels = data.index.str[0] + data.index.str[4] + data.index.str[6]
        plt.xticks(range(0, 95), xlabels, rotation=90, fontsize=10)

        plt.ylabel(ylabel,  fontsize=14)
        plt.title(title,  fontsize=16, y=1.30)
        plt.grid(grid)
        plt.tight_layout()

        return fig

    @staticmethod
    def create_title(suffix):
        base = 'Trinucleotide contexts of somatic mutations'
        return ', '.join([base, suffix])

    @staticmethod
    def create_ylabel(suffix):
        base = 'Mutation Type'
        return '\n'.join([base, suffix])

    @classmethod
    def generate_context_plot(cls, patient_id, contexts):
        title = cls.create_title('counts')
        ylabel = cls.create_ylabel('Counts')

        figure = cls.plot_context(contexts.astype(float), title, ylabel)
        cls.save_fig(figure, patient_id, 'sigs.tricontext.counts.png')

    @classmethod
    def generate_context_plot_normalized(cls, patient_id, contexts):
        title = cls.create_title('normalized')
        ylabel = cls.create_ylabel('Probability')

        data = contexts.astype(float) / contexts.astype(float).max()
        figure = cls.plot_context(data.astype(float), title, ylabel)
        Illustrator.save_fig(figure, patient_id, 'sigs.tricontext.normalized.png')


class ValidationOverlap(Illustrator):
    config_section = 'validation_sequencing'
    min_af = float(CONFIG[config_section]['min_af_for_annotation'])
    min_power = float(CONFIG[config_section]['min_power'])

    feature_type_section = 'feature_types'
    feature_type_mutation = CONFIG[feature_type_section]['mut']

    section = 'validation_sequencing'
    gene = COLNAMES[section]['gene']
    feature_type = COLNAMES[section]['feature_type']
    alt_type = COLNAMES[section]['alt_type']
    alt = COLNAMES[section]['alt']
    tumor_f = COLNAMES[section]['tumor_f']
    coverage = COLNAMES[section]['coverage']
    validation_tumor_f = COLNAMES[section]['validation_tumor_f']
    validation_coverage = COLNAMES[section]['validation_coverage']
    validation_detection_power = COLNAMES[section]['validation_detection_power']

    columns = [tumor_f, coverage,
               validation_tumor_f, validation_coverage, validation_detection_power,
               gene, feature_type, alt]

    @classmethod
    def create_gene_alt_string(cls, data):
        return data[cls.gene].astype(str) + ' ' + data[cls.alt].astype(str)

    @classmethod
    def format_data(cls, df):
        idx = df[df[cls.feature_type] == cls.feature_type_mutation].index
        data = df.loc[idx, cls.columns].fillna(0.0)
        data[cls.coverage] = pd.to_numeric(data[cls.coverage])
        data[cls.tumor_f] = pd.to_numeric(data[cls.tumor_f])
        data[cls.validation_coverage] = pd.to_numeric(data[cls.validation_coverage])
        data[cls.validation_tumor_f] = pd.to_numeric(data[cls.validation_tumor_f])
        data[cls.validation_detection_power] = pd.to_numeric(data[cls.validation_detection_power])
        return data

    @classmethod
    def plot_overlap_af(cls, data, title=''):
        fig = plt.figure(figsize=(7, 7))

        ax = plt.subplot()
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        plt.tick_params(axis="both", which="both",
                        bottom=False, top=False, left=False, right=False,
                        labelbottom=True, labeltop=False, labelleft=True, labelright=False)

        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)

        powered = data[data[cls.validation_detection_power].astype(float) >= cls.min_power]
        lowpower = data[data[cls.validation_detection_power].astype(float) < cls.min_power]
        labels = cls.create_gene_alt_string(data)

        for idx in data.index:
            primary_tumor_f = float(data.loc[idx, cls.tumor_f])
            validation_tumor_f = float(data.loc[idx, cls.validation_tumor_f])
            if (primary_tumor_f >= cls.min_af) & (validation_tumor_f >= cls.min_af):
                ax.annotate(labels[idx], (primary_tumor_f, validation_tumor_f))

        plt.scatter(powered[cls.tumor_f].tolist(), powered[cls.validation_tumor_f].tolist(),
                    color=Illustrator.tableau10['blue'],
                    label=''.join(['Detection power in validation sequencing >= {}'.format(cls.min_power)]))
        plt.scatter(lowpower[cls.tumor_f].tolist(), lowpower[cls.validation_tumor_f].tolist(),
                    color=Illustrator.tableau10['grey'],
                    label=''.join(['Detection power in validation sequencing < {}'.format(cls.min_power)]))

        plt.xlim(-0.01, 1.01)
        plt.ylim(-0.01, 1.01)

        plt.plot(1, 1, ls="--", c=".3")
        plt.grid(False)
        plt.legend(fontsize=10, loc=2, frameon=False)

        plt.xlabel('Allelic fraction in primary',  fontsize=14)
        plt.ylabel('Allelic fraction in validation',  fontsize=14)
        plt.title(title,  fontsize=16)
        plt.tight_layout()

        return fig

    @classmethod
    def generate_dna_rna_plot(cls, df, patient_id):
        data = cls.format_data(df)
        figure = cls.plot_overlap_af(data, title=patient_id)
        Illustrator.save_fig(figure, patient_id, 'validation_overlap.png')


class PreclinicalEfficacy(Illustrator):
    section = 'preclinical'
    mut_values = COLNAMES[section]['mut_values']
    wt_values = COLNAMES[section]['wt_values']
    n_mut = COLNAMES[section]['n_mut_tested']
    n_wt = COLNAMES[section]['n_wt_tested']
    p_value = COLNAMES[section]['pvalue']

    @classmethod
    def convert_figure_base64(cls, fig):
        image = io.BytesIO()
        fig.savefig(image, tight=False, format='png')
        base64_image = base64.encodebytes(image.getvalue()).decode().replace('\n', '')
        return base64_image

    @classmethod
    def create_y_range(cls, values):
        minimum = min(values)
        maximum = max(values)

        yaxis_minimum_log10 = math.floor(math.log10(minimum)) - 0
        yaxis_maximum_log10 = math.ceil(math.log10(maximum)) + 1

        yticks = []
        for item in range(yaxis_minimum_log10, yaxis_maximum_log10):
            yticks.append(10**item)
        return yticks

    @classmethod
    def draw(cls, dictionary, drug, features, patient_id, feature_str):
        figure = cls.plot_boxplots(dictionary, drug, features)
        Illustrator.save_fig(figure, patient_id, '.'.join([feature_str, drug.split(' ')[0], 'png']))
        return cls.convert_figure_base64(figure)

    @classmethod
    def plot_boxplots(cls, dictionary, drug, features):
        def draw_boxplot(data, offset, width, axis, edge_color, alpha):
            for i in np.arange(1, len(data) + 1):
                x = np.random.normal(i, 0.04, size=len(data[i - 1]))
                plt.plot(x + offset, data[i - 1], color=edge_color, marker='.', alpha=alpha, linestyle='None')

            pos = np.arange(1, len(data) + 1) + offset
            bp = axis.boxplot(data, positions=pos, widths=width, manage_ticks=False, showfliers=False)
            for element in ['boxes', 'whiskers', 'caps']:
                plt.setp(bp[element], color=Illustrator.tableau10['grey'], linewidth=2)
            for element in ['medians']:
                plt.setp(bp[element], color=Illustrator.tableau10['orange'], linewidth=2)
            return bp

        def create_legend(color1, color2, label1, label2, font_size, location, frame_on):
            item1 = Line2D([0], [0], color='w', markerfacecolor=color1, label=label1, marker='o', markersize=10)
            item2 = Line2D([0], [0], color='w', markerfacecolor=color2, label=label2, marker='o', markersize=10)
            plt.legend(handles=[item1, item2], fontsize=font_size, loc=location, frameon=frame_on, ncol=2)

        def create_title(string):
            return 'Sensitivity of wild-type and mutant cell lines to\n{} (GDSC)'.format(string)

        wt_values = []
        mut_values = []
        n_wt = []
        n_mut = []
        p_values = []

        feature_length = len(features)
        for feature in features:
            feature_wt_values = dictionary[feature]['comparison'][cls.wt_values]
            feature_mut_values = dictionary[feature]['comparison'][cls.mut_values]
            if 0 in feature_wt_values:
                feature_wt_values.remove(0)
            if 0 in feature_mut_values:
                feature_mut_values.remove(0)

            wt_values.append(feature_wt_values)
            mut_values.append(feature_mut_values)
            n_wt.append(dictionary[feature]['comparison'][cls.n_wt])
            n_mut.append(dictionary[feature]['comparison'][cls.n_mut])
            p_values.append(dictionary[feature]['comparison'][cls.p_value])

        fig = plt.figure(figsize=(10, 7.5))
        ax = fig.add_subplot()
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        plt.tick_params(axis="both", which="both", bottom=False, top=False,
                        labelbottom=True, left=False, right=False, labelleft=True)

        bp_mut = draw_boxplot(mut_values, -0.15, 0.25, ax, Illustrator.tableau10['blue'], 0.50)
        bp_wt = draw_boxplot(wt_values, +0.15, 0.25, ax, Illustrator.tableau10['grey'], 0.40)

        plt.xlim([0.5, feature_length + 0.5])
        plt.yscale('log')
        ax.set_yticks(cls.create_y_range(wt_values[0] + mut_values[0]))
        ax.set_xticks(np.arange(1, len(features) + 1))
        ax.set_xticklabels([s.replace(' ', '\n') for s in features])

        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.ylabel('IC50 (micromolar)', fontsize=18)
        plt.title(create_title(drug),  fontsize=18)
        plt.tight_layout()

        create_legend(Illustrator.tableau10['blue'], Illustrator.tableau10['grey'], 'Mutant', 'Wild type', 16, 9, False)

        return fig
