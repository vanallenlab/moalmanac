import flask
import flask_frozen
import datetime
import tinydb

from config import COLNAMES
from config import CONFIG


class Reporter(object):
    report_section = 'report'
    code = COLNAMES[report_section]['code']
    date = COLNAMES[report_section]['date']
    description = COLNAMES[report_section]['description']
    normal = COLNAMES[report_section]['normal']
    tumor = COLNAMES[report_section]['tumor']
    ontology = COLNAMES[report_section]['ontology']
    patient_id = COLNAMES[report_section]['patient_id']
    stage = COLNAMES[report_section]['stage']
    purity = COLNAMES[report_section]['purity']
    ploidy = COLNAMES[report_section]['ploidy']
    ms_status = COLNAMES[report_section]['ms_status']

    @classmethod
    def drop_double_fusion(cls, dataframe):
        feature_type = COLNAMES[cls.report_section]['feature_type']
        alt = COLNAMES[cls.report_section]['alteration']
        rearrangement = CONFIG['feature_types']['fusion']

        idx_rearrangement = dataframe[dataframe[feature_type].eq(rearrangement)].index
        idx_rearrangement_keep = dataframe.loc[idx_rearrangement, :].drop_duplicates([alt], keep='first').index
        idx_rearrangement_drop = idx_rearrangement.difference(idx_rearrangement_keep)
        idx_keep = dataframe.index.difference(idx_rearrangement_drop)
        return dataframe.loc[idx_keep, :]

    @staticmethod
    def generate_date():
        return datetime.date.today().strftime("%b %d %Y")

    @classmethod
    def generate_version_dictionary(cls):
        version_section = 'versions'
        software_version = CONFIG[version_section]['interpreter']
        database_version = CONFIG[version_section]['database']
        return {
            'software': software_version,
            'database': database_version
        }

    @classmethod
    def generate_dictionary(cls, variants, patient):
        normal_barcode = cls.return_sample_barcode(variants, cls.normal)
        tumor_barcode = cls.return_sample_barcode(variants, cls.tumor)
        return {
            'patient_id': patient[cls.patient_id],
            'code': patient[cls.code].upper(),
            'ontology': patient[cls.ontology].title(),
            'normal_barcode': normal_barcode,
            'tumor_barcode': tumor_barcode,
            'stage': patient[cls.stage].title(),
            'description': patient[cls.description],
            'date': cls.generate_date(),
            'purity': patient[cls.purity],
            'ploidy': patient[cls.ploidy],
            'microsatellite_status': patient[cls.ms_status]
        }

    @classmethod
    def generate_report(cls, actionable, report_dictionary, version_dictionary,
                        preclinical_dictionary, preclinical_dataframe,
                        matchmaking_on, matchmaker, preclinical_reference_dict):
        app = flask.Flask(__name__)
        freezer = flask_frozen.Freezer(app)

        actionable = cls.drop_double_fusion(actionable)
        matches = cls.load_almanac_additional_matches()

        lookup = {}
        if matchmaking_on:
            matchmaker = matchmaker[~matchmaker['comparison'].eq('case-profile')]
            for index in matchmaker.index.tolist()[:10]:
                cell_line = matchmaker.loc[index, 'comparison']
                lookup[index] = preclinical_reference_dict[cell_line]

        @app.route('/')
        def index():
            return flask.render_template('index.html',
                                         df=actionable.fillna(''),
                                         dict=report_dictionary,
                                         version_dict=version_dictionary,
                                         matches=matches,
                                         preclinical_dict=preclinical_dictionary,
                                         preclinical_df=preclinical_dataframe,
                                         matchmaking_on=matchmaking_on,
                                         matchmaker=matchmaker,
                                         lookup=lookup
                                         )
        freezer.freeze()

    @classmethod
    def format_list_elements(cls, series):
        return (series
                .str.replace("\[\'", "")
                .str.replace("\'\]", "")
                .str.replace("\'", "")
                .str.split(',', expand=True)
                .apply(lambda x: ', '.join(x.sort_values().dropna().tolist()), axis=1)
                )

    @staticmethod
    def load_almanac_additional_matches():
        handle = CONFIG['databases']['additional_matches']
        db = tinydb.TinyDB(handle)
        matches = {}
        for table in db.tables():
            if table == '_default':
                continue
            matches[table] = db.table(table).all()
        return matches

    @classmethod
    def return_sample_barcode(cls, variants, column):
        if variants[column].any():
            return variants[column].unique()[0]
        else:
            return None
