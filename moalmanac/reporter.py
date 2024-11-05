import datetime
import flask
import flask_frozen
import logger
import pandas
import os

from config import COLNAMES


class Reporter:
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
    def drop_double_fusion(cls, dataframe, biomarker_type_string):
        feature_type = COLNAMES[cls.report_section]['feature_type']
        alt = COLNAMES[cls.report_section]['alteration']

        idx_rearrangement = dataframe[dataframe[feature_type].eq(biomarker_type_string)].index
        idx_rearrangement_keep = dataframe.loc[idx_rearrangement, :].drop_duplicates([alt], keep='first').index
        idx_rearrangement_drop = idx_rearrangement.difference(idx_rearrangement_keep)
        idx_keep = dataframe.index.difference(idx_rearrangement_drop)
        return dataframe.loc[idx_keep, :]

    @classmethod
    def format_alterations(cls, dataframe, config):
        if dataframe.empty:
            return dataframe

        dataframe = cls.drop_double_fusion(dataframe, biomarker_type_string=config['feature_types']['fusion'])

        lookup = COLNAMES['datasources']
        columns = [lookup['sensitivity'], lookup['resistance'], lookup['prognosis']]
        for column in columns:
            dataframe[column] = cls.format_clinical_columns(dataframe[column], convert_to_float=True)

        columns = [lookup['sensitive_implication'], lookup['resistance_implication'], lookup['prognostic_implication']]
        for column in columns:
            dataframe[column] = cls.format_clinical_columns(dataframe[column], convert_to_float=False)

        columns = [lookup['sensitivity_matches'], lookup['resistance_matches'], lookup['prognostic_matches']]
        for column in columns:
            if column in dataframe.columns:
                dataframe[column] = cls.preallocate_matches_columns(dataframe[column])

        if 'preclinical_efficacy_lookup' in dataframe.columns:
            dataframe.fillna({'preclinical_efficacy_lookup': ''}, inplace=True)

        return dataframe

    @staticmethod
    def format_clinical_columns(series, convert_to_float=False):
        series = (
            series
            .where(series.notna(), None)
            .where(~series.eq(''), None)
        )
        if convert_to_float:
            return series.astype(float)
        else:
            return series

    @classmethod
    def generate_actionability_report(cls, actionable, report_dictionary, config, similarity=None, output_directory=None):
        report = ActionabilityReport()
        report.add_metadata(
            name=report_dictionary['patient_id'],
            code=report_dictionary['code'],
            ontology=report_dictionary['ontology'],
            normal=report_dictionary['normal_barcode'],
            tumor=report_dictionary['tumor_barcode'],
            stage=report_dictionary['stage'],
            description=report_dictionary['description'],
            date=report_dictionary['date'],
            purity=report_dictionary['purity'],
            ploidy=report_dictionary['ploidy'],
            msi=report_dictionary['microsatellite_status']
        )

        versions = cls.generate_version_dictionary(config)
        report.add_versions(
            software=versions['software'],
            database=versions['database']
        )

        actionable = cls.format_alterations(dataframe=actionable, config=config)
        report.add_alterations(actionable)
        report.add_similar_profiles(similarity)

        app = flask.Flask(__name__, static_folder=None)

        output_path = f"/{report.metadata['patient_id']}.report.html"
        logger.Messages.general(message=f"Rendering report to: {output_path}")
        @app.route(output_path)
        def index():
            return flask.render_template('index.html', report=report)

        freezer = flask_frozen.Freezer(app)
        app.config['FREEZER_DESTINATION'] = f"{os.getcwd()}" if output_directory == ("" or None) else output_directory
        app.config['FREEZER_REMOVE_EXTRA_FILES'] = False  # DO NOT REMOVE THIS, FLASK FROZEN WILL DELETE FILES IF TRUE

        @freezer.register_generator
        def index_generator():
            yield flask.url_for('index', report=report)

        freezer.freeze()

    @staticmethod
    def generate_date():
        return datetime.date.today().strftime("%b %d %Y")

    @classmethod
    def generate_version_dictionary(cls, config):
        version_section = 'versions'
        software_version = config[version_section]['interpreter']
        database_version = config[version_section]['database']
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

    @staticmethod
    def preallocate_matches_columns(series):
        return series.apply(lambda x: [] if not isinstance(x, list) and pandas.isna(x) else x)

    @classmethod
    def return_sample_barcode(cls, variants, column):
        if variants[column].any():
            return variants[column].unique()[0]
        else:
            return None


class ActionabilityReport:
    def __init__(self):
        self.metadata = {}
        self.versions = {}
        self.alterations = None
        self.similar_profiles = None

    def add_metadata(self, name, code, ontology, normal, tumor, stage, description, date, purity, ploidy, msi):
        self.metadata['patient_id'] = name
        self.metadata['code'] = code
        self.metadata['ontology'] = ontology
        self.metadata['normal_barcode'] = normal
        self.metadata['tumor_barcode'] = tumor
        self.metadata['stage'] = stage
        self.metadata['description'] = description
        self.metadata['date'] = date
        self.metadata['purity'] = purity
        self.metadata['ploidy'] = ploidy
        self.metadata['microsatellite_status'] = msi

    def add_versions(self, software, database):
        self.versions['software'] = software
        self.versions['database'] = database

    def add_alterations(self, alterations):
        self.alterations = alterations

    def add_similar_profiles(self, similar_profiles):
        self.similar_profiles = similar_profiles
