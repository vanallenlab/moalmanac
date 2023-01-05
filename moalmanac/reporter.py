import flask
import flask_frozen
import datetime
import os
import tinydb

from config import COLNAMES
from config import CONFIG


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
    def drop_double_fusion(cls, dataframe):
        feature_type = COLNAMES[cls.report_section]['feature_type']
        alt = COLNAMES[cls.report_section]['alteration']
        rearrangement = CONFIG['feature_types']['fusion']

        idx_rearrangement = dataframe[dataframe[feature_type].eq(rearrangement)].index
        idx_rearrangement_keep = dataframe.loc[idx_rearrangement, :].drop_duplicates([alt], keep='first').index
        idx_rearrangement_drop = idx_rearrangement.difference(idx_rearrangement_keep)
        idx_keep = dataframe.index.difference(idx_rearrangement_drop)
        return dataframe.loc[idx_keep, :]

    @classmethod
    def generate_actionability_report(cls, actionable, report_dictionary, output_directory=""):
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

        versions = Reporter.generate_version_dictionary()
        report.add_versions(
            software=versions['software'],
            database=versions['database']
        )

        actionable = Reporter.drop_double_fusion(actionable)
        report.add_alterations(actionable)

        app = flask.Flask(__name__, static_folder=None)

        @app.route(f"/{report.metadata['patient_id']}.report.html")
        def index():
            return flask.render_template('index.html', report=report)

        freezer = flask_frozen.Freezer(app)
        app.config['FREEZER_DESTINATION'] = f"{os.getcwd()}" if output_directory == "" else output_directory
        app.config['FREEZER_REMOVE_EXTRA_FILES'] = False  # DO NOT REMOVE THIS, FLASK FROZEN WILL DELETE FILES IF TRUE

        @freezer.register_generator
        def index_generator():
            yield flask.url_for('index', report=report)

        freezer.freeze()

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

    @staticmethod
    def load_almanac_additional_matches(output_folder, patient_id):
        if output_folder == "":
            handle = f"{patient_id}.{CONFIG['databases']['additional_matches']}"
        else:
            handle = f"{output_folder}/{patient_id}.{CONFIG['databases']['additional_matches']}"
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


class ActionabilityReport:
    def __init__(self):
        self.metadata = {}
        self.versions = {}
        self.alterations = None
        self.additional_matches = None

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

    def add_alterations(self, alterations, additional_matches=None):
        self.alterations = alterations.fillna('')
        self.additional_matches = additional_matches
