import os
basedir = os.path.abspath(os.path.dirname(__file__))


class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'hard to guess string'
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    MAIL_SERVER = 'smtp.googlemail.com'
    MAIL_PORT = 587
    MAIL_USE_TLS = True
    MAIL_USERNAME = os.environ.get('MAIL_USERNAME') or 'qqy.receipt@gmail.com'
    MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD') or 'QiyuanQiu211'
    MYBILL_MAIL_SUBJECT_PREFIX = '[MyBill]'
    MYBILL_MAIL_SENDER = 'Flasky Admin <flasky@example.com>'
    MYBILL_ADMIN = os.environ.get('MYBILL_ADMIN')
    MYBILL_POSTS_PER_PAGE = 20
    STATIC_FOLDER = os.environ.get('STATIC_FOLDER') or "/home/qqiu/workspace/inbox/ocr/app/static"

    @staticmethod
    def init_app(app):
        pass


class DevelopmentConfig(Config):
    DEBUG = True
    SQLALCHEMY_DATABASE_URI = os.environ.get('DEV_DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data-dev.sqlite')


class TestingConfig(Config):
    TESTING = True
    SQLALCHEMY_DATABASE_URI = os.environ.get('TEST_DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data-test.sqlite')


class ProductionConfig(Config):
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data.sqlite')


config = {
    'development': DevelopmentConfig,
    'testing': TestingConfig,
    'production': ProductionConfig,

    'default': DevelopmentConfig
}
