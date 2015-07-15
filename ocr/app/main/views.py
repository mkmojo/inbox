from flask import render_template, abort
from flask.ext.login import current_user, login_required
from . import main
from forms import UploadForm
from ..models import User, Pic
from .. import db


@main.route('/', methods=['GET', 'POST'])
def index():
    return render_template('main/index.html')


#@main.route('/', defaults={'path':''})
#@main.route('/<path:path>')
#def catch_all(path):
#    return 'You want path: %s' % path
