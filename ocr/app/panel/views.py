from flask import render_template, abort
from flask.ext.login import current_user, login_required
from . import main
from forms import UploadForm
from run_ocr import run_ocr
from ..models import User, Pic
from .. import db


@panel.route('/user/<user_id>')
@login_required
def user(user_id):
    user = User.query.filter_by(id=int(user_id)).first()
    if user is None:
        abort(404)
    pics = user.pics.order_by(Pic.timestamp.desc())
    return render_template('index.html', user=user, posts=pics)
