from flask import render_template, abort
from flask.ext.login import current_user, login_required
from . import main
from forms import UploadForm
from run_ocr import run_ocr
from ..models import User, Pic
from .. import db


@main.route('/', methods=['GET', 'POST'])
def index():
    form = UploadForm()
    text = ''
    if form.validate_on_submit():
        if current_user.is_anonymous():
            #this is hacky, can work only because the id
            #for sqlite starts from 1
            user = User(id=0)
        else:
            user = User.query.filter_by(id=current_user.get_id()).first_or_404()
        text = run_ocr(user, form)
    return render_template('index.html', form=form, text=text)


@main.route('/user/<user_id>')
@login_required
def user(user_id):
    user = User.query.filter_by(id=int(user_id)).first()
    if user is None:
        abort(404)
    #pics = Pic.query.order_by(Pic.timestamp.desc()).all()
    pics = user.pics.order_by(Pic.timestamp.desc())
    return render_template('user/index.html', user=user, posts=pics)


#@main.route('/', defaults={'path':''})
#@main.route('/<path:path>')
#def catch_all(path):
#    return 'You want path: %s' % path
