from flask import render_template, abort
from flask.ext.login import current_user, login_required
from werkzeug import secure_filename
from . import main
from forms import UploadForm
from time import gmtime, strftime
from run_ocr import run_ocr
from ..models import User, Pic
from .. import db
import os


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

        filename = strftime("%Y-%m-%d-%H:%M:%S_", gmtime()) + \
                secure_filename(form.photo.data.filename)
        photo_path = os.path.join(os.getcwd(), 'uploads/', str(user.id),  filename)
        text_path = os.path.join(os.getcwd(), 'output/', str(user.id), filename)
        photo_dir_name = os.path.dirname(photo_path)
        try:
            if not os.path.exists(photo_dir_name):
                os.mkdir(photo_dir_name)
            form.photo.data.save(photo_path)
        except:
            print ("failed to save picture to: " + photo_path)

        text = run_ocr(user, photo_path, text_path)
    return render_template('index.html', form=form, text=text)


@main.route('/user/<user_id>')
@login_required
def user(user_id):
    user = User.query.filter_by(id=int(user_id)).first()
    if user is None:
        abort(404)
    pics = Pic.query.order_by(Pic.timestamp.desc()).all()
    return render_template('user/user.html', user=user, posts=pics)
