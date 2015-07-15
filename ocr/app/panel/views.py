from flask import render_template, abort
from flask.ext.login import current_user, login_required
from ..models import User, Pic
from . import panel
from run_ocr import run_ocr
from forms import UploadForm

@panel.route('/<user_id>')
@login_required
def user(user_id):
    user = User.query.filter_by(id=int(user_id)).first()
    if user is None or (current_user.id != user.id):
        abort(404)
    pics = user.pics.order_by(Pic.timestamp.desc())
    return render_template("panel/index.html", user=user, posts=pics)

@panel.route('/receipts', methods=['GET', 'POST'])
@login_required
def receipts():
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
    return render_template("panel/receipts.html", form=form, text=text )

@panel.route('/transactionss')
@login_required
def transactions():
    return render_template("panel/transactions.html")
