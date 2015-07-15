from flask import render_template, abort
from flask.ext.login import current_user, login_required
from ..models import User, Pic
from . import panel

@panel.route('/<user_id>')
@login_required
def user(user_id):
    user = User.query.filter_by(id=int(user_id)).first()
    if user is None or (current_user.id != user.id):
        abort(404)
    pics = user.pics.order_by(Pic.timestamp.desc())
    return render_template("panel/index.html", user=user, posts=pics)

@panel.route('/receipts')
@login_required
def receipts():
    return render_template("panel/receipts.html")

@panel.route('/transactionss')
@login_required
def transactions():
    return render_template("panel/transactions.html")
