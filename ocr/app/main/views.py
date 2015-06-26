from flask import render_template
from werkzeug import secure_filename
from . import main
from forms import UploadForm

@main.route('/', methods=['GET', 'POST'])
def index():
    form = UploadForm()
    if form.validate_on_submit():
        filename = secure_filename(form.photo.data.filename)
        form.photo.data.save('uploads/' + filename)
    else:
        filename = None
    return render_template('index.html', form=form)
