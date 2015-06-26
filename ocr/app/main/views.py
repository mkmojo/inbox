from flask import render_template
from werkzeug import secure_filename
from . import main
from forms import UploadForm
import os

def run_ocr(photo_filename, username='anonymous', \
        ocr_engine='tesseract', language = 'eng'):
    photo_path = os.path.join(os.getcwd(), 'uploads/', photo_filename )
    out_path = os.path.join(os.getcwd(), 'output/', username, photo_filename)
    command = ' '.join([ocr_engine, '-l ' + language, photo_path, out_path])
    os.system(command)

    message = []
    with open(out_path+'.txt', 'r') as f:
        for line in f:
            line = line.decode('utf-8').strip()
            if line != '':
                message += [line]

    return message

@main.route('/', methods=['GET', 'POST'])
def index():
    form = UploadForm()
    text = ''
    if form.validate_on_submit():
        filename = secure_filename(form.photo.data.filename)
        form.photo.data.save('uploads/' + filename)

        username = 'test_user' # should be changed to get from database
        text = run_ocr(filename, username)
    return render_template('index.html', form=form, text=text)
