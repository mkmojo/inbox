from flask import render_template
from flask.ext.login import current_user
from werkzeug import secure_filename
from . import main
from forms import UploadForm
from ..models import User, Pic
import subprocess
import os
from .. import db


def run_ocr(user, photo_path, text_path,\
        ocr_engine='tesseract', language = 'eng'):
    """
    photo_path: input/user uploaded file location
    text_path: output/interpreted file location
    !!NOTICE: tesseract appends '.txt' to the end of input filename
    """
    #create directory if not there
    text_dir = os.path.dirname(text_path)
    if not os.path.exists(text_dir):
        os.mkdir(text_dir)

    command = ' '.join([ocr_engine, '-l ' + language, photo_path, text_path])
    try:
        subprocess.call(command)
    except :
        return (['Command Fail:\n' + command])

    #save statistics to database
    #need to add in hash later in the future for better security
    pic = Pic(pic_path = photo_path, text_path = text_path + '.txt')
    db.session.add(pic)
    db.session.commit()

    message = []
    with open(text_path + '.txt', 'r') as f:
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
        if current_user.is_anonymous():
            #this is hacky, can work only because the id
            #for sqlite starts from 1
            user = User(id=0)
        else:
            user = User.query.filter_by(id=current_user.get_id()).first_or_404()

        filename = secure_filename(form.photo.data.filename)
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

